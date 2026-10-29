/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#include <stdint.h>

/*
   Thread safety of library (de)initialization
   -------------------------------------------

   FLA_Init_safe()/FLA_Finalize_safe() bracket every lapack2flame entry point
   (e.g. dgesvd_, dgesdd_). Previously they used an unsynchronized global
   flag and "whoever initialized, finalizes" ownership. When several threads
   call a LAPACK routine concurrently, the first thread to return would call
   FLA_Finalize() and free the global constants (FLA_ONE, ...) and control
   trees while the other threads were still using them, producing random
   aborts ("not yet implemented", "inconsistent datatypes", NULL pointers),
   heap corruption, or hangs.

   The state below fixes this:

   - fla_init_lock: a spinlock that serializes every transition of the
     library state. It is statically initialized, so it is usable before
     FLA_Init() and needs neither pthreads nor FLA_Lock (which requires a
     runtime init and only exists in multithreaded builds).

   - fla_safe_users: number of FLA_Init_safe() calls that have not yet been
     matched by FLA_Finalize_safe(), i.e. library calls in flight.

   - fla_auto_finalize: TRUE when the library must be finalized as soon as
     fla_safe_users drops to zero. It is set when FLA_Init_safe() had to
     initialize the library itself, cleared when the user calls FLA_Init()
     (the user now owns the lifetime), and set again when the user calls
     FLA_Finalize() while library calls are still in flight (finalization
     is deferred to the last in-flight call instead of pulling the state
     out from under it).

   Callers that only need the old behaviour are unaffected: a single-threaded
   program sees exactly the same init/finalize sequence as before.

   Lifetime of library state
   -------------------------

   Every lapack2flame call that finds the library uninitialized initializes
   it and, when it is the last call in flight, finalizes it again. To keep
   that cycle cheap, state is split by how it is stored:

   - Statically allocated state (.bss/.data/.rodata) is initialized once per
     process and is NOT torn down by FLA_Finalize(): the error message table,
     the memory leak counter (count, status and lock), and the global scalar
     constants FLA_ONE, FLA_ZERO, ..., FLA_OVERFLOW_SQUARE_THRES, whose base
     objects and buffers now live in static storage (see
     FLA_Init_constants()). They therefore remain valid after FLA_Finalize().

   - Dynamically allocated state is still created by FLA_Init() and freed by
     FLA_Finalize(): the control trees (FLA_Cntl_init/FLA_Cntl_finalize) and
     the SuperMatrix GPU/HIP device handles.

   Performance note: if no thread holds the library initialized, concurrent
   callers may still init/finalize the dynamic state repeatedly whenever the
   in-flight count touches zero. That is correct but costs time; calling
   FLA_Init() once at program start (FLA_INIT() from Fortran) avoids it.
*/

static FLA_Bool FLA_initialized   = FALSE;
static int      fla_safe_users    = 0;
static FLA_Bool fla_auto_finalize = FALSE;

// TRUE once the process-lifetime static state (error message table, memory
// leak counter) has been initialized. Only accessed with fla_init_lock held.
static FLA_Bool fla_static_state_ready = FALSE;

#if defined(_MSC_VER)

#include <intrin.h>

static volatile long fla_init_lock_word = 0;

static void fla_init_lock( void )
{
  // _InterlockedExchange is a full barrier on MSVC.
  while ( _InterlockedExchange( &fla_init_lock_word, 1 ) != 0 )
    while ( fla_init_lock_word != 0 ) ;
}

static void fla_init_unlock( void )
{
  _InterlockedExchange( &fla_init_lock_word, 0 );
}

static FLA_Bool fla_load_initialized( void )
{
  return *( volatile FLA_Bool* ) &FLA_initialized;
}

static void fla_store_initialized( FLA_Bool val )
{
  *( volatile FLA_Bool* ) &FLA_initialized = val;
}

#else

#include <sched.h>

static int fla_init_lock_word = 0;

static void fla_init_lock( void )
{
  // Test-and-test-and-set: only attempt the atomic exchange when the lock
  // looks free, and yield the CPU while it is held. Initialization is short
  // and rare, so contention is brief.
  while ( __atomic_exchange_n( &fla_init_lock_word, 1, __ATOMIC_ACQUIRE ) != 0 )
    while ( __atomic_load_n( &fla_init_lock_word, __ATOMIC_RELAXED ) != 0 )
      sched_yield();
}

static void fla_init_unlock( void )
{
  __atomic_store_n( &fla_init_lock_word, 0, __ATOMIC_RELEASE );
}

// FLA_initialized is read without the lock by FLA_Initialized(), so
// accesses go through atomics to keep them well defined.
static FLA_Bool fla_load_initialized( void )
{
  return __atomic_load_n( &FLA_initialized, __ATOMIC_ACQUIRE );
}

static void fla_store_initialized( FLA_Bool val )
{
  __atomic_store_n( &FLA_initialized, val, __ATOMIC_RELEASE );
}

#endif

FLA_Obj FLA_THREE = {};
FLA_Obj FLA_TWO = {};
FLA_Obj FLA_ONE = {};
FLA_Obj FLA_ONE_HALF = {};
FLA_Obj FLA_ZERO = {};
FLA_Obj FLA_MINUS_ONE_HALF = {};
FLA_Obj FLA_MINUS_ONE = {};
FLA_Obj FLA_MINUS_TWO = {};
FLA_Obj FLA_MINUS_THREE = {};

FLA_Obj FLA_EPSILON = {};
FLA_Obj FLA_SAFE_MIN = {};
FLA_Obj FLA_SAFE_MIN_SQUARE = {};
FLA_Obj FLA_SAFE_INV_MIN = {};
FLA_Obj FLA_SAFE_INV_MIN_SQUARE = {};
FLA_Obj FLA_UNDERFLOW_THRES = {};
FLA_Obj FLA_OVERFLOW_THRES = {};
FLA_Obj FLA_UNDERFLOW_SQUARE_THRES = {};
FLA_Obj FLA_OVERFLOW_SQUARE_THRES = {};

const float    fzero = 0.0f;
const double   dzero = 0.0;
const scomplex czero = { 0.0f, 0.0f };
const dcomplex zzero = { 0.0 , 0.0  };

/* *************************************************************************

   FLA_Init()

 *************************************************************************** */

// Must be called with fla_init_lock held.
static void fla_init_locked( void )
{
  if ( FLA_initialized == TRUE ) return;

  // Static state: initialized on the first FLA_Init() of the process only.
  // FLA_Finalize() leaves it in place, so later re-initializations skip it.
  if ( fla_static_state_ready == FALSE )
  {
    FLA_Error_messages_init();

    FLA_Memory_leak_counter_init();

    fla_static_state_ready = TRUE;
  }

  // Idempotent: builds the constants on the first call only (except in SCC
  // builds, where they are recreated in shared memory for every session).
  FLA_Init_constants();

  // Dynamic state: rebuilt on every initialization.
  FLA_Cntl_init();

#if FLA_VECTOR_INTRINSIC_TYPE == FLA_SSE_INTRINSICS
  _MM_SET_FLUSH_ZERO_MODE( _MM_FLUSH_ZERO_ON );
#endif

#ifdef FLA_ENABLE_SUPERMATRIX
  FLASH_Queue_init();
#endif

  // Publish only once everything above is in place, so a lock-free
  // FLA_Initialized() never reports a half-built library.
  fla_store_initialized( TRUE );
}

void FLA_Init()
{
  fla_init_lock();

  fla_init_locked();

  // An explicit FLA_Init() means the user now owns the library lifetime:
  // in-flight library calls must not finalize it when they return.
  fla_auto_finalize = FALSE;

  fla_init_unlock();
}

/* *************************************************************************

  FLA_Finalize()

 *************************************************************************** */

// Must be called with fla_init_lock held.
static void fla_finalize_locked( void )
{
  if ( FLA_initialized == FALSE ) return;

  fla_store_initialized( FALSE );

  // No-op except in SCC builds: the constants are static and outlive
  // finalization (see FLA_Init_constants()).
  FLA_Finalize_constants();

  // Dynamic state is released.
  FLA_Cntl_finalize();

#ifdef FLA_ENABLE_SUPERMATRIX
  FLASH_Queue_finalize();
#endif

  // Reports the counter if enabled; the counter and its lock stay
  // initialized (static state).
  FLA_Memory_leak_counter_finalize();
}

void FLA_Finalize()
{
  fla_init_lock();

  if ( fla_safe_users > 0 )
  {
    // Other threads are inside libflame right now. Freeing the global
    // state would crash them, so hand finalization to the last of them.
    fla_auto_finalize = TRUE;
  }
  else
  {
    fla_finalize_locked();
    fla_auto_finalize = FALSE;
  }

  fla_init_unlock();
}

/* *************************************************************************

  FLA_Init_safe()

 *************************************************************************** */

void FLA_Init_safe( FLA_Error* init_result )
{
  fla_init_lock();

  if ( FLA_initialized == FALSE )
  {
    // Nobody holds the library initialized: initialize it on behalf of
    // the in-flight calls and finalize it when the last one returns.
    fla_init_locked();
    fla_auto_finalize = TRUE;
  }

  // Every caller is registered, including when the user initialized the
  // library, so that a concurrent FLA_Finalize() is deferred correctly.
  ++fla_safe_users;

  fla_init_unlock();

  // FLA_SUCCESS now means "registered; pass this to FLA_Finalize_safe()".
  // The value is only ever passed through to FLA_Finalize_safe().
  *init_result = FLA_SUCCESS;
}

/* *************************************************************************

  FLA_Finalize_safe()

 *************************************************************************** */

void FLA_Finalize_safe( FLA_Error init_result )
{
  if ( init_result != FLA_SUCCESS ) return;

  fla_init_lock();

  --fla_safe_users;

  if ( fla_safe_users == 0 && fla_auto_finalize == TRUE )
  {
    fla_finalize_locked();
    fla_auto_finalize = FALSE;
  }

  fla_init_unlock();
}

/* *************************************************************************

   FLA_Initialized()

 *************************************************************************** */

FLA_Bool FLA_Initialized( void )
{
  return fla_load_initialized();
}

/* *************************************************************************

   FLA_Init_constants()

   The global scalar constants (FLA_ONE, ..., FLA_OVERFLOW_SQUARE_THRES) are
   statically allocated: each FLA_Obj points to a base object and a buffer
   in static storage. They are built on the first call and then stay valid
   for the life of the process, including after FLA_Finalize(), so that
   re-initializing the library does not allocate and fill them again.

   Each constant is first built by the regular FLA_Obj_create_constant[_ext]
   routines, so every field and stored value (int, float, double, scomplex,
   dcomplex) is exactly what the library always produced, then its base
   object and buffer are copied into static storage and the heap copies are
   freed (leaving the memory leak counter balanced).

   In SCC builds the constant buffers must live in shared memory
   (FLA_shmalloc), so there the constants keep the original behaviour:
   created by every FLA_Init() and freed by every FLA_Finalize().

 *************************************************************************** */

// Must be called with fla_init_lock held (it is, from fla_init_locked()).
static FLA_Bool fla_constants_ready = FALSE;

#ifndef FLA_ENABLE_SCC

#define FLA_N_STATIC_CONSTANTS 18

// FLA_Obj_create_ext() may pad the 1x1 element to the configured alignment
// boundary (FLA_align_ldim() rounds the element up to a multiple of it), so
// the padded size never exceeds FLA_CONSTANT_SIZE + boundary.
#ifdef FLA_ENABLE_MEMORY_ALIGNMENT
  #define FLA_STATIC_CONSTANT_CAPACITY \
          ( FLA_CONSTANT_SIZE + FLA_MEMORY_ALIGNMENT_BOUNDARY )
  #if defined(__GNUC__)
    #define FLA_STATIC_CONSTANT_ALIGN \
            __attribute__(( aligned( FLA_MEMORY_ALIGNMENT_BOUNDARY ) ))
  #else
    #define FLA_STATIC_CONSTANT_ALIGN
  #endif
#else
  #define FLA_STATIC_CONSTANT_CAPACITY FLA_CONSTANT_SIZE
  #define FLA_STATIC_CONSTANT_ALIGN
#endif

// The union members give at least the alignment malloc() provides for every
// element type stored in a constant buffer.
typedef union
{
  char        bytes[ FLA_STATIC_CONSTANT_CAPACITY ];
  double      d;
  dcomplex    z;
  long double ld;
} FLA_STATIC_CONSTANT_ALIGN fla_static_constant_buf_t;

static FLA_Base_obj              fla_static_constant_base[ FLA_N_STATIC_CONSTANTS ];
static fla_static_constant_buf_t fla_static_constant_buf[ FLA_N_STATIC_CONSTANTS ];

// Compile-time assertion usable as a statement.
// C11 and later: _Static_assert. C99 (libflame's default -std=c99): a bit-field
// whose width is negative when cond is false; a bit-field width must also be
// an integer constant expression, so a non-constant cond is rejected as well.
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
  #define FLA_STATIC_ASSERT( cond, msg ) _Static_assert( cond, msg )
#else
  #define FLA_STATIC_ASSERT( cond, msg ) \
          ( void ) sizeof( struct { int fla_static_assert_failed : ( cond ) ? 1 : -1; } )
#endif

// Move the freshly created constant *obj into static slot idx.
//
// fla_constant_make_static() is a macro so that idx, a literal at every call
// site, can be checked at compile time: a slot index outside
// [0, FLA_N_STATIC_CONSTANTS-1] (or a non-constant index) fails to compile.
// A function parameter cannot be checked statically inside the function body.
#define fla_constant_make_static( obj, idx ) \
  do { \
    FLA_STATIC_ASSERT( ( idx ) >= 0 && ( idx ) <= FLA_N_STATIC_CONSTANTS - 1, \
                       "fla_constant_make_static: idx must be in " \
                       "[0, FLA_N_STATIC_CONSTANTS-1]" ); \
    fla_constant_make_static_impl( obj, idx ); \
  } while ( 0 )

static void fla_constant_make_static_impl( FLA_Obj* obj, int idx )
{
  FLA_Base_obj* heap_base = obj->base;
  FLA_Base_obj* base      = &fla_static_constant_base[ idx ];
  char*         buf       = fla_static_constant_buf[ idx ].bytes;
  size_t        n_bytes   = ( size_t ) heap_base->n_elem_alloc *
                            ( size_t ) FLA_CONSTANT_SIZE;

  // Should not happen given the capacity above; if a configuration ever
  // breaks the assumption, keep the heap object (still never freed, so the
  // constant still outlives FLA_Finalize()).
  if ( n_bytes > sizeof( fla_static_constant_buf[ idx ].bytes ) )
    return;
#ifdef FLA_ENABLE_MEMORY_ALIGNMENT
  if ( ( ( uintptr_t ) buf ) % ( uintptr_t ) FLA_MEMORY_ALIGNMENT_BOUNDARY != 0 )
    return;
#endif

  *base        = *heap_base;
  base->buffer = buf;
  base->id     = ( unsigned long ) base;

  // Copy the five typed values (not the raw bytes) into a zeroed buffer, so
  // the padding between them is deterministic instead of heap leftovers.
  // *obj still views the heap object here; the FLA_*_PTR macros resolve the
  // per-type offsets inside a constant's buffer.
  memset( buf, 0, n_bytes );
  {
    FLA_Obj tmp = *obj;

    tmp.base = base;

    *FLA_INT_PTR( tmp )            = *FLA_INT_PTR( *obj );
    *FLA_FLOAT_PTR( tmp )          = *FLA_FLOAT_PTR( *obj );
    *FLA_DOUBLE_PTR( tmp )         = *FLA_DOUBLE_PTR( *obj );
    *FLA_COMPLEX_PTR( tmp )        = *FLA_COMPLEX_PTR( *obj );
    *FLA_DOUBLE_COMPLEX_PTR( tmp ) = *FLA_DOUBLE_COMPLEX_PTR( *obj );
  }

  // Same release path as FLA_Obj_free() (non-SCC), without touching the
  // view fields of *obj, which stay those of the 1x1 constant.
  FLA_buff_free( heap_base->buffer );
  FLA_free( ( void * ) heap_base );

  obj->base = base;
}

#define FLA_CONSTANT_STATIC( value, obj, idx ) \
  do { FLA_Obj_create_constant( value, obj ); \
       fla_constant_make_static( obj, idx ); } while ( 0 )

#define FLA_CONSTANT_EXT_STATIC( value_s, value_d, obj, idx ) \
  do { FLA_Obj_create_constant_ext( value_s, value_d, obj ); \
       fla_constant_make_static( obj, idx ); } while ( 0 )

#else /* FLA_ENABLE_SCC: original per-session heap / shared-memory constants */

#define FLA_CONSTANT_STATIC( value, obj, idx ) \
  FLA_Obj_create_constant( value, obj )

#define FLA_CONSTANT_EXT_STATIC( value_s, value_d, obj, idx ) \
  FLA_Obj_create_constant_ext( value_s, value_d, obj )

#endif

void FLA_Init_constants()
{
  if ( fla_constants_ready == TRUE ) return;

  FLA_CONSTANT_STATIC(  3.0, &FLA_THREE,          0 );
  FLA_CONSTANT_STATIC(  2.0, &FLA_TWO,            1 );
  FLA_CONSTANT_STATIC(  1.0, &FLA_ONE,            2 );
  FLA_CONSTANT_STATIC(  0.5, &FLA_ONE_HALF,       3 );
  FLA_CONSTANT_STATIC(  0.0, &FLA_ZERO,           4 );
  FLA_CONSTANT_STATIC( -0.5, &FLA_MINUS_ONE_HALF, 5 );
  FLA_CONSTANT_STATIC( -1.0, &FLA_MINUS_ONE,      6 );
  FLA_CONSTANT_STATIC( -2.0, &FLA_MINUS_TWO,      7 );
  FLA_CONSTANT_STATIC( -3.0, &FLA_MINUS_THREE,    8 );


  { 
    float  
      eps_f, 
      sfmin_f = FLT_MIN, sfmin_f2,
      small_f = ( 1.0F / FLT_MAX ), 
      under_f = FLT_MIN, 
      over_f  = FLT_MAX;

    double 
      eps_d, 
      sfmin_d = DBL_MIN, sfmin_d2,
      small_d = ( 1.0  / DBL_MAX ), 
      under_d = DBL_MIN,
      over_d  = DBL_MAX;

    if ( FLT_ROUNDS == 1 )
    {
      eps_f = FLT_EPSILON*0.5F;
      eps_d = DBL_EPSILON*0.5;
    }
    else 
    {
      eps_f = FLT_EPSILON;
      eps_d = DBL_EPSILON;
    }

    if ( small_f >= sfmin_f ) sfmin_f = small_f * ( 1.0F + eps_f );
    if ( small_d >= sfmin_d ) sfmin_d = small_d * ( 1.0  + eps_d );

    sfmin_f  = sfmin_f/eps_f;
    sfmin_d  = sfmin_d/eps_d;

    sfmin_f2 = sqrt( sfmin_f );
    sfmin_d2 = sqrt( sfmin_d );

    FLA_CONSTANT_EXT_STATIC( eps_f,           eps_d,           &FLA_EPSILON,                 9 );

    FLA_CONSTANT_EXT_STATIC( sfmin_f,         sfmin_d,         &FLA_SAFE_MIN,               10 );
    FLA_CONSTANT_EXT_STATIC( 1.0F/sfmin_f,    1.0/sfmin_d,     &FLA_SAFE_INV_MIN,           11 );

    FLA_CONSTANT_EXT_STATIC( sfmin_f2,        sfmin_d2,        &FLA_SAFE_MIN_SQUARE,        12 );
    FLA_CONSTANT_EXT_STATIC( 1.0F/sfmin_f2,   1.0/sfmin_d2,    &FLA_SAFE_INV_MIN_SQUARE,    13 );
    
    FLA_CONSTANT_EXT_STATIC( under_f,         under_d,         &FLA_UNDERFLOW_THRES,        14 );
    FLA_CONSTANT_EXT_STATIC( over_f,          over_d,          &FLA_OVERFLOW_THRES,         15 );

    FLA_CONSTANT_EXT_STATIC( sqrt( under_f ), sqrt( under_d ), &FLA_UNDERFLOW_SQUARE_THRES, 16 );
    FLA_CONSTANT_EXT_STATIC( sqrt( over_f ),  sqrt( over_d ),  &FLA_OVERFLOW_SQUARE_THRES,  17 );
  } 

  fla_constants_ready = TRUE;
}

/* *************************************************************************

   FLA_Finalize_constants()

   The constants are statically allocated and remain valid after
   FLA_Finalize() (see FLA_Init_constants()), so this is a no-op, except in
   SCC builds where they are freed as before and rebuilt by the next
   FLA_Init_constants().

 *************************************************************************** */

void FLA_Finalize_constants()
{
#ifdef FLA_ENABLE_SCC
  FLA_Obj_free( &FLA_THREE );
  FLA_Obj_free( &FLA_TWO );
  FLA_Obj_free( &FLA_ONE );
  FLA_Obj_free( &FLA_ONE_HALF );
  FLA_Obj_free( &FLA_ZERO );
  FLA_Obj_free( &FLA_MINUS_ONE_HALF );
  FLA_Obj_free( &FLA_MINUS_ONE );
  FLA_Obj_free( &FLA_MINUS_TWO );
  FLA_Obj_free( &FLA_MINUS_THREE );

  FLA_Obj_free( &FLA_EPSILON );
  FLA_Obj_free( &FLA_SAFE_MIN );
  FLA_Obj_free( &FLA_SAFE_MIN_SQUARE );
  FLA_Obj_free( &FLA_SAFE_INV_MIN );
  FLA_Obj_free( &FLA_SAFE_INV_MIN_SQUARE );
  FLA_Obj_free( &FLA_UNDERFLOW_THRES );
  FLA_Obj_free( &FLA_OVERFLOW_THRES );  
  FLA_Obj_free( &FLA_UNDERFLOW_SQUARE_THRES );
  FLA_Obj_free( &FLA_OVERFLOW_SQUARE_THRES );

  fla_constants_ready = FALSE;
#endif
}

