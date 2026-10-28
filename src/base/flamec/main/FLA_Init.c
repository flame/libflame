/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

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

   Performance note: if no thread holds the library initialized, concurrent
   callers may still init/finalize repeatedly whenever the in-flight count
   touches zero. That is correct but costs time; calling FLA_Init() once at
   program start (FLA_INIT() from Fortran) avoids it.
*/

static FLA_Bool FLA_initialized   = FALSE;
static int      fla_safe_users    = 0;
static FLA_Bool fla_auto_finalize = FALSE;

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

  FLA_Error_messages_init();

  FLA_Memory_leak_counter_init();

  FLA_Init_constants();

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

  FLA_Finalize_constants();

  FLA_Cntl_finalize();

#ifdef FLA_ENABLE_SUPERMATRIX
  FLASH_Queue_finalize();
#endif

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

 *************************************************************************** */

void FLA_Init_constants()
{
  FLA_Obj_create_constant(  3.0, &FLA_THREE );
  FLA_Obj_create_constant(  2.0, &FLA_TWO );
  FLA_Obj_create_constant(  1.0, &FLA_ONE );
  FLA_Obj_create_constant(  0.5, &FLA_ONE_HALF );
  FLA_Obj_create_constant(  0.0, &FLA_ZERO );
  FLA_Obj_create_constant( -0.5, &FLA_MINUS_ONE_HALF );
  FLA_Obj_create_constant( -1.0, &FLA_MINUS_ONE );
  FLA_Obj_create_constant( -2.0, &FLA_MINUS_TWO );
  FLA_Obj_create_constant( -3.0, &FLA_MINUS_THREE );


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

    FLA_Obj_create_constant_ext( eps_f,           eps_d,           &FLA_EPSILON );

    FLA_Obj_create_constant_ext( sfmin_f,         sfmin_d,         &FLA_SAFE_MIN );
    FLA_Obj_create_constant_ext( 1.0F/sfmin_f,    1.0/sfmin_d,     &FLA_SAFE_INV_MIN );

    FLA_Obj_create_constant_ext( sfmin_f2,        sfmin_d2,        &FLA_SAFE_MIN_SQUARE );
    FLA_Obj_create_constant_ext( 1.0F/sfmin_f2,   1.0/sfmin_d2,    &FLA_SAFE_INV_MIN_SQUARE );
    
    FLA_Obj_create_constant_ext( under_f,         under_d,         &FLA_UNDERFLOW_THRES );
    FLA_Obj_create_constant_ext( over_f,          over_d,          &FLA_OVERFLOW_THRES );

    FLA_Obj_create_constant_ext( sqrt( under_f ), sqrt( under_d ), &FLA_UNDERFLOW_SQUARE_THRES );
    FLA_Obj_create_constant_ext( sqrt( over_f ),  sqrt( over_d ),  &FLA_OVERFLOW_SQUARE_THRES );
  } 
}

/* *************************************************************************

   FLA_Finalize_constants()

 *************************************************************************** */

void FLA_Finalize_constants()
{
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
}

