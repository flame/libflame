/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

static TLS_CLASS_SPEC integer      fla_mem_leak_counter;
static TLS_CLASS_SPEC FLA_Bool fla_mem_leak_counter_status;
#ifdef FLA_ENABLE_MULTITHREADING
static TLS_CLASS_SPEC FLA_Lock fla_mem_leak_counter_lock;
#endif


/* *************************************************************************

   FLA_Memory_leak_counter_init()

 *************************************************************************** */

void FLA_Memory_leak_counter_init( void )
{
  // Initialize the memory leak counter to zero.
  fla_mem_leak_counter = 0;

  // Initialize the memory leak counter status to whatever was requested at
  // configure-time.
#ifdef FLA_ENABLE_MEMORY_LEAK_COUNTER
  fla_mem_leak_counter_status = TRUE;
#else
  fla_mem_leak_counter_status = FALSE;
#endif

  // Initialize the memory leak counter lock, but only if we have locks in
  // the first place (ie: only if multithreading is enabled).
#ifdef FLA_ENABLE_MULTITHREADING
  FLA_Lock_init( &fla_mem_leak_counter_lock );
#endif
}

/* *************************************************************************

   FLA_Memory_leak_counter_finalize()

 *************************************************************************** */

void FLA_Memory_leak_counter_finalize( void )
{
  // Output the memory leak counter, but only if it's currently enabled.
  if ( FLA_Memory_leak_counter_status() == TRUE )
  {
    fprintf( stderr, "libflame: memory leak counter: %d\n", fla_mem_leak_counter );
    fflush( stderr );
  }

  // Destroy the memory leak counter lock, but only if we have locks in
  // the first place (ie: only if multithreading is enabled).
#ifdef FLA_ENABLE_MULTITHREADING
  FLA_Lock_destroy( &fla_mem_leak_counter_lock );
#endif

  // We leave the fla_mem_leak_counter_status variable alone.

  // Reset the counter, just for good measure.
  fla_mem_leak_counter = 0;
}

/* *************************************************************************

   FLA_Memory_leak_counter_status()

 *************************************************************************** */

FLA_Bool FLA_Memory_leak_counter_status( void )
{
  return fla_mem_leak_counter_status;
}

/* *************************************************************************

   FLA_Memory_leak_counter_set()

 *************************************************************************** */

FLA_Bool FLA_Memory_leak_counter_set( FLA_Bool new_status )
{
  FLA_Bool old_status;

  // Grab the current status.
  old_status = fla_mem_leak_counter_status;

  // Only make the change if the status is boolean. If the user provides us
  // with garbage, we do nothing.
  if ( new_status == TRUE || new_status == FALSE )
    fla_mem_leak_counter_status = new_status;

  return old_status;
}

/* **************************************************************************

   FLA_memset ()

 *************************************************************************** */

void* FLA_memset( void* str, integer c, uinteger len )                               // The memalign perf issue will not come here because we are taking void* (casting is expensive)
{
  unsigned char* ptr = str;
  while( len-- )
    *ptr++ = (unsigned char)c;

  return str;
}

/* ***************************************************************************

   FLA_malloc()

 *************************************************************************** */

void* FLA_malloc( size_t size )
{
  void*     ptr = NULL;
  FLA_Error e_val;
#ifdef FLA_ENABLE_MEMORY_ALIGNMENT
  integer       r_val;
#endif

  // In practice, the size argument should very rarely be zero. However, if the
  // calling code does request a memory region of zero length, we short-circut
  // the actual allocation request and just return NULL. Hopefully, the calling
  // code is written such that the pointer is never dereferenced. At free()-time
  // everything will be fine, as calling free() with a NULL pointer is safe.
  // Also note that we do NOT increment the memory leak counter before returning.
  // (Likewise, we will not decrement the counter when a NULL pointer is freed.)
  if ( size == 0 ) return NULL;

#ifdef FLA_ENABLE_MEMORY_ALIGNMENT

  // Allocate size bytes of memory. Here, we call posix_memalign() if
  // memory alignment was requested at configure-time, providing the
  // alignment boundary value given by the user. posix_memalign() also
  // returns an error code, which is how it signals that something
  // went wrong. Compare to malloc(), which does this by simply returning
  // a NULL pointer.
  r_val = posix_memalign( &ptr, ( size_t ) FLA_MEMORY_ALIGNMENT_BOUNDARY, size );

  // Check the return value of posix_memalign() for evidence that the
  // request failed.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
  {
    e_val = FLA_Check_posix_memalign_failure( r_val );
    FLA_Check_error_code( e_val );
  }

#else

  // Allocate size bytes of memory. Note that malloc() only guarantees 8-byte
  // alignment.
  ptr = malloc( size );

  // It may not seem useful to have a check for a null pointer here, given
  // that such an occurance would cause the file and line of the error to
  // be reported as the below line of the current file instead of the file
  // and line number of the calling code. However, consider that in the
  // unlikely event that malloc() does return a null pointer, the user will
  // have much bigger problems on his hands (e.g. an exhausted memory heap)
  // than needing to know exactly what line in the library triggered error.
  // Note that such a line in the application code is likely not the root
  // source of the problem anyway (ie: not the reason why the heap is full).
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
  {
    e_val = FLA_Check_malloc_pointer( ptr );
    FLA_Check_error_code( e_val );
  }

#endif

  // Update the memory leak counter if it is enabled, and do so thread-safely
  // if multithreading is enabled.
  if ( FLA_Memory_leak_counter_status() == TRUE )
  {
#ifdef FLA_ENABLE_MULTITHREADING
    FLA_Lock_acquire( &fla_mem_leak_counter_lock );
    fla_mem_leak_counter += 1;
    FLA_Lock_release( &fla_mem_leak_counter_lock );
#else
    fla_mem_leak_counter += 1;
#endif
  }
  
  // Return the pointer to the new memory region returned by malloc().
  return ptr;
}

/* ***************************************************************************

   FLA_realloc()

 *************************************************************************** */

void* FLA_realloc( void* old_ptr, size_t size )
{
  FLA_Error e_val;
  void*     new_ptr;

  // We can't do much if size is zero. To emulate realloc(), we must
  // return a NULL pointer, regardless of the value of old_ptr.
  if ( size == 0 )
  {
    // If the pointer is valid, free() it.
    if ( old_ptr != NULL )
      FLA_free( old_ptr );

    // If size is zero, we should return a NULL pointer.
    new_ptr = NULL;
  }
  else
  {
    // If old_ptr is NULL, allocate size bytes as if it were a first-time
    // FLA_malloc() request. Otherwise, proceed to realloc() the memory.
    if ( old_ptr == NULL )
    {
      new_ptr = FLA_malloc( size );
    }
    else
    {
      // At this point, we know that size is non-zero and old_ptr is valid.

      // Since we may need aligned addresses, we don't really want to call
      // realloc(), since it does not guarantee arbitrary aligned pointers.
      // But we can't implement it ourselves either, because we don't know
      // how large the original buffer is, therefor we don't know how much
      // to copy over after the new buffer is allocated. So we're stuck with
      // the system implementation.
      new_ptr = realloc( old_ptr, size );

      if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
      {
        e_val = FLA_Check_malloc_pointer( new_ptr );
        FLA_Check_error_code( e_val );
      }
    }
  }

  // Return the pointer (either NULL, or the return value from FLA_malloc()
  // or realloc()).
  return new_ptr;
}

/* ***************************************************************************

   FLA_free()

 *************************************************************************** */

void FLA_free( void* ptr )
{
  // We don't want to decrement the counter if the buffer is NULL.
  // This is because it's likely that the buffer was never allocated
  // a valid pointer to begin with, which means FLA_malloc() was never
  // called and thus the counter was never incremented. Or, it means
  // memory was allocated but the address has since been lost, which
  // means that we can't free it anyway.
  if ( ptr != NULL )
  {
    // Free the memory addressed by ptr.
    free( ptr );

    if ( FLA_Memory_leak_counter_status() == TRUE )
    {
#ifdef FLA_ENABLE_MULTITHREADING
      FLA_Lock_acquire( &fla_mem_leak_counter_lock );
      fla_mem_leak_counter -= 1;
      FLA_Lock_release( &fla_mem_leak_counter_lock );
#else
      fla_mem_leak_counter -= 1;
#endif
    }
  }
}

