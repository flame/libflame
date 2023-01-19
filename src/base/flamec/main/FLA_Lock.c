/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#ifdef FLA_ENABLE_MULTITHREADING


#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
#ifdef FLA_ENABLE_TIDSP
#include <ti/omp/omp.h>
#else
#include <omp.h>
#endif
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
#include <pthread.h>
#endif


void FLA_Lock_init( FLA_Lock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_init

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_init_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_mutex_init( &(fla_lock_ptr->lock), NULL );
#endif
}


void FLA_Lock_acquire( FLA_Lock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_acquire

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_set_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_mutex_lock( &(fla_lock_ptr->lock) );
#endif
}


void FLA_Lock_release( FLA_Lock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_release

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_unset_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_mutex_unlock( &(fla_lock_ptr->lock) );
#endif
}


void FLA_Lock_destroy( FLA_Lock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_destroy

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_destroy_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_mutex_destroy( &(fla_lock_ptr->lock) );
#endif
}


void FLA_RWLock_init( FLA_RWLock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_init

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_init_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_rwlock_init( &(fla_lock_ptr->lock), NULL );
#endif
}


void FLA_RWLock_write_acquire( FLA_RWLock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_acquire

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_set_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_rwlock_wrlock( &(fla_lock_ptr->lock) );
#endif
}


void FLA_RWLock_read_acquire( FLA_RWLock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_read_acquire

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_set_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_rwlock_rdlock( &(fla_lock_ptr->lock) );
#endif
}


void FLA_RWLock_release( FLA_RWLock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_release

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_unset_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_rwlock_unlock( &(fla_lock_ptr->lock) );
#endif
}


void FLA_RWLock_destroy( FLA_RWLock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_destroy

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_destroy_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_rwlock_destroy( &(fla_lock_ptr->lock) );
#endif
}

#endif // FLA_ENABLE_MULTITHREADING

