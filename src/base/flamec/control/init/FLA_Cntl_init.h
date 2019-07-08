/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Cntl_init_flamec.h"
#include "FLA_Cntl_init_flash.h"


#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES

typedef struct FLA_cntl_init_s 
{
  FLA_Cntl_init_flamec_s *FLA_Cntl_init_flamec_i;
  FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i;

  FLA_Bool FLA_initialized;
  unsigned int fla_error_checking_level;

#ifdef FLA_ENABLE_SUPERMATRIX

  FLASH_Queue           _tq;

  FLA_Bool       flash_queue_initialized;
 
  int            flash_queue_n_read_blocks;
  int            flash_queue_n_write_blocks;
 
  FLASH_Verbose  flash_queue_verbose;
  FLA_Bool       flash_queue_sorting;
  FLA_Bool       flash_queue_caching;
  FLA_Bool       flash_queue_work_stealing;
  FLASH_Data_aff flash_queue_data_affinity;
 
  double         flash_queue_total_time;
  double         flash_queue_parallel_time;
 
  dim_t          flash_queue_block_size;
  dim_t          flash_queue_cache_size2;
  dim_t          flash_queue_cache_line_size;
 
  int            flash_queue_cores_per_cache;
#endif

  unsigned int   flash_queue_stack;
  FLA_Bool       flash_queue_enabled;
  unsigned int   flash_queue_n_threads;
}FLA_cntl_init_s;

void FLA_Cntl_init_ts( FLA_cntl_init_s* FLA_cntl_ts);
void FLA_Cntl_finalize_ts( FLA_cntl_init_s* FLA_cntl_ts);

#endif

void FLA_Cntl_init( void );
void FLA_Cntl_finalize( void );

