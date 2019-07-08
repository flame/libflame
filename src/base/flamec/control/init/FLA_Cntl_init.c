/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES 
void FLA_Cntl_init_ts( FLA_cntl_init_s* FLA_cntl_ts)
{
  FLA_cntl_ts->fla_error_checking_level = FLA_INTERNAL_ERROR_CHECKING_LEVEL;

  FLA_cntl_ts->FLA_Cntl_init_flamec_i = (FLA_Cntl_init_flamec_s *) FLA_malloc(sizeof(FLA_Cntl_init_flamec_s));
  FLA_cntl_ts->FLA_Cntl_init_flash_i= (FLA_Cntl_init_flash_s *) FLA_malloc(sizeof(FLA_Cntl_init_flash_s));

  FLA_Cntl_init_flamec_ts(FLA_cntl_ts->FLA_Cntl_init_flamec_i);
  FLA_Cntl_init_flash_ts(FLA_cntl_ts->FLA_Cntl_init_flash_i);

#ifdef FLA_ENABLE_SUPERMATRIX
  
  FLA_cntl_ts->flash_queue_initialized     = FALSE;

  FLA_cntl_ts->flash_queue_n_read_blocks   = 0;
  FLA_cntl_ts->flash_queue_n_write_blocks  = 0;

  FLA_cntl_ts->flash_queue_verbose         = FLASH_QUEUE_VERBOSE_NONE;
  FLA_cntl_ts->flash_queue_sorting         = FALSE;
  FLA_cntl_ts->flash_queue_caching         = FALSE;
  FLA_cntl_ts->flash_queue_work_stealing   = FALSE;
  FLA_cntl_ts->flash_queue_data_affinity   = FLASH_QUEUE_AFFINITY_NONE;
  
  FLA_cntl_ts->flash_queue_total_time      = 0.0;
  FLA_cntl_ts->flash_queue_parallel_time   = 0.0;
  FLA_cntl_ts->
  FLA_cntl_ts->flash_queue_block_size      = 0;
  FLA_cntl_ts->flash_queue_cache_size      = 2 * 1024 * 1024;
  FLA_cntl_ts->flash_queue_cache_line_size = 64;

  FLA_cntl_ts->flash_queue_cores_per_cache = 1;
  FLA_cntl_ts->flash_queue_cores_per_queue = 0;

#endif

  FLA_cntl_ts->flash_queue_stack           = 0;
  FLA_cntl_ts->flash_queue_enabled         = TRUE;

  FLA_cntl_ts->flash_queue_n_threads       = 1;
}

void FLA_Cntl_finalize_ts( FLA_cntl_init_s* FLA_cntl_ts)
{
  FLA_Cntl_finalize_flamec_ts(FLA_cntl_ts->FLA_Cntl_init_flamec_i);
  FLA_Cntl_finalize_flash_ts(FLA_cntl_ts->FLA_Cntl_init_flash_i);

  FLA_free(FLA_cntl_ts->FLA_Cntl_init_flamec_i);
  FLA_free(FLA_cntl_ts->FLA_Cntl_init_flash_i);
}
#endif

void FLA_Cntl_init()
{
  FLA_Cntl_init_flamec();
  FLA_Cntl_init_flash();
}

void FLA_Cntl_finalize()
{
  FLA_Cntl_finalize_flamec();
  FLA_Cntl_finalize_flash();
}

