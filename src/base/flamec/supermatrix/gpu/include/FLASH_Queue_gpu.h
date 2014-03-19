/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef FLASH_QUEUE_GPU_H
#define FLASH_QUEUE_GPU_H

#ifdef FLA_ENABLE_GPU


void           FLASH_Queue_init_gpu( void );
void           FLASH_Queue_finalize_gpu( void );

FLA_Error      FLASH_Queue_enable_gpu( void );
FLA_Error      FLASH_Queue_disable_gpu( void );
FLA_Bool       FLASH_Queue_get_enabled_gpu( void );


// --- helper functions -------------------------------------------------------

void           FLASH_Queue_set_gpu_num_blocks( dim_t n_blocks );
dim_t          FLASH_Queue_get_gpu_num_blocks( void );

FLA_Error      FLASH_Queue_bind_gpu( int thread );
FLA_Error      FLASH_Queue_alloc_gpu( dim_t size, FLA_Datatype datatype, void** buffer_gpu );
FLA_Error      FLASH_Queue_free_gpu( void* buffer_gpu );
FLA_Error      FLASH_Queue_write_gpu( FLA_Obj obj, void* buffer_gpu );
FLA_Error      FLASH_Queue_read_gpu( FLA_Obj obj, void* buffer_gpu );

void           FLASH_Queue_exec_task_gpu( FLASH_Task* t, void** input_arg, void** output_arg );


#endif

#endif // FLASH_QUEUE_GPU_H
