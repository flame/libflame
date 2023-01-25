/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef FLASH_QUEUE_HIP_H
#define FLASH_QUEUE_HIP_H

#ifdef FLA_ENABLE_HIP


void           FLASH_Queue_init_hip( void );
void           FLASH_Queue_finalize_hip( void );

FLA_Error      FLASH_Queue_enable_hip( void );
FLA_Error      FLASH_Queue_disable_hip( void );
FLA_Bool       FLASH_Queue_get_enabled_hip( void );


// --- helper functions -------------------------------------------------------

FLA_Error      FLASH_Queue_available_devices_hip( int* device_count );

FLA_Error      FLASH_Queue_enable_malloc_managed_hip( void );
FLA_Error      FLASH_Queue_disable_malloc_managed_hip( void );
FLA_Bool       FLASH_Queue_get_malloc_managed_enabled_hip( void );

void           FLASH_Queue_set_hip_num_blocks( dim_t n_blocks );
dim_t          FLASH_Queue_get_hip_num_blocks( void );

FLA_Error      FLASH_Queue_bind_hip( int thread );
FLA_Error      FLASH_Queue_alloc_hip( dim_t size, FLA_Datatype datatype, void** buffer_hip );
FLA_Error      FLASH_Queue_free_async_hip( void* buffer_hip );
FLA_Error      FLASH_Queue_write_hip( FLA_Obj obj, void* buffer_hip );
FLA_Error      FLASH_Queue_read_hip( int thread, FLA_Obj obj, void* buffer_hip );
FLA_Error      FLASH_Queue_read_async_hip( int thread, FLA_Obj obj, void* buffer_hip );
FLA_Error      FLASH_Queue_sync_device_hip( int device );
FLA_Error      FLASH_Queue_sync_hip( );

void           FLASH_Queue_exec_task_hip( FLASH_Task* t, void** input_arg, void** output_arg );


#endif

#endif // FLASH_QUEUE_HIP_H
