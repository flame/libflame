/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#ifdef FLA_ENABLE_HIP

#include "hip/hip_runtime.h"
#include "rocblas/rocblas.h"


static FLA_Bool flash_queue_enabled_hip  = FALSE;
static FLA_Bool flash_malloc_managed_hip = FALSE;
static dim_t    flash_queue_hip_n_blocks = 128;
static rocblas_handle* handles;

void FLASH_Queue_init_hip( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_init_hip

----------------------------------------------------------------------------*/
{

   int device_count;
   hipError_t err = hipGetDeviceCount( &device_count );
   if ( err != hipSuccess )
   {
     fprintf( stderr, "Failure to get device count: %d\n", err);
     return;
   }

   err = hipHostMalloc( (void**) &handles,
                        sizeof( rocblas_handle ) * device_count,
                        hipHostMallocDefault );
   if (err != hipSuccess )
   {
     fprintf( stderr, "Failure to malloc rocBLAS handles: %d\n", err);
     return;
   }
   for ( int i  = 0; i < device_count; i++ )
   {
      handles[i] = NULL;
   }

   return;
}


void FLASH_Queue_finalize_hip( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_finalize_hip

----------------------------------------------------------------------------*/
{

   int device_count;
   hipError_t err = hipGetDeviceCount( &device_count );
   if ( err != hipSuccess )
   {
     fprintf( stderr, "Failure to get device count: %d\n", err);
     return;
   }

   for ( int i  = 0; i < device_count; i++ )
   {
      // destroy rocBLAS handle
      if ( handles[i] != NULL )
         rocblas_destroy_handle( handles[i] );
   }

   hipHostFree( handles );

   return;
}


FLA_Error FLASH_Queue_available_devices_hip( int* device_count )
/*----------------------------------------------------------------------------

   FLASH_Queue_available_devices_hip

----------------------------------------------------------------------------*/
{

   hipError_t err = hipGetDeviceCount( device_count );
   if ( err != hipSuccess )
   {
     fprintf( stderr, "Failure to get device count: %d\n", err);
     return FLA_FAILURE;
   }

   return FLA_SUCCESS;
}


FLA_Error FLASH_Queue_enable_hip( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_enable_hip

----------------------------------------------------------------------------*/
{
   if ( FLASH_Queue_stack_depth() == 0 && FLASH_Queue_get_enabled() )
   {
      // Enable if not begin parallel region yet and SuperMatrix is enabled.
      flash_queue_enabled_hip = TRUE;
      return FLA_SUCCESS;
   }
   else
   {
      // Cannot change status during parallel region.
      return FLA_FAILURE;
   }
}


FLA_Error FLASH_Queue_disable_hip( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_disable_hip

----------------------------------------------------------------------------*/
{
   if ( FLASH_Queue_stack_depth() == 0 )
   {
      // Disable if not begin parallel region yet.
      flash_queue_enabled_hip = FALSE;
      return FLA_SUCCESS;
   }
   else
   {
      // Cannot change status during parallel region.
      return FLA_FAILURE;
   }
}


FLA_Bool FLASH_Queue_get_enabled_hip( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_enabled_hip

----------------------------------------------------------------------------*/
{
   // Return if SuperMatrix is enabled, but always false if not.
   if ( FLASH_Queue_get_enabled() )
      return flash_queue_enabled_hip;
   else
      return FALSE;
}

FLA_Error FLASH_Queue_enable_malloc_managed_hip( void )
/*---------------------------------------------------------------------------

   FLASH_Queue_enable_malloc_managed_hip

----------------------------------------------------------------------------*/
{
   if ( FLASH_Queue_stack_depth() == 0 && FLASH_Queue_get_enabled() )
   {
      // Enable if not begin parallel region yet and SuperMatrix is enabled.
      flash_malloc_managed_hip = TRUE;
      return FLA_SUCCESS;
   }
   else
   {
      // Cannot change status during parallel region.
      return FLA_FAILURE;
   }
}

FLA_Error FLASH_Queue_disable_malloc_managed_hip( void )
/*---------------------------------------------------------------------------

   FLASH_Queue_disable_malloc_managed_hip

----------------------------------------------------------------------------*/
{
   if ( FLASH_Queue_stack_depth() == 0 && FLASH_Queue_get_enabled() )
   {
      // Enable if not begin parallel region yet and SuperMatrix is enabled.
      flash_malloc_managed_hip = FALSE;
      return FLA_SUCCESS;
   }
   else
   {
      // Cannot change status during parallel region.
      return FLA_FAILURE;
   }
}

FLA_Bool FLASH_Queue_get_malloc_managed_enabled_hip( void )
/*---------------------------------------------------------------------------

   FLASH_Queue_get_malloc_managed_enabled_hip

----------------------------------------------------------------------------*/
{
   // Return if SuperMatrix is enabled, but always false if not.
   if ( FLASH_Queue_get_enabled() )
      return flash_malloc_managed_hip;
   else
      return FALSE;
}

void FLASH_Queue_set_hip_num_blocks( dim_t n_blocks )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_hip_num_blocks

----------------------------------------------------------------------------*/
{
   flash_queue_hip_n_blocks = n_blocks;

   return;
}


dim_t FLASH_Queue_get_hip_num_blocks( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_hip_num_blocks

----------------------------------------------------------------------------*/
{
   return flash_queue_hip_n_blocks;
}


// --- helper functions --- ===================================================


FLA_Error FLASH_Queue_bind_hip( int thread )
/*----------------------------------------------------------------------------

   FLASH_Queue_bind_hip

----------------------------------------------------------------------------*/
{

   // Bind a HIP device to this thread.
   if ( hipSetDevice( thread ) != hipSuccess ) return FLA_FAILURE;

   // initialize its rocBLAS handle
   if ( handles[thread] == NULL )
   {
      //hipStream_t stream;
      //if ( hipStreamCreate(&stream) != hipSuccess ) return FLA_FAILURE;
      if ( rocblas_create_handle( &(handles[thread]) ) != rocblas_status_success ) return FLA_FAILURE;
      //if ( rocblas_set_stream( handles[thread], stream ) != rocblas_status_success ) return FLA_FAILURE;
   }

   return FLA_SUCCESS;
}


FLA_Error FLASH_Queue_alloc_async_hip( int thread,
                                 dim_t size,
                                 FLA_Datatype datatype,
                                 void** buffer_hip )
/*----------------------------------------------------------------------------

   FLASH_Queue_alloc_hip

----------------------------------------------------------------------------*/
{

   // Allocate memory for a block on HIP.
   hipStream_t stream;
   rocblas_get_stream( handles[thread], &stream );
   //fprintf(stdout, "Trying to allocate %ld bytes on device %d\n", size, thread);
   hipError_t status = hipMallocAsync( buffer_hip,
                                       size * FLA_Obj_datatype_size( datatype ),
                                       stream );

   // Check to see if the allocation was successful.
   if ( status != hipSuccess )
   {
      int device;
      hipGetDevice( &device );
      fprintf(stderr, "Error err=%d on device=%d\n", status, device);
      FLA_Check_error_code( FLA_MALLOC_GPU_RETURNED_NULL_POINTER );
   }

   //fprintf( stdout, "allocating on thread %d for pointer %p\n", thread, *buffer_hip );

   return FLA_SUCCESS;
}


FLA_Error FLASH_Queue_free_async_hip( int thread, void* buffer_hip )
/*----------------------------------------------------------------------------

   FLASH_Queue_free_async_hip

----------------------------------------------------------------------------*/
{
   // Free memory for a block on HIP.
   hipStream_t stream;
   rocblas_get_stream( handles[thread], &stream );
   hipFreeAsync( stream, buffer_hip );

   return FLA_SUCCESS;
}


FLA_Error FLASH_Queue_write_async_hip( int thread,
                                       FLA_Obj obj,
                                       void* buffer_hip )
/*----------------------------------------------------------------------------

   FLASH_Queue_write_async_hip

----------------------------------------------------------------------------*/
{
   if ( flash_malloc_managed_hip )
   {
     return FLA_SUCCESS; // HMM will take care of getting the memory over
   }
   // Write the contents of a block in main memory to HIP.
   const size_t count = FLA_Obj_elem_size( obj )
                          * FLA_Obj_col_stride( obj )
                          * FLA_Obj_width( obj );
   hipStream_t stream;
   rocblas_get_stream( handles[thread], &stream );
   const hipError_t err = hipMemcpyAsync( buffer_hip,
                                          FLA_Obj_buffer_at_view( obj ),
                                          count,
                                          hipMemcpyHostToDevice,
                                          stream );

   if ( err != hipSuccess )
   {
     fprintf( stderr,
              "Failure to write block to HIP device. Size=%ld, err=%d\n",
              count, err );
     return FLA_FAILURE;
   }

   return FLA_SUCCESS;
}


FLA_Error FLASH_Queue_read_hip( int thread, FLA_Obj obj, void* buffer_hip )
/*----------------------------------------------------------------------------

   FLASH_Queue_read_hip

----------------------------------------------------------------------------*/
{

   FLA_Error err1 = FLASH_Queue_read_async_hip( thread, obj, buffer_hip );
   if ( err1 != FLA_SUCCESS ) return err1;

   return FLASH_Queue_sync_stream_hip( thread );
}

FLA_Error FLASH_Queue_read_async_hip( int thread, FLA_Obj obj, void* buffer_hip )
/*----------------------------------------------------------------------------

   FLASH_Queue_read_async_hip

----------------------------------------------------------------------------*/
{
   if ( flash_malloc_managed_hip )
   {
     return FLA_SUCCESS;
   }

   // Read the memory of a block on HIP to main memory.
   hipSetDevice( thread );
   const size_t count = FLA_Obj_elem_size( obj )
                          * FLA_Obj_col_stride( obj )
                          * FLA_Obj_width( obj );
   hipStream_t stream;
   rocblas_get_stream( handles[thread], &stream );
   const hipError_t err = hipMemcpyAsync( FLA_Obj_buffer_at_view( obj ),
                                          buffer_hip,
                                          count,
                                          hipMemcpyDeviceToHost,
                                          stream );

   if ( err != hipSuccess )
   {
     fprintf( stderr,
              "Failure to async read block from HIP device. Size=%ld, err=%d\n",
              count, err );
     return FLA_FAILURE;
   }

   return FLA_SUCCESS;
}

FLA_Error FLASH_Queue_sync_device_hip( int device )
/*----------------------------------------------------------------------------

   FLASH_Queue_sync_device_hip

----------------------------------------------------------------------------*/
{
   hipSetDevice( device );
   const hipError_t err = hipDeviceSynchronize( );
   if ( err != hipSuccess )
   {
     fprintf( stderr,
              "Failure to sync HIP device. Device=%d, err=%d\n",
              device, err );
     return FLA_FAILURE;
   }

   return FLA_SUCCESS;
}

FLA_Error FLASH_Queue_sync_stream_hip( int thread )
/*----------------------------------------------------------------------------

   FLASH_Queue_sync_stream_hip

----------------------------------------------------------------------------*/
{
   hipStream_t stream;
   rocblas_get_stream( handles[thread], &stream );
   const hipError_t err = hipStreamSynchronize( stream );
   if ( err != hipSuccess )
   {
     fprintf( stderr,
              "Failure to sync HIP stream. Thread=%d, err=%d\n",
              thread, err );
     return FLA_FAILURE;
   }

   return FLA_SUCCESS;
}


FLA_Error FLASH_Queue_sync_hip( )
/*----------------------------------------------------------------------------

   FLASH_Queue_sync_hip

----------------------------------------------------------------------------*/
{
   const hipError_t err = hipStreamSynchronize( (hipStream_t) 0 );
   if ( err != hipSuccess )
   {
     fprintf( stderr,
              "Failure to sync HIP null stream. Err=%d\n",
              err );
     return FLA_FAILURE;
   }

   return FLA_SUCCESS;
}


void FLASH_Queue_exec_task_hip( FLASH_Task* t, 
                                void** input_arg, 
                                void** output_arg )
/*----------------------------------------------------------------------------

   FLASH_Queue_exec_task_hip

----------------------------------------------------------------------------*/
{
   // Define local function pointer types.

   // LAPACK
   typedef FLA_Error(*flash_chol_hip_p)(rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip );
   typedef FLA_Error(*flash_trinv_hip_p)(rocblas_handle handle, FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A, void* A_hip );
   typedef FLA_Error(*flash_eig_gest_hip_p)(rocblas_handle handle, FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
   typedef FLA_Error(*flash_lu_piv_hip_p)(rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj p );
   typedef FLA_Error(*flash_lu_piv_copy_hip_p)(rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj p, FLA_Obj U, void* U_hip );
   typedef FLA_Error(*flash_sa_fs_hip_p)(rocblas_handle handle, FLA_Obj L, FLA_Obj D, void* D_hip, FLA_Obj p, FLA_Obj C, void* C_hip,FLA_Obj E, void* E_hip, dim_t nb_alg );
   typedef FLA_Error(*flash_trsm_piv_hip_p)(rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj p );

   // Level-3 BLAS
   typedef FLA_Error(*flash_gemm_hip_p)(rocblas_handle handle, FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj beta, FLA_Obj C, void* C_hip);
   typedef FLA_Error(*flash_hemm_hip_p)(rocblas_handle handle, FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj beta, FLA_Obj C, void* C_hip);
   typedef FLA_Error(*flash_herk_hip_p)(rocblas_handle handle, FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj beta, FLA_Obj C, void* C_hip);
   typedef FLA_Error(*flash_her2k_hip_p)(rocblas_handle handle, FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj beta, FLA_Obj C, void* C_hip);
   typedef FLA_Error(*flash_symm_hip_p)(rocblas_handle handle, FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj beta, FLA_Obj C, void* C_hip);
   typedef FLA_Error(*flash_syrk_hip_p)(rocblas_handle handle, FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj beta, FLA_Obj C, void* C_hip);
   typedef FLA_Error(*flash_syr2k_hip_p)(rocblas_handle handle, FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj beta, FLA_Obj C, void* C_hip);
   typedef FLA_Error(*flash_trmm_hip_p)(rocblas_handle handle, FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj C, void* C_hip);
   typedef FLA_Error(*flash_trsm_hip_p)(rocblas_handle handle, FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj C, void* C_hip);

   // Level-2 BLAS
   typedef FLA_Error(*flash_gemv_hip_p)(rocblas_handle handle, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj x, void* x_hip, FLA_Obj beta, FLA_Obj y, void* y_hip);
   typedef FLA_Error(*flash_trsv_hip_p)(rocblas_handle handle, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, void* A_hip, FLA_Obj x, void* x_hip);

   // Level-1 BLAS
   typedef FLA_Error(*flash_axpy_hip_p)(rocblas_handle handle, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip);
   typedef FLA_Error(*flash_copy_hip_p)(rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip);
   typedef FLA_Error(*flash_copyr_hip_p)(rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip);
   typedef FLA_Error(*flash_scal_hip_p)(rocblas_handle handle, FLA_Obj alpha, FLA_Obj A, void* A_hip);
   typedef FLA_Error(*flash_scalr_hip_p)(rocblas_handle handle, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, void* A_hip);

   // Only execute task if it is not NULL.
   if ( t == NULL )
      return;

   // Now "switch" between the various possible task functions.
   hipSetDevice( t->thread );
   rocblas_handle handle = handles[t->thread];
   // FLA_Chol
   if ( t-> func == (void*) FLA_Chol_task )
   {
      flash_chol_hip_p func;
      func = (flash_chol_hip_p) FLA_Chol_blk_external_hip;

      func(
                            handle,
            ( FLA_Uplo    ) t->int_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Trinv
   else if ( t-> func == (void*) FLA_Trinv_task )
   {
      flash_trinv_hip_p func;
      func = (flash_trinv_hip_p) FLA_Trinv_blk_external_hip;

      func(
                            handle,
            ( FLA_Uplo    ) t->int_arg[0],
            ( FLA_Diag    ) t->int_arg[1],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Eig_gest
   else if ( t->func == (void *) FLA_Eig_gest_task )
   {
      flash_eig_gest_hip_p func;
      func = (flash_eig_gest_hip_p) FLA_Eig_gest_blk_external_hip;

      func(                     handle,
            ( FLA_Inv         ) t->int_arg[0],
            ( FLA_Uplo        ) t->int_arg[1],
                                t->output_arg[1],
                                output_arg[1],
                                t->input_arg[0],
                                input_arg[0] );
   }
   // FLA_LU_piv
   else if ( t->func == (void *) FLA_LU_piv_task )
   {
      flash_lu_piv_hip_p func;
      func = (flash_lu_piv_hip_p) FLA_LU_piv_blk_external_hip;

      func(               handle,
                          t->output_arg[0],
                          output_arg[0],
                          t->fla_arg[0] );
   }
   // FLA_LU_piv_copy
   else if ( t->func == (void *) FLA_LU_piv_copy_task )
   {
      flash_lu_piv_copy_hip_p func;
      func = (flash_lu_piv_copy_hip_p) FLA_LU_piv_copy_external_hip;

      func(               handle,
                          t->output_arg[0],
                          output_arg[0],
                          t->fla_arg[0],
                          t->output_arg[1],
                          output_arg[1] );
   }
   // FLA_SA_FS
   else if ( t->func == (void *) FLA_SA_FS_task )
   {
      flash_sa_fs_hip_p func;
      func = (flash_sa_fs_hip_p) FLA_SA_FS_blk_hip;

      func(                 handle,
                            t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->fla_arg[1],
                            t->output_arg[1],
                            output_arg[1],
                            t->output_arg[0],
                            output_arg[0],
                            t->int_arg[0] );
   }
   // FLA_Trsm_piv
   else if ( t->func == (void *) FLA_Trsm_piv_task )
   {
      flash_trsm_piv_hip_p func;
      func = (flash_trsm_piv_hip_p) FLA_Trsm_piv_external_hip;

      func(               handle,
                          t->input_arg[0],
                          input_arg[0],
                          t->output_arg[0],
                          output_arg[0],
                          t->fla_arg[0]
      );
   }
   // FLA_Gemm
   else if ( t->func == (void *) FLA_Gemm_task )
   {
      flash_gemm_hip_p func;
      func = (flash_gemm_hip_p) FLA_Gemm_external_hip;
      
      func(                 handle,
            ( FLA_Trans   ) t->int_arg[0],
            ( FLA_Trans   ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->input_arg[1],
                            input_arg[1],
                            t->fla_arg[1],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Hemm
   else if ( t->func == (void *) FLA_Hemm_task )
   {
      flash_hemm_hip_p func;
      func = (flash_hemm_hip_p) FLA_Hemm_external_hip;
      
      func(                 handle,
            ( FLA_Side    ) t->int_arg[0],
            ( FLA_Uplo    ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->input_arg[1],
                            input_arg[1],
                            t->fla_arg[1],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Herk
   else if ( t->func == (void *) FLA_Herk_task )
   {
      flash_herk_hip_p func;
      func = (flash_herk_hip_p) FLA_Herk_external_hip;
      
      func(                 handle,
            ( FLA_Uplo    ) t->int_arg[0],
            ( FLA_Trans   ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->fla_arg[1],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Her2k
   else if ( t->func == (void *) FLA_Her2k_task )
   {
      flash_her2k_hip_p func;
      func = (flash_her2k_hip_p) FLA_Her2k_external_hip;
      
      func(                  handle,
            ( FLA_Uplo     ) t->int_arg[0],
            ( FLA_Trans    ) t->int_arg[1],
                             t->fla_arg[0],
                             t->input_arg[0],
                             input_arg[0],
                             t->input_arg[1],
                             input_arg[1],
                             t->fla_arg[1],
                             t->output_arg[0],
                             output_arg[0] );
   }
   // FLA_Symm
   else if ( t->func == (void *) FLA_Symm_task )
   {
      flash_symm_hip_p func;
      func = (flash_symm_hip_p) FLA_Symm_external_hip;
      
      func(                 handle,
            ( FLA_Side    ) t->int_arg[0],
            ( FLA_Uplo    ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->input_arg[1],
                            input_arg[1],
                            t->fla_arg[1],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Syrk
   else if ( t->func == (void *) FLA_Syrk_task )
   {
      flash_syrk_hip_p func;
      func = (flash_syrk_hip_p) FLA_Syrk_external_hip;
      
      func(                 handle,
            ( FLA_Uplo    ) t->int_arg[0],
            ( FLA_Trans   ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->fla_arg[1],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Syr2k
   else if ( t->func == (void *) FLA_Syr2k_task )
   {
      flash_syr2k_hip_p func;
      func = (flash_syr2k_hip_p) FLA_Syr2k_external_hip;
      
      func(                  handle,
            ( FLA_Uplo     ) t->int_arg[0],
            ( FLA_Trans    ) t->int_arg[1],
                             t->fla_arg[0],
                             t->input_arg[0],
                             input_arg[0],
                             t->input_arg[1],
                             input_arg[1],
                             t->fla_arg[1],
                             t->output_arg[0],
                             output_arg[0] );
   }
   // FLA_Trmm
   else if ( t->func == (void *) FLA_Trmm_task )
   {
      flash_trmm_hip_p func;
      func = (flash_trmm_hip_p) FLA_Trmm_external_hip;
      
      func(                 handle,
            ( FLA_Side    ) t->int_arg[0],
            ( FLA_Uplo    ) t->int_arg[1],
            ( FLA_Trans   ) t->int_arg[2],
            ( FLA_Diag    ) t->int_arg[3],
                            t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Trsm
   else if ( t->func == (void *) FLA_Trsm_task )
   {
      flash_trsm_hip_p func;
      func = (flash_trsm_hip_p) FLA_Trsm_external_hip;
      
      func(                 handle,
            ( FLA_Side    ) t->int_arg[0],
            ( FLA_Uplo    ) t->int_arg[1],
            ( FLA_Trans   ) t->int_arg[2],
            ( FLA_Diag    ) t->int_arg[3],
                            t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Gemv
   else if ( t->func == (void *) FLA_Gemv_task )
   {
      flash_gemv_hip_p func;
      func = (flash_gemv_hip_p) FLA_Gemv_external_hip;
      
      func(                 handle,
            ( FLA_Trans   ) t->int_arg[0],
                            t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->input_arg[1],
                            input_arg[1],
                            t->fla_arg[1],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Trsv
   else if ( t->func == (void *) FLA_Trsv_task )
   {
      flash_trsv_hip_p func;
      func = (flash_trsv_hip_p) FLA_Trsv_external_hip;
      
      func(                 handle,
            ( FLA_Uplo    ) t->int_arg[0],
            ( FLA_Trans   ) t->int_arg[1],
            ( FLA_Diag    ) t->int_arg[2],
                            t->input_arg[0],
                            input_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Axpy
   else if ( t->func == (void *) FLA_Axpy_task )
   {
      flash_axpy_hip_p func;
      func = (flash_axpy_hip_p) FLA_Axpy_external_hip;
         
      func(                 handle,
                            t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Copy
   else if ( t->func == (void *) FLA_Copy_task )
   {
      flash_copy_hip_p func;
      func = (flash_copy_hip_p) FLA_Copy_external_hip;
         
      func(                 handle,
                            t->input_arg[0],
                            input_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Copyr
   else if ( t->func == (void *) FLA_Copyr_task )
   {
      flash_copyr_hip_p func;
      func = (flash_copyr_hip_p) FLA_Copyr_external_hip;

      func(                  handle,
            ( FLA_Uplo     ) t->int_arg[0],
                             t->input_arg[0],
                             input_arg[0],
                             t->output_arg[0],
                             output_arg[0] );
   }
   // FLA_Scal
   else if ( t->func == (void *) FLA_Scal_task )
   {
      flash_scal_hip_p func;
      func = (flash_scal_hip_p) FLA_Scal_external_hip;

      func(                 handle,
                            t->fla_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Scalr
   else if ( t->func == (void *) FLA_Scalr_task )
   {
      flash_scalr_hip_p func;
      func = (flash_scalr_hip_p) FLA_Scalr_external_hip;

      func(                 handle,
            ( FLA_Uplo    ) t->int_arg[0],
                            t->fla_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   else
   {
      FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
   }

   return;
}


#endif // FLA_ENABLE_HIP
