/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#ifdef FLA_ENABLE_GPU

#include "cublas.h"


static FLA_Bool flash_queue_enabled_gpu  = FALSE;
static dim_t    flash_queue_gpu_n_blocks = 128;


void FLASH_Queue_init_gpu( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_init_gpu

----------------------------------------------------------------------------*/
{
   cublasInit();

   return;
}


void FLASH_Queue_finalize_gpu( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_finalize_gpu

----------------------------------------------------------------------------*/
{
   cublasShutdown();

   return;
}


FLA_Error FLASH_Queue_enable_gpu( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_enable_gpu

----------------------------------------------------------------------------*/
{
   if ( FLASH_Queue_stack_depth() == 0 && FLASH_Queue_get_enabled() )
   {
      // Enable if not begin parallel region yet and SuperMatrix is enabled.
      flash_queue_enabled_gpu = TRUE;
      return FLA_SUCCESS;
   }
   else
   {
      // Cannot change status during parallel region.
      return FLA_FAILURE;
   }
}


FLA_Error FLASH_Queue_disable_gpu( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_disable_gpu

----------------------------------------------------------------------------*/
{
   if ( FLASH_Queue_stack_depth() == 0 )
   {
      // Disable if not begin parallel region yet.
      flash_queue_enabled_gpu = FALSE;
      return FLA_SUCCESS;
   }
   else
   {
      // Cannot change status during parallel region.
      return FLA_FAILURE;
   }
}


FLA_Bool FLASH_Queue_get_enabled_gpu( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_enabled_gpu

----------------------------------------------------------------------------*/
{
   // Return if SuperMatrix is enabled, but always false if not.
   if ( FLASH_Queue_get_enabled() )
      return flash_queue_enabled_gpu;
   else
      return FALSE;
}


void FLASH_Queue_set_gpu_num_blocks( dim_t n_blocks )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_gpu_num_blocks

----------------------------------------------------------------------------*/
{
   flash_queue_gpu_n_blocks = n_blocks;

   return;
}


dim_t FLASH_Queue_get_gpu_num_blocks( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_gpu_num_blocks

----------------------------------------------------------------------------*/
{
   return flash_queue_gpu_n_blocks;
}


// --- helper functions --- ===================================================


FLA_Error FLASH_Queue_bind_gpu( integer thread )
/*----------------------------------------------------------------------------

   FLASH_Queue_bind_gpu

----------------------------------------------------------------------------*/
{
   // Bind a GPU to this thread.
   cudaSetDevice( thread );

   return FLA_SUCCESS;
}


FLA_Error FLASH_Queue_alloc_gpu( dim_t size,
                                 FLA_Datatype datatype,
                                 void** buffer_gpu )
/*----------------------------------------------------------------------------

   FLASH_Queue_alloc_gpu

----------------------------------------------------------------------------*/
{
   cublasStatus status;

   // Allocate memory for a block on GPU.
   status = cublasAlloc( size,
                         FLA_Obj_datatype_size( datatype ),
                         buffer_gpu );

   // Check to see if the allocation was successful.
   if ( status != CUBLAS_STATUS_SUCCESS )
      FLA_Check_error_code( FLA_MALLOC_GPU_RETURNED_NULL_POINTER );

   return FLA_SUCCESS;
}


FLA_Error FLASH_Queue_free_gpu( void* buffer_gpu )
/*----------------------------------------------------------------------------

   FLASH_Queue_free_gpu

----------------------------------------------------------------------------*/
{
   // Free memory for a block on GPU.
   cublasFree( buffer_gpu );

   return FLA_SUCCESS;
}


FLA_Error FLASH_Queue_write_gpu( FLA_Obj obj, void* buffer_gpu )
/*----------------------------------------------------------------------------

   FLASH_Queue_write_gpu

----------------------------------------------------------------------------*/
{
   // Write the contents of a block in main memory to GPU.
   cublasSetMatrix( FLA_Obj_length( obj ),
                    FLA_Obj_width( obj ),
                    FLA_Obj_datatype_size( FLA_Obj_datatype( obj ) ),
                    FLA_Obj_buffer_at_view( obj ), 
                    FLA_Obj_col_stride( obj ),
                    buffer_gpu,
                    FLA_Obj_length( obj ) );

   return FLA_SUCCESS;
}


FLA_Error FLASH_Queue_read_gpu( FLA_Obj obj, void* buffer_gpu )
/*----------------------------------------------------------------------------

   FLASH_Queue_read_gpu

----------------------------------------------------------------------------*/
{
   // Read the memory of a block on GPU to main memory.
   cublasGetMatrix( FLA_Obj_length( obj ),
                    FLA_Obj_width( obj ),
                    FLA_Obj_datatype_size( FLA_Obj_datatype( obj ) ),
                    buffer_gpu,
                    FLA_Obj_length( obj ),
                    FLA_Obj_buffer_at_view( obj ),
                    FLA_Obj_col_stride( obj ) );

   return FLA_SUCCESS;
}


void FLASH_Queue_exec_task_gpu( FLASH_Task* t, 
                                void** input_arg, 
                                void** output_arg )
/*----------------------------------------------------------------------------

   FLASH_Queue_exec_task_gpu

----------------------------------------------------------------------------*/
{
   // Define local function pointer types.

   // Level-3 BLAS
   typedef FLA_Error(*flash_gemm_gpu_p)(FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu, FLA_Obj beta, FLA_Obj C, void* C_gpu);
   typedef FLA_Error(*flash_hemm_gpu_p)(FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu, FLA_Obj beta, FLA_Obj C, void* C_gpu);
   typedef FLA_Error(*flash_herk_gpu_p)(FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj beta, FLA_Obj C, void* C_gpu);
   typedef FLA_Error(*flash_her2k_gpu_p)(FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu, FLA_Obj beta, FLA_Obj C, void* C_gpu);
   typedef FLA_Error(*flash_symm_gpu_p)(FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu, FLA_Obj beta, FLA_Obj C, void* C_gpu);
   typedef FLA_Error(*flash_syrk_gpu_p)(FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj beta, FLA_Obj C, void* C_gpu);
   typedef FLA_Error(*flash_syr2k_gpu_p)(FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu, FLA_Obj beta, FLA_Obj C, void* C_gpu);
   typedef FLA_Error(*flash_trmm_gpu_p)(FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj C, void* C_gpu);
   typedef FLA_Error(*flash_trsm_gpu_p)(FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj C, void* C_gpu);

   // Level-2 BLAS
   typedef FLA_Error(*flash_gemv_gpu_p)(FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu, FLA_Obj beta, FLA_Obj y, void* y_gpu);
   typedef FLA_Error(*flash_trsv_gpu_p)(FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu);

   // Level-1 BLAS
   typedef FLA_Error(*flash_axpy_gpu_p)(FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu);
   typedef FLA_Error(*flash_copy_gpu_p)(FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu);
   typedef FLA_Error(*flash_scal_gpu_p)(FLA_Obj alpha, FLA_Obj A, void* A_gpu);
   typedef FLA_Error(*flash_scalr_gpu_p)(FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, void* A_gpu);

   // Only execute task if it is not NULL.
   if ( t == NULL )
      return;

   // Now "switch" between the various possible task functions.

   // FLA_Gemm
   if ( t->func == (void *) FLA_Gemm_task )
   {
      flash_gemm_gpu_p func;
      func = (flash_gemm_gpu_p) FLA_Gemm_external_gpu;
      
      func( ( FLA_Trans   ) t->int_arg[0],
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
      flash_hemm_gpu_p func;
      func = (flash_hemm_gpu_p) FLA_Hemm_external_gpu;
      
      func( ( FLA_Side    ) t->int_arg[0],
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
      flash_herk_gpu_p func;
      func = (flash_herk_gpu_p) FLA_Herk_external_gpu;
      
      func( ( FLA_Uplo    ) t->int_arg[0],
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
      flash_her2k_gpu_p func;
      func = (flash_her2k_gpu_p) FLA_Her2k_external_gpu;
      
      func( ( FLA_Uplo     ) t->int_arg[0],
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
      flash_symm_gpu_p func;
      func = (flash_symm_gpu_p) FLA_Symm_external_gpu;
      
      func( ( FLA_Side    ) t->int_arg[0],
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
      flash_syrk_gpu_p func;
      func = (flash_syrk_gpu_p) FLA_Syrk_external_gpu;
      
      func( ( FLA_Uplo    ) t->int_arg[0],
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
      flash_syr2k_gpu_p func;
      func = (flash_syr2k_gpu_p) FLA_Syr2k_external_gpu;
      
      func( ( FLA_Uplo     ) t->int_arg[0],
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
      flash_trmm_gpu_p func;
      func = (flash_trmm_gpu_p) FLA_Trmm_external_gpu;
      
      func( ( FLA_Side    ) t->int_arg[0],
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
      flash_trsm_gpu_p func;
      func = (flash_trsm_gpu_p) FLA_Trsm_external_gpu;
      
      func( ( FLA_Side    ) t->int_arg[0],
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
      flash_gemv_gpu_p func;
      func = (flash_gemv_gpu_p) FLA_Gemv_external_gpu;
      
      func( ( FLA_Trans   ) t->int_arg[0],
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
      flash_trsv_gpu_p func;
      func = (flash_trsv_gpu_p) FLA_Trsv_external_gpu;
      
      func( ( FLA_Uplo    ) t->int_arg[0],
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
      flash_axpy_gpu_p func;
      func = (flash_axpy_gpu_p) FLA_Axpy_external_gpu;
         
      func(                 t->fla_arg[0],
                            t->input_arg[0],
                            input_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Copy
   else if ( t->func == (void *) FLA_Copy_task )
   {
      flash_copy_gpu_p func;
      func = (flash_copy_gpu_p) FLA_Copy_external_gpu;
         
      func(                 t->input_arg[0],
                            input_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Scal
   else if ( t->func == (void *) FLA_Scal_task )
   {
      flash_scal_gpu_p func;
      func = (flash_scal_gpu_p) FLA_Scal_external_gpu;

      func(                 t->fla_arg[0],
                            t->output_arg[0],
                            output_arg[0] );
   }
   // FLA_Scalr
   else if ( t->func == (void *) FLA_Scalr_task )
   {
      flash_scalr_gpu_p func;
      func = (flash_scalr_gpu_p) FLA_Scalr_external_gpu;

      func( ( FLA_Uplo    ) t->int_arg[0],
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


#endif // FLA_ENABLE_GPU
