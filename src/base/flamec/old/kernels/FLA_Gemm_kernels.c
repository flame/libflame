/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Gemm_initialized = FALSE;

void* blas_memory_alloc( void );
void* blas_memory_free( void* );

double* FLA_Work_buffer_aligned_A;
double* FLA_Work_buffer_aligned_B;

void FLA_Gemm_init( int mb, int kb )
{
  if( FLA_Gemm_initialized == TRUE )
    return;
 
  FLA_Work_buffer_aligned_A = ( double* ) blas_memory_alloc();
  FLA_Work_buffer_aligned_B = FLA_Work_buffer_aligned_A + mb * kb;
  //  FLA_Work_buffer_aligned_B = ( double* ) blas_memory_alloc();
 
  FLA_Gemm_initialized = TRUE;
}

void FLA_Gemm_finish()
{
  if( FLA_Gemm_initialized == FALSE )
    return;

  blas_memory_free( FLA_Work_buffer_aligned_A );
  //  blas_memory_free( FLA_Work_buffer_aligned_B );

  FLA_Gemm_initialized = FALSE;
}

void FLA_Gemm_pack_C( FLA_Trans transC, FLA_Obj C, FLA_Obj *packed_C )
{
  *packed_C = C;
}

void FLA_Gemm_unpack_andor_scale_C( FLA_Trans transC, FLA_Obj alpha, 
                                    FLA_Obj C, FLA_Obj* packed_C )
{
}

void FLA_Gemm_pack_andor_scale_B( FLA_Trans transB, FLA_Obj alpha,
                                  FLA_Obj B, FLA_Obj* packed_B )
{
  int m, n, ldim_B;
  double *buff_packed_B, *buff_packed_aligned_B, *buff_B;

  m      = FLA_Obj_length( B );
  n      = FLA_Obj_width ( B );
  ldim_B = FLA_Obj_ldim  ( B );
  buff_B = ( double* ) FLA_Obj_buffer_at_view( B );

  dgemm_oncopy( m, n, buff_B, ldim_B, FLA_Work_buffer_aligned_B );

  FLA_Obj_create_without_buffer( FLA_DOUBLE, m, n, packed_B );
  //  FLA_Obj_attach_buffer( buff_packed_B, m, packed_B );
}

void FLA_Gemm_release_pack_B( FLA_Trans transB, FLA_Obj* packed_B )
{
  //  blas_memory_free( FLA_Work_buffer_aligned_B );
}

// static int new_buffer_A = TRUE;

void FLA_Gemm_pack_andor_scale_A( FLA_Trans transA, FLA_Obj alpha,
                                  FLA_Obj A, FLA_Obj* packed_A )
{
  int m, n, ldim_A;
  double *buff_packed_A, *buff_packed_aligned_A, *buff_A;

  m      = FLA_Obj_length( A );
  n      = FLA_Obj_width ( A );
  ldim_A = FLA_Obj_ldim  ( A );
  buff_A = ( double* ) FLA_Obj_buffer_at_view( A );

  dgemm_itcopy( m, n, buff_A, ldim_A, FLA_Work_buffer_aligned_A );

  FLA_Obj_create_without_buffer( FLA_DOUBLE, m, n, packed_A );
  //  FLA_Obj_attach_buffer( buff_packed_A, n, packed_A );
}

void FLA_Gemm_release_pack_A( FLA_Trans transA, FLA_Obj* packed_A )
{
  //  blas_memory_free( FLA_Work_buffer_aligned_A );
}

void FLA_Gemm_kernel( FLA_Obj alpha, FLA_Obj packed_A, 
                      FLA_Obj packed_B, FLA_Obj packed_C )
{
  int m, n, k, ldim_C;
  double alpha_value, 
    *buff_A, *buff_B, *buff_C,
    *buff_aligned_A,   *buff_aligned_B;

  m      = FLA_Obj_length( packed_C );
  n      = FLA_Obj_width ( packed_C );
  k      = FLA_Obj_length( packed_A ); 
  ldim_C = FLA_Obj_ldim  ( packed_C );

  //  buff_A = ( double * ) FLA_Obj_buffer_at_view( packed_A );
  //  buff_B = ( double * ) FLA_Obj_buffer_at_view( packed_B );
  buff_C = ( double* ) FLA_Obj_buffer_at_view( packed_C );

  alpha_value = FLA_DOUBLE_VALUE( alpha );

  dgemm_kernel( m, n, k, alpha_value, 
                FLA_Work_buffer_aligned_A, 
                FLA_Work_buffer_aligned_B, buff_C, ldim_C );
}
