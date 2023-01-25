/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include "hip/hip_runtime_api.h"
#include "rocblas/rocblas.h"
#include "rocsolver/rocsolver.h"

FLA_Error FLA_Bidiag_apply_U_external_hip( rocblas_handle handle, FLA_Side side, FLA_Trans trans, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip, FLA_Obj B, void* B_hip )
{
  FLA_Datatype datatype;
  dim_t        n_A;
  dim_t        m_B, n_B;
  dim_t        cs_A;
  dim_t        cs_B;
  dim_t        k_t;
  rocblas_storev blas_vect = rocblas_column_wise;

  //if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
  //  FLA_Apply_Q_check( side, trans, storev, A, t, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  // m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  cs_B     = FLA_Obj_col_stride( B );

  if ( blas_vect == rocblas_column_wise ) k_t = FLA_Obj_vector_dim( t );
  else                    k_t = FLA_Obj_vector_dim( t ) + 1;

  if ( FLA_Obj_is_real( A ) && trans == FLA_CONJ_TRANSPOSE )
    trans = FLA_TRANSPOSE;

  void* A_vecs = NULL;
  void* t_scals = NULL;
  void* B_mat = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_vecs = FLA_Obj_buffer_at_view( A );
    t_scals = FLA_Obj_buffer_at_view( t );
    B_mat = FLA_Obj_buffer_at_view( B );
  }
  else
  {
    A_vecs = A_hip;
    t_scals = t_hip;
    B_mat = B_hip;
  }

  FLA_Trans trans_a_corr = trans;
  FLA_Bool conj_no_trans_a = FALSE;

  void* A_mat_corr = NULL;
  if ( FLA_Obj_is_complex( A ) && trans == FLA_CONJ_NO_TRANSPOSE )
  {
    // must correct by copying to temporary buffer and conjugating there
    trans_a_corr = FLA_NO_TRANSPOSE;
    conj_no_trans_a = TRUE;

    dim_t elem_size = FLA_Obj_elem_size( A );
    size_t count = elem_size * cs_A * n_A;
    hipStream_t stream;
    rocblas_get_stream( handle, &stream );
    hipMallocAsync( &A_mat_corr, count, stream );
    FLA_Copyconj_general_external_hip( handle, A, A_hip, A_mat_corr );
    A_vecs = A_mat_corr;
  }

  rocblas_side blas_side = FLA_Param_map_flame_to_rocblas_side( side );
  rocblas_operation blas_trans = FLA_Param_map_flame_to_rocblas_trans( trans_a_corr, FLA_Obj_is_real( A ) );

  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A    = ( float * ) A_vecs;
    float *buff_t    = ( float * ) t_scals;
    float *buff_B    = ( float * ) B_mat;

    rocsolver_sormbr( handle,
                      blas_vect,
                      blas_side,
                      blas_trans,
                      m_B,
                      n_B,
                      k_t,
                      buff_A,    cs_A,
                      buff_t,
                      buff_B,    cs_B );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A    = ( double * ) A_vecs;
    double *buff_t    = ( double * ) t_scals;
    double *buff_B    = ( double * ) B_mat;

    rocsolver_dormbr( handle,
                      blas_vect,
                      blas_side,
                      blas_trans,
                      m_B,
                      n_B,
                      k_t,
                      buff_A,    cs_A,
                      buff_t,
                      buff_B,    cs_B );

    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex *buff_A    = ( rocblas_float_complex * ) A_vecs;
    rocblas_float_complex *buff_t    = ( rocblas_float_complex * ) t_scals;
    rocblas_float_complex *buff_B    = ( rocblas_float_complex * ) B_mat;

    rocsolver_cunmbr( handle,
                      blas_vect,
                      blas_side,
                      blas_trans,
                      m_B,
                      n_B,
                      k_t,
                      buff_A,    cs_A,
                      buff_t,
                      buff_B,    cs_B );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex *buff_A    = ( rocblas_double_complex * ) A_vecs;
    rocblas_double_complex *buff_t    = ( rocblas_double_complex * ) t_scals;
    rocblas_double_complex *buff_B    = ( rocblas_double_complex * ) B_mat;

    rocsolver_zunmbr( handle,
                      blas_vect,
                      blas_side,
                      blas_trans,
                      m_B,
                      n_B,
                      k_t,
                      buff_A,    cs_A,
                      buff_t,
                      buff_B,    cs_B );

    break;
  }

  }

  if( conj_no_trans_a )
  {
    hipStream_t stream;
    rocblas_get_stream( handle, &stream );
    hipFreeAsync( stream, A_mat_corr );
  }

  return FLA_SUCCESS;
}

#endif
