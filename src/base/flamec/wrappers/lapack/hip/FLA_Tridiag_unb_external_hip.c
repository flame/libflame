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

FLA_Error FLA_Tridiag_unb_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip )
{
  FLA_Datatype datatype;
  int          n_A, cs_A;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Tridiag_check( uplo, A, t );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  void* buff_d;
  void* buff_e;
  hipStream_t stream;
  rocblas_get_stream( handle, &stream );
  hipMallocAsync( &buff_d, FLA_Obj_datatype_size ( FLA_Obj_datatype_proj_to_real( A ) ) * n_A, stream );
  hipMallocAsync( &buff_e, FLA_Obj_datatype_size ( FLA_Obj_datatype_proj_to_real( A ) ) * ( n_A - 1 ), stream );

  rocblas_fill blas_uplo = FLA_Param_map_flame_to_rocblas_uplo( uplo );

  void* A_mat = NULL;
  void* t_vec = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
    t_vec = FLA_Obj_buffer_at_view( t );
  }
  else
  {
    A_mat = A_hip;
    t_vec = t_hip;
  }

  switch( datatype ){

  case FLA_FLOAT:
  {
    float* buff_A    = ( float * ) A_mat;
    float* buff_t    = ( float * ) t_vec;

    rocsolver_ssytd2( handle,
                      blas_uplo,
                      n_A,
                      buff_A, cs_A,
                      buff_d,
                      buff_e,
                      buff_t );

    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_A    = ( double * ) A_mat;
    double* buff_t    = ( double * ) t_vec;

    rocsolver_dsytd2( handle,
                      blas_uplo,
                      n_A,
                      buff_A, cs_A,
                      buff_d,
                      buff_e,
                      buff_t );

    break;
  } 

  case FLA_COMPLEX:
  {
    rocblas_float_complex* buff_A    = A_mat;
    rocblas_float_complex* buff_t    = t_vec;

    rocsolver_chetd2( handle,
                      blas_uplo,
                      n_A,
                      buff_A, cs_A,
                      buff_d,
                      buff_e,
                      buff_t );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex* buff_A    = A_mat;
    rocblas_double_complex* buff_t    = t_vec;

    rocsolver_zhetd2( handle,
                      blas_uplo,
                      n_A,
                      buff_A, cs_A,
                      buff_d,
                      buff_e,
                      buff_t );

    break;
  } 

  }

  hipFreeAsync( stream, buff_d );
  hipFreeAsync( stream, buff_e );

  return FLA_SUCCESS;
}

FLA_Error FLA_Tridiag_unb_ext_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip )
{
  return FLA_Tridiag_unb_external_hip( handle, uplo, A, A_hip, t, t_hip );
}
#endif
