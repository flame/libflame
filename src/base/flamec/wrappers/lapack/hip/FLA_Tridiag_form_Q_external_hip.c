/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include "hip/hip_runtime_api.h"
#include "rocblas/rocblas.h"
#include "rocsolver/rocsolver.h"

FLA_Error FLA_Tridiag_form_Q_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip )
{
  FLA_Datatype datatype;
  int          m_A;
  int          cs_A;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Tridiag_form_Q_check( uplo, A, t );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  cs_A     = FLA_Obj_col_stride( A );

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
    float*    buff_A    = ( float    * ) A_mat;
    float*    buff_t    = ( float    * ) t_vec;

    rocsolver_sorgtr( handle,
                      blas_uplo,
                      m_A,
                      buff_A, cs_A,
                      buff_t );

    break;
  }

  case FLA_DOUBLE:
  {
    double*   buff_A    = ( double   * ) A_mat;
    double*   buff_t    = ( double   * ) t_vec;

    rocsolver_dorgtr( handle,
                      blas_uplo,
                      m_A,
                      buff_A, cs_A,
                      buff_t );

    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex* buff_A    = ( rocblas_float_complex * ) A_mat;
    rocblas_float_complex* buff_t    = ( rocblas_float_complex * ) t_vec;

    rocsolver_cungtr( handle,
                      blas_uplo,
                      m_A,
                      buff_A, cs_A,
                      buff_t );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex *buff_A    = ( rocblas_double_complex * ) A_mat;
    rocblas_double_complex *buff_t    = ( rocblas_double_complex * ) t_vec;

    rocsolver_zungtr( handle,
                      blas_uplo,
                      m_A,
                      buff_A, cs_A,
                      buff_t );

    break;
  }

  }

  return FLA_SUCCESS;
}
#endif
