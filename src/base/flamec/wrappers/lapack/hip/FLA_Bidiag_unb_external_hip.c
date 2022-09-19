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

FLA_Error FLA_Bidiag_unb_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj tu, void* tu_hip, FLA_Obj tv, void* tv_hip )
{
  int          info = 0;
  FLA_Datatype datatype;
  int          m_A, n_A, cs_A;
  int          min_m_n;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Bidiag_check( A, tu, tv );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  min_m_n  = FLA_Obj_min_dim( A );
  cs_A     = FLA_Obj_col_stride( A );

  void* buff_d;
  void* buff_e;
  hipMalloc( &buff_d, FLA_Obj_datatype_size ( FLA_Obj_datatype_proj_to_real( A ) ) * min_m_n );
  hipMalloc( &buff_e, FLA_Obj_datatype_size ( FLA_Obj_datatype_proj_to_real( A ) ) * min_m_n );

  void* A_mat = NULL;
  void* tu_scals = NULL;
  void* tv_scals = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
    tu_scals = FLA_Obj_buffer_at_view( tu );
    tv_scals = FLA_Obj_buffer_at_view( tv );
  }
  else
  {
    A_mat = A_hip;
    tu_scals = tu_hip;
    tv_scals = tv_hip;
  }


  switch( datatype ){

  case FLA_FLOAT:
  {
    float* buff_A    = ( float * ) A_mat;
    float* buff_tu   = ( float * ) tu_scals;
    float* buff_tv   = ( float * ) tv_scals;

    rocsolver_sgebd2( handle,
                      m_A,
                      n_A,
                      buff_A, cs_A,
                      buff_d,
                      buff_e,
                      buff_tu,
                      buff_tv );

    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_A    = ( double * ) A_mat;
    double* buff_tu   = ( double * ) tu_scals;
    double* buff_tv   = ( double * ) FLA_DOUBLE_PTR( tv );

    rocsolver_dgebd2( handle,
                      m_A,
                      n_A,
                      buff_A, cs_A,
                      buff_d,
                      buff_e,
                      buff_tu,
                      buff_tv );

    break;
  } 

  case FLA_COMPLEX:
  {
    rocblas_float_complex* buff_A    = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( A );
    rocblas_float_complex* buff_tu   = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( tu );
    rocblas_float_complex* buff_tv   = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( tv );

    rocsolver_cgebd2( handle,
                      m_A,
                      n_A,
                      buff_A, cs_A,
                      buff_d,
                      buff_e,
                      buff_tu,
                      buff_tv );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex* buff_A    = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    rocblas_double_complex* buff_tu   = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( tu );
    rocblas_double_complex* buff_tv   = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( tv );

    rocsolver_zgebd2( handle,
                      m_A,
                      n_A,
                      buff_A, cs_A,
                      buff_d,
                      buff_e,
                      buff_tu,
                      buff_tv );

    break;
  } 

  }

  hipFree( buff_d );
  hipFree( buff_e );

  return info;
}

FLA_Error FLA_Bidiag_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj tu, void* tu_hip, FLA_Obj tv, void* tv_hip )
{
  return FLA_Bidiag_unb_external_hip( handle, A, A_hip, tu, tu_hip, tv, tv_hip );
}

#endif
