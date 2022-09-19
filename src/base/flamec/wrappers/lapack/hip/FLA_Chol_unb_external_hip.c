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

FLA_Error FLA_Chol_l_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip )
{
  return FLA_Chol_unb_external_hip( handle, FLA_LOWER_TRIANGULAR, A, A_hip ); 
}

FLA_Error FLA_Chol_u_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip )
{
  return FLA_Chol_unb_external_hip( handle, FLA_UPPER_TRIANGULAR, A, A_hip );
}

FLA_Error FLA_Chol_unb_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip )
{
  FLA_Datatype datatype;
  int          n_A;
  int          ldim_A;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Chol_check( uplo, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  ldim_A   = FLA_Obj_col_stride( A );

  rocblas_fill blas_uplo = FLA_Param_map_flame_to_rocblas_uplo( uplo );
  rocblas_int* info;
  hipMalloc( (void**) &info, sizeof( rocblas_int ) );

  void* A_mat = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
  }
  else
  {
    A_mat = A_hip;
  }

  switch( datatype ){

  case FLA_FLOAT:
  {
    rocsolver_spotf2( handle,
                      blas_uplo,
                      n_A,
                      ( float * ) A_mat, ldim_A,
                      info );
    
    break;
  }

  case FLA_DOUBLE:
  {
    rocsolver_dpotf2( handle,
                      blas_uplo,
                      n_A,
                      ( double * ) A_mat, ldim_A,
                      info );

    break;
  }

  case FLA_COMPLEX:
  {
    rocsolver_cpotf2( handle,
                      blas_uplo,
                      n_A,
                      ( rocblas_float_complex * ) A_mat, ldim_A,
                      info );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocsolver_zpotf2( handle,
                      blas_uplo,
                      n_A,
                      ( rocblas_double_complex * ) A_mat, ldim_A,
                      info );

    break;
  }

  }
  hipFree(info);

  return FLA_SUCCESS;
}

#endif
