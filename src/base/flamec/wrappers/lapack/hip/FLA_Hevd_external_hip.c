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
#include "rocblas.h"
#include "rocsolver.h"

FLA_Error FLA_Hevd_external_hip( rocblas_handle handle, FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj e, void* e_hip )
{
  FLA_Datatype datatype;
  int          n_A, cs_A;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Hevdd_check( jobz, uplo, A, e );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  rocblas_evect blas_jobz = FLA_Param_map_flame_to_rocblas_evd_type( jobz );
  rocblas_fill blas_uplo = FLA_Param_map_flame_to_rocblas_uplo( uplo );

  rocblas_int* info;
  hipMalloc( (void**) &info, sizeof( rocblas_int ) );
  void* buff_work;
  hipMalloc( &buff_work, n_A * FLA_Obj_datatype_size ( FLA_Obj_datatype( A ) ) );

  switch( datatype ) {

    case FLA_FLOAT:
    {
      float* buff_A     = ( float* ) FLA_FLOAT_PTR( A );
      float* buff_e     = ( float* ) FLA_FLOAT_PTR( e );

      rocsolver_ssyev( handle,
                       blas_jobz,
                       blas_uplo,
                       n_A,
                       buff_A,     cs_A,
                       buff_e,
                       buff_work,
                       info );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A     = ( double* ) FLA_DOUBLE_PTR( A );
      double* buff_e     = ( double* ) FLA_DOUBLE_PTR( e );
  
      rocsolver_dsyev( handle,
                       blas_jobz,
                       blas_uplo,
                       n_A,
                       buff_A,     cs_A,
                       buff_e,
                       buff_work,
                       info );
  
      break;
    } 
  
    case FLA_COMPLEX:
    {
      rocblas_float_complex* buff_A     = ( rocblas_float_complex* ) FLA_COMPLEX_PTR( A );
      float*    buff_e     = ( float*    ) FLA_FLOAT_PTR( e );
  
      rocsolver_cheev( handle,
                       blas_jobz,
                       blas_uplo,
                       n_A,
                       buff_A,     cs_A,
                       buff_e,
                       buff_work,
                       info );
  
      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      rocblas_double_complex* buff_A     = ( rocblas_double_complex* ) FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_e     = ( double*   ) FLA_DOUBLE_PTR( e );
  
      rocsolver_zheev( handle,
                       blas_jobz,
                       blas_uplo,
                       n_A,
                       buff_A,     cs_A,
                       buff_e,
                       buff_work,
                       info );
  
      break;
    } 

  }

  int rval = ( info == hipSuccess ) ? FLA_SUCCESS : FLA_FAILURE;

  hipFree( info );
  hipFree( buff_work );

  return rval;
}

#endif
