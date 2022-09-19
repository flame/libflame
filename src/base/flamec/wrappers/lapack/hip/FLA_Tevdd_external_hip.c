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

FLA_Error FLA_Tevdd_external_hip( rocblas_handle handle, FLA_Evd_type jobz, FLA_Obj d, void* d_hip, FLA_Obj e, void* e_hip, FLA_Obj A, void* A_hip )
{
  FLA_Datatype datatype;
  int          n_A, cs_A;

  //if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
  //  FLA_Tevdd_check( jobz, d, e, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  rocblas_evect blas_jobz = FLA_Param_map_flame_to_rocblas_evd_type( jobz );

  rocblas_int* info;
  hipMalloc( (void**) &info, sizeof( rocblas_int ) );

  void* d_vec = NULL;
  void* e_vec = NULL;
  void* A_mat = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    d_vec = FLA_Obj_buffer_at_view( d );
    e_vec = FLA_Obj_buffer_at_view( e );
    A_mat = FLA_Obj_buffer_at_view( A );
  }
  else
  {
    d_vec = d_hip;
    e_vec = e_hip;
    A_mat = A_hip;
  }

  switch( datatype ) {

    case FLA_FLOAT:
    {
      float* buff_d     = ( float* ) d_vec;
      float* buff_e     = ( float* ) e_vec;
      float* buff_A     = ( float* ) A_mat;

      rocsolver_sstedc( handle,
                        blas_jobz,
                        n_A,
                        buff_d,
                        buff_e,
                        buff_A,     cs_A,
                        info );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_d     = ( double* ) d_vec;
      double* buff_e     = ( double* ) e_vec;
      double* buff_A     = ( double* ) A_mat;
  
      rocsolver_dstedc( handle,
                        blas_jobz,
                        n_A,
                        buff_d,
                        buff_e,
                        buff_A,     cs_A,
                        info );
  
      break;
    } 
  
    case FLA_COMPLEX:
    {
      float*    buff_d     = ( float*    ) d_vec;
      float*    buff_e     = ( float*    ) e_vec;
      rocblas_float_complex* buff_A     = ( rocblas_float_complex* ) A_mat;
  
      rocsolver_cstedc( handle,
                        blas_jobz,
                        n_A,
                        buff_d,
                        buff_e,
                        buff_A,     cs_A,
                        info );
  
      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      double*   buff_d     = ( double*   ) d_vec;
      double*   buff_e     = ( double*   ) e_vec;
      rocblas_double_complex* buff_A     = ( rocblas_double_complex* ) A_mat;
  
      rocsolver_zstedc( handle,
                        blas_jobz,
                        n_A,
                        buff_d,
                        buff_e,
                        buff_A,     cs_A,
                        info );
  
      break;
    } 
  }

  int rval = ( info == hipSuccess ) ? FLA_SUCCESS : FLA_FAILURE;

  hipFree( info );

  return rval;
}

#endif
