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

FLA_Error FLA_Svd_external_hip( rocblas_handle handle, FLA_Svd_type jobu, FLA_Svd_type jobv, FLA_Obj A, void* A_hip, FLA_Obj s, void* s_hip, FLA_Obj U, void* U_hip, FLA_Obj V, void* V_hip )
{
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  int          m_A, n_A, cs_A;
  int          cs_U;
  int          cs_V;
  int          min_m_n;
  rocblas_workmode fast_alg = rocblas_outofplace;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Svd_check( jobu, jobv, A, s, U, V );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );
  dt_real  = FLA_Obj_datatype_proj_to_real( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  cs_U     = FLA_Obj_col_stride( U );

  cs_V     = FLA_Obj_col_stride( V );

  min_m_n  = min( m_A, n_A );

  rocblas_svect blas_jobu = FLA_Param_map_flame_to_rocblas_svd_type( jobu );
  rocblas_svect blas_jobv = FLA_Param_map_flame_to_rocblas_svd_type( jobv );

  rocblas_int* info;
  hipMalloc( (void**) &info, sizeof( rocblas_int ) );
  void* buff_work;
  hipMalloc( &buff_work, FLA_Obj_datatype_size ( dt_real ) * min_m_n );

  void* A_mat = NULL;
  void* s_vec = NULL;
  void* U_mat = NULL;
  void* V_mat = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
    s_vec = FLA_Obj_buffer_at_view( s );
    U_mat = FLA_Obj_buffer_at_view( U );
    V_mat = FLA_Obj_buffer_at_view( V );
  }
  else
  {
    A_mat = A_hip;
    s_vec = s_hip;
    U_mat = U_hip;
    V_mat = V_hip;
  }

  switch( datatype ) {

    case FLA_FLOAT:
    {
      float*    buff_A     = ( float*    ) A_mat;
      float*    buff_s     = ( float*    ) s_vec;
      float*    buff_U     = ( float*    ) U_mat;
      float*    buff_V     = ( float*    ) V_mat;
  
      rocsolver_sgesvd( handle,
                        blas_jobu,
                        blas_jobv,
                        m_A,
                        n_A,
                        buff_A,    cs_A,
                        buff_s,
                        buff_U,    cs_U,
                        buff_V,    cs_V,
                        buff_work,
                        fast_alg,
                        info );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A     = ( double*   ) A_mat;
      double*   buff_s     = ( double*   ) s_vec;
      double*   buff_U     = ( double*   ) U_mat;
      double*   buff_V     = ( double*   ) V_mat;
  
      rocsolver_dgesvd( handle,
                        blas_jobu,
                        blas_jobv,
                        m_A,
                        n_A,
                        buff_A,    cs_A,
                        buff_s,
                        buff_U,    cs_U,
                        buff_V,    cs_V,
                        buff_work,
                        fast_alg,
                        info );
  
      break;
    } 
  
    case FLA_COMPLEX:
    {
      rocblas_float_complex* buff_A     = ( rocblas_float_complex* ) A_mat;
      float*    buff_s     = ( float*    ) s_vec;
      rocblas_float_complex* buff_U     = ( rocblas_float_complex* ) U_mat;
      rocblas_float_complex* buff_V     = ( rocblas_float_complex* ) V_mat;
  
      rocsolver_cgesvd( handle,
                        blas_jobu,
                        blas_jobv,
                        m_A,
                        n_A,
                        buff_A,    cs_A,
                        buff_s,
                        buff_U,    cs_U,
                        buff_V,    cs_V,
                        buff_work,
                        fast_alg,
                        info );
  
      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      rocblas_double_complex* buff_A     = ( rocblas_double_complex* ) A_mat;
      double*   buff_s     = ( double*   ) s_vec;
      rocblas_double_complex* buff_U     = ( rocblas_double_complex* ) U_mat;
      rocblas_double_complex* buff_V     = ( rocblas_double_complex* ) V_mat;
  
      rocsolver_zgesvd( handle,
                        blas_jobu,
                        blas_jobv,
                        m_A,
                        n_A,
                        buff_A,    cs_A,
                        buff_s,
                        buff_U,    cs_U,
                        buff_V,    cs_V,
                        buff_work,
                        fast_alg,
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

