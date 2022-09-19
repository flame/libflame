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

FLA_Error FLA_Bsvd_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj d, void* d_hip, FLA_Obj e, void* e_hip, FLA_Obj U, void* U_hip, FLA_Obj V, void* V_hip )
{
  FLA_Datatype datatype;
  int          m_U, cs_U;
  int          n_V, cs_V;
  int          n_C, cs_C;
  int          min_m_n;

  //if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
  //  FLA_Hevd_check( jobz, uplo, A, e );

  if ( FLA_Obj_has_zero_dim( d ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( U );

  m_U      = FLA_Obj_length( U );
  cs_U     = FLA_Obj_col_stride( U );

  n_V      = FLA_Obj_length( V );
  cs_V     = FLA_Obj_col_stride( V );

  n_C      = 0;
  cs_C     = 1;

  min_m_n  = FLA_Obj_vector_dim( d );

  rocblas_fill blas_uplo = FLA_Param_map_flame_to_rocblas_uplo( uplo );

  void* d_vec = NULL;
  void* e_vec = NULL;
  void* U_mat = NULL;
  void* V_mat = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    d_vec = FLA_Obj_buffer_at_view( d );
    e_vec = FLA_Obj_buffer_at_view( e );
    U_mat = FLA_Obj_buffer_at_view( U );
    V_mat = FLA_Obj_buffer_at_view( V );
  }
  else
  {
    d_vec = d_hip;
    e_vec = e_hip;
    U_mat = U_hip;
    V_mat = V_hip;
  }

  rocblas_int* info;
  hipMalloc( (void**) &info, sizeof( rocblas_int ) );

  switch( datatype ) {

    case FLA_FLOAT:
    {
      float*    buff_d     = ( float * ) d_vec;
      float*    buff_e     = ( float * ) e_vec;
      float*    buff_U     = ( float * ) U_mat;
      float*    buff_V     = ( float * ) V_mat;
      float*    buff_C     = ( float * ) NULL;
  
      rocsolver_sbdsqr( handle,
                        blas_uplo,
                        min_m_n,
                        n_V,
                        m_U,
                        n_C,
                        buff_d,
                        buff_e,
                        buff_V, cs_V,
                        buff_U, cs_U,
                        buff_C, cs_C,
                        info );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_d     = ( double * ) FLA_DOUBLE_PTR( d );
      double*   buff_e     = ( double * ) FLA_DOUBLE_PTR( e );
      double*   buff_U     = ( double * ) FLA_DOUBLE_PTR( U );
      double*   buff_V     = ( double * ) FLA_DOUBLE_PTR( V );
      double*   buff_C     = ( double * ) NULL;
  
      rocsolver_dbdsqr( handle,
                        blas_uplo,
                        min_m_n,
                        n_V,
                        m_U,
                        n_C,
                        buff_d,
                        buff_e,
                        buff_V, cs_V,
                        buff_U, cs_U,
                        buff_C, cs_C,
                        info );

      break;
    } 
  
    case FLA_COMPLEX:
    {
      float*    buff_d     = ( float    * ) FLA_FLOAT_PTR( d );
      float*    buff_e     = ( float    * ) FLA_FLOAT_PTR( e );
      rocblas_float_complex* buff_U     = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( U );
      rocblas_float_complex* buff_V     = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( V );
      rocblas_float_complex* buff_C     = ( rocblas_float_complex * ) NULL;
  
      rocsolver_cbdsqr( handle,
                        blas_uplo,
                        min_m_n,
                        n_V,
                        m_U,
                        n_C,
                        buff_d,
                        buff_e,
                        buff_V, cs_V,
                        buff_U, cs_U,
                        buff_C, cs_C,
                        info );

      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      double*   buff_d     = ( double   * ) FLA_DOUBLE_PTR( d );
      double*   buff_e     = ( double   * ) FLA_DOUBLE_PTR( e );
      rocblas_double_complex* buff_U     = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( U );
      rocblas_double_complex* buff_V     = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( V );
      rocblas_double_complex* buff_C     = ( rocblas_double_complex * ) NULL;
  
      rocsolver_zbdsqr( handle,
                        blas_uplo,
                        min_m_n,
                        n_V,
                        m_U,
                        n_C,
                        buff_d,
                        buff_e,
                        buff_V, cs_V,
                        buff_U, cs_U,
                        buff_C, cs_C,
                        info );

      break;
    } 

  }

  int rval = ( info == hipSuccess ) ? FLA_SUCCESS : FLA_FAILURE;

  hipFree(info);

  return rval;
}

#endif
