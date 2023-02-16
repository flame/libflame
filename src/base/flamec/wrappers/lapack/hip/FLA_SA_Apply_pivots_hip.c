/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2023, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include "hip/hip_runtime_api.h"
#include "rocblas/rocblas.h"

FLA_Error FLA_SA_Apply_pivots_hip( rocblas_handle handle, FLA_Obj C, void* C_hip, FLA_Obj E, void* E_hip, FLA_Obj p )
{
  FLA_Datatype datatype;
  int          m_C, n_C, cs_C;
  int                    cs_E;
  // int                    rs_C;
  // int                    rs_E;
  int          m_p;
  int          i;
  int*         buff_p;

  if ( FLA_Obj_has_zero_dim( C ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( C );

  m_C    = FLA_Obj_length( C );
  n_C    = FLA_Obj_width( C );
  cs_C   = FLA_Obj_col_stride( C );
  // rs_C   = FLA_Obj_row_stride( C );

  cs_E   = FLA_Obj_col_stride( E );
  // rs_E   = FLA_Obj_row_stride( E );

  m_p    = FLA_Obj_length( p );
  
  buff_p = ( int * ) FLA_INT_PTR( p );

  void* C_mat = NULL;
  void* E_mat = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    C_mat = FLA_Obj_buffer_at_view( C );
    E_mat = FLA_Obj_buffer_at_view( E );
  }
  else
  {
    C_mat = C_hip;
    E_mat = E_hip;
  }

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_C = ( float * ) C_mat;
    float* buff_E = ( float * ) E_mat;

    for ( i = 0; i < m_p; ++i )
    {
      if ( buff_p[ i ] != 0 ) 
        rocblas_sswap( handle, n_C,
                       buff_C + 0*cs_C + i,                         cs_C, 
                       buff_E + 0*cs_E + buff_p[ i ] - ( m_C - i ), cs_E );
    }
    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_C = ( double * ) C_mat;
    double* buff_E = ( double * ) E_mat;

    for ( i = 0; i < m_p; ++i )
    {
      if ( buff_p[ i ] != 0 ) 
        rocblas_dswap( handle, n_C, 
                       buff_C + 0*cs_C + i,                         cs_C, 
                       buff_E + 0*cs_E + buff_p[ i ] - ( m_C - i ), cs_E );
    }
    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex* buff_C = ( rocblas_float_complex * ) C_mat;
    rocblas_float_complex* buff_E = ( rocblas_float_complex * ) E_mat;

    for ( i = 0; i < m_p; ++i )
    {
      if ( buff_p[ i ] != 0 ) 
        rocblas_cswap( handle, n_C,
                       buff_C + 0*cs_C + i,                         cs_C, 
                       buff_E + 0*cs_E + buff_p[ i ] - ( m_C - i ), cs_E );
    }
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex* buff_C = ( rocblas_double_complex * ) C_mat;
    rocblas_double_complex* buff_E = ( rocblas_double_complex * ) E_mat;

    for ( i = 0; i < m_p; ++i )
    {
      if ( buff_p[ i ] != 0 ) 
        rocblas_zswap( handle, n_C,
                       buff_C + 0*cs_C + i,                         cs_C, 
                       buff_E + 0*cs_E + buff_p[ i ] - ( m_C - i ), cs_E );
    }
    break;
  }

  }

  return FLA_SUCCESS;
}

#endif
