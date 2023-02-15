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
#include "rocsolver/rocsolver.h"

FLA_Error FLA_Apply_pivots_unb_external_hip( rocblas_handle handle, FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, void* A_hip )
{
  FLA_Datatype datatype;
  int          n_A, cs_A;
  int          m_p;
  int          inc_p;
  int*         buff_p;
  int          k1_1, k2_1;
  int*         pivots_lapack;
  int          i;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Apply_pivots_check( side, trans, p, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_p    = FLA_Obj_vector_inc( p );
  m_p      = FLA_Obj_vector_dim( p );

  buff_p   = FLA_INT_PTR( p );

  // Use one-based indices for LAPACK.
  k1_1     = 1;
  k2_1     = m_p;

  // Translate FLAME pivot indices to LAPACK-compatible indices. It is
  // important to note that this conversion, unlike the one done by
  // FLA_Shift_pivots_to(), is NOT in-place, but rather done separately
  // in a temporary buffer.
  hipMallocManaged( (void**) &pivots_lapack, m_p * sizeof( int ), hipMemAttachGlobal );
  // this implies a sync

  for ( i = 0; i < m_p; i++ )
  {
    pivots_lapack[ i ] = buff_p[ i ] + i + 1;
  }

  void* A_mat = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
  }
  else
  {
    A_mat = A_hip;
  }

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_A = ( float * ) A_mat;

    rocsolver_slaswp( handle,
                      n_A,
                      buff_A, cs_A,
                      k1_1, 
                      k2_1,
                      pivots_lapack,
                      inc_p );
    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_A = ( double * ) A_mat;

    rocsolver_dlaswp( handle,
                      n_A,
                      buff_A, cs_A,
                      k1_1, 
                      k2_1,
                      pivots_lapack,
                      inc_p );
    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex* buff_A = ( rocblas_float_complex * ) A_mat;

    rocsolver_claswp( handle,
                      n_A,
                      buff_A, cs_A,
                      k1_1, 
                      k2_1,
                      pivots_lapack,
                      inc_p );
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex* buff_A = ( rocblas_double_complex * ) A_mat; 

    rocsolver_zlaswp( handle,
                      n_A,
                      buff_A, cs_A,
                      k1_1, 
                      k2_1,
                      pivots_lapack,
                      inc_p );
    break;
  }

  }

  hipFree( pivots_lapack );

  return FLA_SUCCESS;
}

FLA_Error FLA_Apply_pivots_ln_unb_ext_hip( rocblas_handle handle, FLA_Obj p, FLA_Obj A, void* A_hip )
{
  return FLA_Apply_pivots_unb_external_hip( handle, FLA_LEFT, FLA_NO_TRANSPOSE, p, A, A_hip );
}

#endif
