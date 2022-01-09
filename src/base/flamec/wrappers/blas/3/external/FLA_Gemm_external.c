/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Gemm_external( FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Datatype datatype;
  integer          k_AB;
  integer          m_A, n_A;
  integer          m_C, n_C;
  integer          rs_A, cs_A;
  integer          rs_B, cs_B;
  integer          rs_C, cs_C;
  trans1_t      blis_transa;
  trans1_t      blis_transb;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Gemm_check( transa, transb, alpha, A, B, beta, C );

  if ( FLA_Obj_has_zero_dim( C ) ) return FLA_SUCCESS;

  if ( FLA_Obj_has_zero_dim( A ) || FLA_Obj_has_zero_dim( B ) )
  {
    FLA_Scal_external( beta, C );
    return FLA_SUCCESS;
  }

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );

  m_C      = FLA_Obj_length( C );
  n_C      = FLA_Obj_width( C );
  rs_C     = FLA_Obj_row_stride( C );
  cs_C     = FLA_Obj_col_stride( C );

  if ( transa == FLA_NO_TRANSPOSE || transa == FLA_CONJ_NO_TRANSPOSE )
    k_AB = n_A;
  else
    k_AB = m_A;

  FLA_Param_map_flame_to_blis_trans( transa, &blis_transa );
  FLA_Param_map_flame_to_blis_trans( transb, &blis_transb );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_B     = ( float * ) FLA_FLOAT_PTR( B );
    float *buff_C     = ( float * ) FLA_FLOAT_PTR( C );
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    bl1_sgemm( blis_transa, 
               blis_transb,
               m_C,
               k_AB, 
               n_C, 
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_B, rs_B, cs_B,
               buff_beta,
               buff_C, rs_C, cs_C );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_B     = ( double * ) FLA_DOUBLE_PTR( B );
    double *buff_C     = ( double * ) FLA_DOUBLE_PTR( C );
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );
                         
    bl1_dgemm( blis_transa, 
               blis_transb,
               m_C,
               k_AB, 
               n_C, 
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_B, rs_B, cs_B,
               buff_beta,
               buff_C, rs_C, cs_C );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_B     = ( scomplex * ) FLA_COMPLEX_PTR( B );
    scomplex *buff_C     = ( scomplex * ) FLA_COMPLEX_PTR( C );
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );
    scomplex *buff_beta  = ( scomplex * ) FLA_COMPLEX_PTR( beta );
                         
    bl1_cgemm( blis_transa, 
               blis_transb,
               m_C,
               k_AB, 
               n_C, 
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_B, rs_B, cs_B,
               buff_beta,
               buff_C, rs_C, cs_C );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_B     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
    dcomplex *buff_C     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( C );
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    dcomplex *buff_beta  = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    bl1_zgemm( blis_transa, 
               blis_transb,
               m_C,
               k_AB, 
               n_C, 
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_B, rs_B, cs_B,
               buff_beta,
               buff_C, rs_C, cs_C );

    break;
  }

  }

  return FLA_SUCCESS;
}

