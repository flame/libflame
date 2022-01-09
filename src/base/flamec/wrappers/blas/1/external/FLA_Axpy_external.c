/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Axpy_external( FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype datatype;
  integer          m_B, n_B;
  integer          rs_A, cs_A;
  integer          rs_B, cs_B;
  trans1_t      blis_trans;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Axpy_check( alpha, A, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );

  if ( FLA_Obj_is_conformal_to( FLA_NO_TRANSPOSE, A, B ) )
    FLA_Param_map_flame_to_blis_trans( FLA_NO_TRANSPOSE, &blis_trans );
  else // if ( FLA_Obj_is_conformal_to( FLA_TRANSPOSE, A, B ) )
    FLA_Param_map_flame_to_blis_trans( FLA_TRANSPOSE, &blis_trans );

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_B     = ( float * ) FLA_FLOAT_PTR( B );

    bl1_saxpymt( blis_trans,
                 m_B,
                 n_B,
                 buff_alpha,
                 buff_A, rs_A, cs_A,
                 buff_B, rs_B, cs_B );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_B     = ( double * ) FLA_DOUBLE_PTR( B );

    bl1_daxpymt( blis_trans,
                 m_B,
                 n_B,
                 buff_alpha,
                 buff_A, rs_A, cs_A,
                 buff_B, rs_B, cs_B );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );
    scomplex *buff_A =     ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_B =     ( scomplex * ) FLA_COMPLEX_PTR( B );

    bl1_caxpymt( blis_trans,
                 m_B,
                 n_B,
                 buff_alpha,
                 buff_A, rs_A, cs_A,
                 buff_B, rs_B, cs_B );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_B     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

    bl1_zaxpymt( blis_trans,
                 m_B,
                 n_B,
                 buff_alpha,
                 buff_A, rs_A, cs_A,
                 buff_B, rs_B, cs_B );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

