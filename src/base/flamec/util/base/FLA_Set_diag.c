/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Set_diag( FLA_Obj alpha, FLA_Obj A )
{
  FLA_Datatype datatype;
  integer          m_A, n_A;
  integer          rs_A, cs_A;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Set_diag_check( alpha, A );

  datatype = FLA_Obj_datatype( A );
  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  switch ( datatype ){

  case FLA_INT:
  {
    integer *buff_A     = ( integer * ) FLA_INT_PTR( A );
    integer *buff_alpha = ( integer * ) FLA_INT_PTR( alpha );

    bl1_isetdiag( 0,
                  m_A,
                  n_A,
                  buff_alpha,
                  buff_A, rs_A, cs_A );

    break;
  }

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    bl1_ssetdiag( 0,
                  m_A,
                  n_A,
                  buff_alpha,
                  buff_A, rs_A, cs_A );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    bl1_dsetdiag( 0,
                  m_A,
                  n_A,
                  buff_alpha,
                  buff_A, rs_A, cs_A );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );

    bl1_csetdiag( 0,
                  m_A,
                  n_A,
                  buff_alpha,
                  buff_A, rs_A, cs_A );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );

    bl1_zsetdiag( 0,
                  m_A,
                  n_A,
                  buff_alpha,
                  buff_A, rs_A, cs_A );

    break;
  }

  }

  return FLA_SUCCESS;
}
