/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Conjugate( FLA_Obj A )
{
  FLA_Datatype datatype;
  integer          m_A, n_A;
  integer          rs_A, cs_A;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING ) 
    FLA_Conjugate_check( A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  if ( FLA_Obj_is_real( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );


  switch ( datatype ){

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );

    bl1_cconjm( m_A,
                n_A,
                buff_A, rs_A, cs_A );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  { 
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

    bl1_zconjm( m_A,
                n_A,
                buff_A, rs_A, cs_A );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

