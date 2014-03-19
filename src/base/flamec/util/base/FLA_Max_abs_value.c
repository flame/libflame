/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Max_abs_value( FLA_Obj A, FLA_Obj maxabs )
{
  FLA_Datatype datatype;
  FLA_Datatype dt_maxabs;
  dim_t        m_A, n_A;
  dim_t        rs_A, cs_A;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Max_abs_value_check( A, maxabs );

  datatype  = FLA_Obj_datatype( A );
  dt_maxabs = FLA_Obj_datatype( maxabs );

  m_A       = FLA_Obj_length( A );
  n_A       = FLA_Obj_width( A );
  rs_A      = FLA_Obj_row_stride( A );
  cs_A      = FLA_Obj_col_stride( A );
 
 
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float*    buff_A      = ( float    * ) FLA_FLOAT_PTR( A );
    float*    buff_maxabs = ( float    * ) FLA_FLOAT_PTR( maxabs );

    bl1_smaxabsm( m_A,
                  n_A,
                  buff_A, rs_A, cs_A,
                  buff_maxabs );

    break;
  }

  case FLA_DOUBLE:
  {
    double*   buff_A      = ( double   * ) FLA_DOUBLE_PTR( A );
    double*   buff_maxabs = ( double   * ) FLA_DOUBLE_PTR( maxabs );

    bl1_dmaxabsm( m_A,
                  n_A,
                  buff_A, rs_A, cs_A,
                  buff_maxabs );

    break;
  }

  case FLA_COMPLEX:
  {
    if ( dt_maxabs == FLA_FLOAT )
    {
      scomplex* buff_A      = ( scomplex * ) FLA_COMPLEX_PTR( A );
      float*    buff_maxabs = ( float    * ) FLA_FLOAT_PTR( maxabs );

      bl1_cmaxabsm( m_A,
                    n_A,
                    buff_A, rs_A, cs_A,
                    buff_maxabs );
    }
    else
    {
      scomplex* buff_A      = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex* buff_maxabs = ( scomplex * ) FLA_COMPLEX_PTR( maxabs );

      bl1_cmaxabsm( m_A,
                    n_A,
                    buff_A, rs_A, cs_A,
                    &(buff_maxabs->real) );

      buff_maxabs->imag = 0.0;
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    if ( dt_maxabs == FLA_DOUBLE )
    {
      dcomplex* buff_A      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_maxabs = ( double   * ) FLA_DOUBLE_PTR( maxabs );

      bl1_zmaxabsm( m_A,
                    n_A,
                    buff_A, rs_A, cs_A,
                    buff_maxabs );
    }
    else
    {
      dcomplex* buff_A      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_maxabs = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( maxabs );

      bl1_zmaxabsm( m_A,
                    n_A,
                    buff_A, rs_A, cs_A,
                    &(buff_maxabs->real) );

      buff_maxabs->imag = 0.0;
    }

    break;
  }

  }

  return FLA_SUCCESS;
}

