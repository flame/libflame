/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Max_abs_value_herm( FLA_Uplo uplo, FLA_Obj A, FLA_Obj maxabs )
{
  FLA_Datatype datatype;
  dim_t        n_A;
  dim_t        rs_A, cs_A;
  uplo1_t       blis_uplo;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Max_abs_value_herm_check( uplo, A, maxabs );

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );
 
  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );
 
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float*    buff_A      = ( float    * ) FLA_FLOAT_PTR( A );
    float*    buff_maxabs = ( float    * ) FLA_FLOAT_PTR( maxabs );

    bl1_smaxabsmr( blis_uplo,
                   n_A,
                   n_A,
                   buff_A, rs_A, cs_A,
                   buff_maxabs );

    break;
  }

  case FLA_DOUBLE:
  {
    double*   buff_A      = ( double   * ) FLA_DOUBLE_PTR( A );
    double*   buff_maxabs = ( double   * ) FLA_DOUBLE_PTR( maxabs );

    bl1_dmaxabsmr( blis_uplo,
                   n_A,
                   n_A,
                   buff_A, rs_A, cs_A,
                   buff_maxabs );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A      = ( scomplex * ) FLA_COMPLEX_PTR( A );
    float    *buff_maxabs = ( float    * ) FLA_FLOAT_PTR( maxabs );

    bl1_cmaxabsmr( blis_uplo,
                   n_A,
                   n_A,
                   buff_A, rs_A, cs_A,
                   buff_maxabs );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    double   *buff_maxabs = ( double   * ) FLA_DOUBLE_PTR( maxabs );

    bl1_zmaxabsmr( blis_uplo,
                   n_A,
                   n_A,
                   buff_A, rs_A, cs_A,
                   buff_maxabs );

    break;
  }

  }

  return FLA_SUCCESS;
}

