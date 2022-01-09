/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Trinv_lu_opt_var3( FLA_Obj A )
{
  FLA_Datatype datatype;
  integer          mn_A;
  integer          rs_A, cs_A;

  datatype = FLA_Obj_datatype( A );

  mn_A     = FLA_Obj_length( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );

      FLA_Trinv_lu_ops_var3( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Trinv_lu_opd_var3( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Trinv_lu_opc_var3( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Trinv_lu_opz_var3( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_lu_ops_var3( integer mn_A,
                                 float* buff_A, integer rs_A, integer cs_A )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  integer       i;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    integer       mn_ahead  = mn_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    bl1_sscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a21, rs_A );

    // FLA_Ger_external( FLA_ONE, a21, a10t, A20 );
    bl1_sger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              mn_ahead,
              mn_behind,
              buff_1,
              a21, rs_A,
              a10t, cs_A,
              A20, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_lu_opd_var3( integer mn_A,
                                 double* buff_A, integer rs_A, integer cs_A )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  integer       i;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    integer       mn_ahead  = mn_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    bl1_dscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a21, rs_A );

    // FLA_Ger_external( FLA_ONE, a21, a10t, A20 );
    bl1_dger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              mn_ahead,
              mn_behind,
              buff_1,
              a21, rs_A,
              a10t, cs_A,
              A20, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_lu_opc_var3( integer mn_A,
                                 scomplex* buff_A, integer rs_A, integer cs_A )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  integer       i;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    integer       mn_ahead  = mn_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    bl1_cscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a21, rs_A );

    // FLA_Ger_external( FLA_ONE, a21, a10t, A20 );
    bl1_cger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              mn_ahead,
              mn_behind,
              buff_1,
              a21, rs_A,
              a10t, cs_A,
              A20, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_lu_opz_var3( integer mn_A,
                                 dcomplex* buff_A, integer rs_A, integer cs_A )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  integer       i;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    integer       mn_ahead  = mn_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    bl1_zscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a21, rs_A );

    // FLA_Ger_external( FLA_ONE, a21, a10t, A20 );
    bl1_zger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              mn_ahead,
              mn_behind,
              buff_1,
              a21, rs_A,
              a10t, cs_A,
              A20, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

