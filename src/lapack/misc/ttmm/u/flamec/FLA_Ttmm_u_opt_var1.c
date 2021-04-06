/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Ttmm_u_opt_var1( FLA_Obj A )
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

      FLA_Ttmm_u_ops_var1( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Ttmm_u_opd_var1( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Ttmm_u_opc_var1( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Ttmm_u_opz_var1( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_u_ops_var1( integer mn_A,
                               float* buff_A, integer rs_A, integer cs_A )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  integer       i;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Her_external( FLA_UPPER_TRIANGULAR, FLA_ONE, a01, A00 );
    bl1_ssyr( BLIS1_UPPER_TRIANGULAR,
              mn_behind,
              buff_1,
              a01, rs_A,
              A00, rs_A, cs_A );

    // FLA_Scal_external( alpha11, a01 );
    bl1_sscalv( BLIS1_NO_CONJUGATE,
                mn_behind,
                alpha11,
                a01, rs_A );

    // FLA_Absolute_square( alpha11 );
    bl1_sabsqr( alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_u_opd_var1( integer mn_A,
                               double* buff_A, integer rs_A, integer cs_A )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  integer       i;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Her_external( FLA_UPPER_TRIANGULAR, FLA_ONE, a01, A00 );
    bl1_dsyr( BLIS1_UPPER_TRIANGULAR,
              mn_behind,
              buff_1,
              a01, rs_A,
              A00, rs_A, cs_A );

    // FLA_Scal_external( alpha11, a01 );
    bl1_dscalv( BLIS1_NO_CONJUGATE,
                mn_behind,
                alpha11,
                a01, rs_A );

    // FLA_Absolute_square( alpha11 );
    bl1_dabsqr( alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_u_opc_var1( integer mn_A,
                               scomplex* buff_A, integer rs_A, integer cs_A )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  integer       i;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Her_external( FLA_UPPER_TRIANGULAR, FLA_ONE, a01, A00 );
    bl1_cher( BLIS1_UPPER_TRIANGULAR,
              BLIS1_NO_CONJUGATE,
              mn_behind,
              buff_1,
              a01, rs_A,
              A00, rs_A, cs_A );

    // FLA_Scal_external( alpha11, a01 );
    bl1_cscalv( BLIS1_NO_CONJUGATE,
                mn_behind,
                alpha11,
                a01, rs_A );

    // FLA_Absolute_square( alpha11 );
    bl1_cabsqr( alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_u_opz_var1( integer mn_A,
                               dcomplex* buff_A, integer rs_A, integer cs_A )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  integer       i;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Her_external( FLA_UPPER_TRIANGULAR, FLA_ONE, a01, A00 );
    bl1_zher( BLIS1_UPPER_TRIANGULAR,
              BLIS1_NO_CONJUGATE,
              mn_behind,
              buff_1,
              a01, rs_A,
              A00, rs_A, cs_A );

    // FLA_Scal_external( alpha11, a01 );
    bl1_zscalv( BLIS1_NO_CONJUGATE,
                mn_behind,
                alpha11,
                a01, rs_A );

    // FLA_Absolute_square( alpha11 );
    bl1_zabsqr( alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
