/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Ttmm_u_opt_var3( FLA_Obj A )
{
  FLA_Datatype datatype;
  int          mn_A;
  int          rs_A, cs_A;

  datatype = FLA_Obj_datatype( A );

  mn_A     = FLA_Obj_length( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );

      FLA_Ttmm_u_ops_var3( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Ttmm_u_opd_var3( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Ttmm_u_opc_var3( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Ttmm_u_opz_var3( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_u_ops_var3( int mn_A,
                               float* buff_A, int rs_A, int cs_A )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Absolute_square( alpha11 );
    bl1_sabsqr( alpha11 );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_ONE, a12t, a12t, FLA_ONE, alpha11 );
    bl1_sdots( BLIS1_CONJUGATE,
               mn_ahead,
               buff_1,
               a12t, cs_A,
               a12t, cs_A,
               buff_1,
               alpha11 );

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bl1_strmv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_u_opd_var3( int mn_A,
                               double* buff_A, int rs_A, int cs_A )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Absolute_square( alpha11 );
    bl1_dabsqr( alpha11 );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_ONE, a12t, a12t, FLA_ONE, alpha11 );
    bl1_ddots( BLIS1_CONJUGATE,
               mn_ahead,
               buff_1,
               a12t, cs_A,
               a12t, cs_A,
               buff_1,
               alpha11 );

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bl1_dtrmv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_u_opc_var3( int mn_A,
                               scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Absolute_square( alpha11 );
    bl1_cabsqr( alpha11 );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_ONE, a12t, a12t, FLA_ONE, alpha11 );
    bl1_cdots( BLIS1_CONJUGATE,
               mn_ahead,
               buff_1,
               a12t, cs_A,
               a12t, cs_A,
               buff_1,
               alpha11 );

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bl1_ctrmv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_u_opz_var3( int mn_A,
                               dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Absolute_square( alpha11 );
    bl1_zabsqr( alpha11 );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_ONE, a12t, a12t, FLA_ONE, alpha11 );
    bl1_zdots( BLIS1_CONJUGATE,
               mn_ahead,
               buff_1,
               a12t, cs_A,
               a12t, cs_A,
               buff_1,
               alpha11 );

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bl1_ztrmv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
