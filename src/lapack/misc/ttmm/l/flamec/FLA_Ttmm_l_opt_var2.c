
#include "FLAME.h"

FLA_Error FLA_Ttmm_l_opt_var2( FLA_Obj A )
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

      FLA_Ttmm_l_ops_var2( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Ttmm_l_opd_var2( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Ttmm_l_opc_var2( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Ttmm_l_opz_var2( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_l_ops_var2( int mn_A,
                               float* buff_A, int rs_A, int cs_A )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( alpha11, a10t );
    bl1_sscalv( BLIS1_NO_CONJUGATE,
                mn_behind,
                alpha11,
                a10t, cs_A );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A20, a21, FLA_ONE, a10t );
    bl1_sgemv( BLIS1_TRANSPOSE,
               BLIS1_CONJUGATE,
               mn_ahead,
               mn_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               a10t, cs_A );

    // FLA_Absolute_square( alpha11 );
    bl1_sabsqr( alpha11 );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_ONE, a21, a21, FLA_ONE, alpha11 );
    bl1_sdots( BLIS1_CONJUGATE,
               mn_ahead,
               buff_1,
               a21, rs_A,
               a21, rs_A,
               buff_1,
               alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_l_opd_var2( int mn_A,
                               double* buff_A, int rs_A, int cs_A )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( alpha11, a10t );
    bl1_dscalv( BLIS1_NO_CONJUGATE,
                mn_behind,
                alpha11,
                a10t, cs_A );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A20, a21, FLA_ONE, a10t );
    bl1_dgemv( BLIS1_TRANSPOSE,
               BLIS1_CONJUGATE,
               mn_ahead,
               mn_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               a10t, cs_A );

    // FLA_Absolute_square( alpha11 );
    bl1_dabsqr( alpha11 );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_ONE, a21, a21, FLA_ONE, alpha11 );
    bl1_ddots( BLIS1_CONJUGATE,
               mn_ahead,
               buff_1,
               a21, rs_A,
               a21, rs_A,
               buff_1,
               alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_l_opc_var2( int mn_A,
                               scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( alpha11, a10t );
    bl1_cscalv( BLIS1_NO_CONJUGATE,
                mn_behind,
                alpha11,
                a10t, cs_A );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A20, a21, FLA_ONE, a10t );
    bl1_cgemv( BLIS1_TRANSPOSE,
               BLIS1_CONJUGATE,
               mn_ahead,
               mn_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               a10t, cs_A );

    // FLA_Absolute_square( alpha11 );
    bl1_cabsqr( alpha11 );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_ONE, a21, a21, FLA_ONE, alpha11 );
    bl1_cdots( BLIS1_CONJUGATE,
               mn_ahead,
               buff_1,
               a21, rs_A,
               a21, rs_A,
               buff_1,
               alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Ttmm_l_opz_var2( int mn_A,
                               dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( alpha11, a10t );
    bl1_zscalv( BLIS1_NO_CONJUGATE,
                mn_behind,
                alpha11,
                a10t, cs_A );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A20, a21, FLA_ONE, a10t );
    bl1_zgemv( BLIS1_TRANSPOSE,
               BLIS1_CONJUGATE,
               mn_ahead,
               mn_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               a10t, cs_A );

    // FLA_Absolute_square( alpha11 );
    bl1_zabsqr( alpha11 );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_ONE, a21, a21, FLA_ONE, alpha11 );
    bl1_zdots( BLIS1_CONJUGATE,
               mn_ahead,
               buff_1,
               a21, rs_A,
               a21, rs_A,
               buff_1,
               alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

