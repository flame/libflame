
#include "FLAME.h"

FLA_Error FLA_Trinv_ln_opt_var3( FLA_Obj A )
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

      FLA_Trinv_ln_ops_var3( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Trinv_ln_opd_var3( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Trinv_ln_opc_var3( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Trinv_ln_opz_var3( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_ln_ops_var3( int mn_A,
                                 float* buff_A, int rs_A, int cs_A )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float     alpha11_m1;
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

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_sneg2( alpha11, &alpha11_m1 );
    bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
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

    // FLA_Inv_scal_external( alpha11, a10t );
    bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                   mn_behind,
                   alpha11,
                   a10t, cs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bl1_sinverts( BLIS1_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_ln_opd_var3( int mn_A,
                                 double* buff_A, int rs_A, int cs_A )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double    alpha11_m1;
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

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_dneg2( alpha11, &alpha11_m1 );
    bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
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

    // FLA_Inv_scal_external( alpha11, a10t );
    bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                   mn_behind,
                   alpha11,
                   a10t, cs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bl1_dinverts( BLIS1_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_ln_opc_var3( int mn_A,
                                 scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex  alpha11_m1;
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

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_cneg2( alpha11, &alpha11_m1 );
    bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
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

    // FLA_Inv_scal_external( alpha11, a10t );
    bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                   mn_behind,
                   alpha11,
                   a10t, cs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bl1_cinverts( BLIS1_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_ln_opz_var3( int mn_A,
                                 dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex  alpha11_m1;
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

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_zneg2( alpha11, &alpha11_m1 );
    bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
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

    // FLA_Inv_scal_external( alpha11, a10t );
    bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                   mn_behind,
                   alpha11,
                   a10t, cs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bl1_zinverts( BLIS1_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

