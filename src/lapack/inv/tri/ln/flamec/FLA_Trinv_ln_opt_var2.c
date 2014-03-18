
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Trinv_ln_opt_var2( FLA_Obj A )
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

      FLA_Trinv_ln_ops_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Trinv_ln_opd_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Trinv_ln_opc_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Trinv_ln_opz_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_ln_ops_var2( int mn_A,
                                 float* buff_A, int rs_A, int cs_A )
{
  float     alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a21 );
    bl1_strsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a21, rs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_sneg2( alpha11, &alpha11_m1 );
    bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
                   a21, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bl1_sinverts( BLIS1_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_ln_opd_var2( int mn_A,
                                 double* buff_A, int rs_A, int cs_A )
{
  double    alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a21 );
    bl1_dtrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a21, rs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_dneg2( alpha11, &alpha11_m1 );
    bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
                   a21, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bl1_dinverts( BLIS1_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_ln_opc_var2( int mn_A,
                                 scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex  alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a21 );
    bl1_ctrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a21, rs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_cneg2( alpha11, &alpha11_m1 );
    bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
                   a21, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bl1_cinverts( BLIS1_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_ln_opz_var2( int mn_A,
                                 dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex  alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a21 );
    bl1_ztrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a21, rs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_zneg2( alpha11, &alpha11_m1 );
    bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
                   a21, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bl1_zinverts( BLIS1_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
