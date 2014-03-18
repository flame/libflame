
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Trinv_uu_opt_var2( FLA_Obj A )
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

      FLA_Trinv_uu_ops_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Trinv_uu_opd_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Trinv_uu_opc_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Trinv_uu_opz_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_uu_ops_var2( int mn_A,
                                 float* buff_A, int rs_A, int cs_A )
{
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_UNIT_DIAG, A22, a12t );
    bl1_strsv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    bl1_sscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_uu_opd_var2( int mn_A,
                                 double* buff_A, int rs_A, int cs_A )
{
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_UNIT_DIAG, A22, a12t );
    bl1_dtrsv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    bl1_dscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_uu_opc_var2( int mn_A,
                                 scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_UNIT_DIAG, A22, a12t );
    bl1_ctrsv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    bl1_cscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_uu_opz_var2( int mn_A,
                                 dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_UNIT_DIAG, A22, a12t );
    bl1_ztrsv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    bl1_zscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
