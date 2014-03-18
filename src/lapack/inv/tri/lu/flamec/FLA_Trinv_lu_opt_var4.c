
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Trinv_lu_opt_var4( FLA_Obj A )
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

      FLA_Trinv_lu_ops_var4( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Trinv_lu_opd_var4( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Trinv_lu_opc_var4( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Trinv_lu_opz_var4( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_lu_ops_var4( int mn_A,
                                 float* buff_A, int rs_A, int cs_A )
{
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A22, a21 );
    bl1_sscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a21, rs_A );
    bl1_strsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a21, rs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a21, a10t, A20 );
    bl1_sger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              mn_ahead,
              mn_behind,
              buff_m1,
              a21, rs_A,
              a10t, cs_A,
              A20, rs_A, cs_A );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_UNIT_DIAG, A00, a10t );
    bl1_strmv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_lu_opd_var4( int mn_A,
                                 double* buff_A, int rs_A, int cs_A )
{
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A22, a21 );
    bl1_dscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a21, rs_A );
    bl1_dtrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a21, rs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a21, a10t, A20 );
    bl1_dger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              mn_ahead,
              mn_behind,
              buff_m1,
              a21, rs_A,
              a10t, cs_A,
              A20, rs_A, cs_A );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_UNIT_DIAG, A00, a10t );
    bl1_dtrmv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_lu_opc_var4( int mn_A,
                                 scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A22, a21 );
    bl1_cscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a21, rs_A );
    bl1_ctrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a21, rs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a21, a10t, A20 );
    bl1_cger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              mn_ahead,
              mn_behind,
              buff_m1,
              a21, rs_A,
              a10t, cs_A,
              A20, rs_A, cs_A );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_UNIT_DIAG, A00, a10t );
    bl1_ctrmv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_lu_opz_var4( int mn_A,
                                 dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a21 );
    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A22, a21 );
    bl1_zscalv( BLIS1_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a21, rs_A );
    bl1_ztrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a21, rs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a21, a10t, A20 );
    bl1_zger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              mn_ahead,
              mn_behind,
              buff_m1,
              a21, rs_A,
              a10t, cs_A,
              A20, rs_A, cs_A );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_UNIT_DIAG, A00, a10t );
    bl1_ztrmv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
