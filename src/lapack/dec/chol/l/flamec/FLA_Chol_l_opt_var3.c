
#include "FLAME.h"

FLA_Error FLA_Chol_l_opt_var3( FLA_Obj A )
{
  FLA_Error    r_val = FLA_SUCCESS;
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

      r_val = FLA_Chol_l_ops_var3( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      r_val = FLA_Chol_l_opd_var3( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      r_val = FLA_Chol_l_opc_var3( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      r_val = FLA_Chol_l_opz_var3( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }
  }

  return r_val;
}



FLA_Error FLA_Chol_l_ops_var3( int mn_A,
                               float* buff_A, int rs_A, int cs_A )
{
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_ssqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a21, rs_A );

    // FLA_Her_external( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, a21, A22 );
    bl1_ssyr( BLIS1_LOWER_TRIANGULAR,
              mn_ahead,
              buff_m1,
              a21, rs_A,
              A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_l_opd_var3( int mn_A,
                               double* buff_A, int rs_A, int cs_A )
{
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_dsqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a21, rs_A );

    // FLA_Her_external( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, a21, A22 );
    bl1_dsyr( BLIS1_LOWER_TRIANGULAR,
              mn_ahead,
              buff_m1,
              a21, rs_A,
              A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_l_opc_var3( int mn_A,
                               scomplex* buff_A, int rs_A, int cs_A )
{
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_csqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a21, rs_A );

    // FLA_Her_external( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, a21, A22 );
    bl1_cher( BLIS1_LOWER_TRIANGULAR,
              BLIS1_NO_CONJUGATE,
              mn_ahead,
              buff_m1,
              a21, rs_A,
              A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_l_opz_var3( int mn_A,
                               dcomplex* buff_A, int rs_A, int cs_A )
{
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_zsqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a21, rs_A );

    // FLA_Her_external( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, a21, A22 );
    bl1_zher( BLIS1_LOWER_TRIANGULAR,
              BLIS1_NO_CONJUGATE,
              mn_ahead,
              buff_m1,
              a21, rs_A,
              A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

