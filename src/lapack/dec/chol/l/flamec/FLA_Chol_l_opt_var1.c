
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Chol_l_opt_var1( FLA_Obj A )
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

      r_val = FLA_Chol_l_ops_var1( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      r_val = FLA_Chol_l_opd_var1( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      r_val = FLA_Chol_l_opc_var1( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      r_val = FLA_Chol_l_opz_var1( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }
  }

  return r_val;
}



FLA_Error FLA_Chol_l_ops_var1( int mn_A,
                               float* buff_A, int rs_A, int cs_A )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a10t );
    bl1_strsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a10t, a10t, FLA_ONE, alpha11 );
    bl1_sdots( BLIS1_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a10t, cs_A,
               buff_1,
               alpha11 );

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_ssqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_l_opd_var1( int mn_A,
                               double* buff_A, int rs_A, int cs_A )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a10t );
    bl1_dtrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a10t, a10t, FLA_ONE, alpha11 );
    bl1_ddots( BLIS1_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a10t, cs_A,
               buff_1,
               alpha11 );

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_dsqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_l_opc_var1( int mn_A,
                               scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a10t );
    bl1_ctrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a10t, a10t, FLA_ONE, alpha11 );
    bl1_cdots( BLIS1_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a10t, cs_A,
               buff_1,
               alpha11 );

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_csqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_l_opz_var1( int mn_A,
                               dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a10t );
    bl1_ztrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a10t, a10t, FLA_ONE, alpha11 );
    bl1_zdots( BLIS1_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a10t, cs_A,
               buff_1,
               alpha11 );

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_zsqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
