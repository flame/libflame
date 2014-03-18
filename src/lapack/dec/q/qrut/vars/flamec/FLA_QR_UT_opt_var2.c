
#include "FLAME.h"

FLA_Error FLA_QR_UT_opt_var2( FLA_Obj A, FLA_Obj T )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          rs_T, cs_T;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_T = FLA_FLOAT_PTR( T );

      FLA_QR_UT_ops_var2( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_QR_UT_opd_var2( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_QR_UT_opc_var2( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_QR_UT_opz_var2( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR_UT_ops_var2( int m_A,
                              int n_A,
                              float* buff_A, int rs_A, int cs_A, 
                              float* buff_T, int rs_T, int cs_T )
{
  float* buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  int    min_m_n = min( m_A, n_A );
  int    i;

  for ( i = 0; i < min_m_n; ++i )
  {
    float* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    float* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    float* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int    m_ahead  = m_A - i - 1;
    int    n_ahead  = n_A - i - 1;
    int    n_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    FLA_Househ2_UT_l_ops( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21, a12t,
    //                                        A22 );
    FLA_Apply_H2_UT_l_ops_var1( m_ahead,
                                n_ahead,
                                tau11,
                                a21, rs_A,
                                a12t, cs_A,
                                A22, rs_A, cs_A );

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    bl1_scopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01, rs_T );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR_UT_opd_var2( int m_A,
                              int n_A,
                              double* buff_A, int rs_A, int cs_A, 
                              double* buff_T, int rs_T, int cs_T )
{
  double* buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  int     min_m_n = min( m_A, n_A );
  int     i;

  for ( i = 0; i < min_m_n; ++i )
  {
    double* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    double* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    double* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int     m_ahead  = m_A - i - 1;
    int     n_ahead  = n_A - i - 1;
    int     n_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    FLA_Househ2_UT_l_opd( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21, a12t,
    //                                        A22 );
    FLA_Apply_H2_UT_l_opd_var1( m_ahead,
                                n_ahead,
                                tau11,
                                a21, rs_A,
                                a12t, cs_A,
                                A22, rs_A, cs_A );

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    bl1_dcopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01, rs_T );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR_UT_opc_var2( int m_A,
                              int n_A,
                              scomplex* buff_A, int rs_A, int cs_A, 
                              scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    scomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    FLA_Househ2_UT_l_opc( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21, a12t,
    //                                        A22 );
    FLA_Apply_H2_UT_l_opc_var1( m_ahead,
                                n_ahead,
                                tau11,
                                a21, rs_A,
                                a12t, cs_A,
                                A22, rs_A, cs_A );

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    bl1_ccopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01, rs_T );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR_UT_opz_var2( int m_A,
                              int n_A,
                              dcomplex* buff_A, int rs_A, int cs_A, 
                              dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    dcomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    FLA_Househ2_UT_l_opz( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21, a12t,
    //                                        A22 );
    FLA_Apply_H2_UT_l_opz_var1( m_ahead,
                                n_ahead,
                                tau11,
                                a21, rs_A,
                                a12t, cs_A,
                                A22, rs_A, cs_A );

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    bl1_zcopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01, rs_T );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

