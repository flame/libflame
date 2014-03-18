
#include "FLAME.h"

FLA_Error FLA_LQ_UT_opt_var2( FLA_Obj A, FLA_Obj T )
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

      FLA_LQ_UT_ops_var2( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_LQ_UT_opd_var2( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_LQ_UT_opc_var2( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_LQ_UT_opz_var2( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LQ_UT_ops_var2( int m_A,
                              int n_A,
                              float* buff_A, int rs_A, int cs_A,
                              float* buff_T, int rs_T, int cs_T )
{
  float* buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  int    min_m_n = min( m_A, n_A );
  int    i;

  for ( i = 0; i < min_m_n; ++i )
  {
    float* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    float* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    float* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int    m_ahead  = m_A - i - 1;
    int    n_ahead  = n_A - i - 1;
    int    m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_RIGHT, alpha11, a12t
    //                 tau11 );
    FLA_Househ2_UT_r_ops( n_ahead,
                          alpha11,
                          a12t, cs_A,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_RIGHT, tau11, a12t, a21, A22 );
    FLA_Apply_H2_UT_r_ops_var1( m_ahead,
                                n_ahead,
                                tau11,
                                a12t, cs_A,
                                a21, rs_A,
                                A22, rs_A, cs_A );

    // FLA_Copyt_external( FLA_CONJ_NO_TRANSPOSE, a01, t01 );
    bl1_scopyv( BLIS1_CONJUGATE,
                m_behind,
                a01, rs_A,
                t01, rs_T );

    // FLA_Gemvc_external( FLA_CONJ_NO_TRANSPOSE, FLA_NO_CONJUGATE,
    //                     FLA_ONE, A02, a12t, FLA_ONE, t01 );
    bl1_sgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_1,
               A02, rs_A, cs_A,
               a12t, cs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LQ_UT_opd_var2( int m_A,
                              int n_A,
                              double* buff_A, int rs_A, int cs_A,
                              double* buff_T, int rs_T, int cs_T )
{
  double* buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  int     min_m_n = min( m_A, n_A );
  int     i;

  for ( i = 0; i < min_m_n; ++i )
  {
    double* a01     = buff_A + (i  )*cs_A + (0  )*rs_A;
    double* alpha11 = buff_A + (i  )*cs_A + (i  )*rs_A;
    double* a21     = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double* A02     = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double* a12t    = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double* A22     = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double* tau11   = buff_T + (i  )*cs_T + (i  )*rs_T;
    double* t01     = buff_T + (i  )*cs_T + (0  )*rs_T;

    int     m_ahead  = m_A - i - 1;
    int     n_ahead  = n_A - i - 1;
    int     m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_RIGHT, alpha11, a12t
    //                 tau11 );
    FLA_Househ2_UT_r_opd( n_ahead,
                        alpha11,
                        a12t, cs_A,
                        tau11 );

    // FLA_Apply_H2_UT( FLA_RIGHT, tau11, a12t, a21, A22 );
    FLA_Apply_H2_UT_r_opd_var1( m_ahead,
                                n_ahead,
                                tau11,
                                a12t, cs_A,
                                a21, rs_A,
                                A22, rs_A, cs_A );

    // FLA_Copyt_external( FLA_CONJ_NO_TRANSPOSE, a01, t01 );
    bl1_dcopyv( BLIS1_CONJUGATE,
                m_behind,
                a01, rs_A,
                t01, rs_T );

    // FLA_Gemvc_external( FLA_CONJ_NO_TRANSPOSE, FLA_NO_CONJUGATE,
    //                     FLA_ONE, A02, a12t, FLA_ONE, t01 );
    bl1_dgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_1,
               A02, rs_A, cs_A,
               a12t, cs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LQ_UT_opc_var2( int m_A,
                              int n_A,
                              scomplex* buff_A, int rs_A, int cs_A,
                              scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    scomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_RIGHT, alpha11, a12t
    //                 tau11 );
    FLA_Househ2_UT_r_opc( n_ahead,
                          alpha11,
                          a12t, cs_A,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_RIGHT, tau11, a12t, a21, A22 );
    FLA_Apply_H2_UT_r_opc_var1( m_ahead,
                                n_ahead,
                                tau11,
                                a12t, cs_A,
                                a21, rs_A,
                                A22, rs_A, cs_A );

    // FLA_Copyt_external( FLA_CONJ_NO_TRANSPOSE, a01, t01 );
    bl1_ccopyv( BLIS1_CONJUGATE,
                m_behind,
                a01, rs_A,
                t01, rs_T );

    // FLA_Gemvc_external( FLA_CONJ_NO_TRANSPOSE, FLA_NO_CONJUGATE,
    //                     FLA_ONE, A02, a12t, FLA_ONE, t01 );
    bl1_cgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_1,
               A02, rs_A, cs_A,
               a12t, cs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LQ_UT_opz_var2( int m_A,
                              int n_A,
                              dcomplex* buff_A, int rs_A, int cs_A,
                              dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    dcomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_RIGHT, alpha11, a12t
    //                 tau11 );
    FLA_Househ2_UT_r_opz( n_ahead,
                          alpha11,
                          a12t, cs_A,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_RIGHT, tau11, a12t, a21, A22 );
    FLA_Apply_H2_UT_r_opz_var1( m_ahead,
                                n_ahead,
                                tau11,
                                a12t, cs_A,
                                a21, rs_A,
                                A22, rs_A, cs_A );

    // FLA_Copyt_external( FLA_CONJ_NO_TRANSPOSE, a01, t01 );
    bl1_zcopyv( BLIS1_CONJUGATE,
                m_behind,
                a01, rs_A,
                t01, rs_T );

    // FLA_Gemvc_external( FLA_CONJ_NO_TRANSPOSE, FLA_NO_CONJUGATE,
    //                     FLA_ONE, A02, a12t, FLA_ONE, t01 );
    bl1_zgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_1,
               A02, rs_A, cs_A,
               a12t, cs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

