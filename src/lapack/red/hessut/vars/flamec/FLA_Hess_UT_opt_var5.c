
#include "FLAME.h"

FLA_Error FLA_Hess_UT_opt_var5( FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;
  FLA_Obj   U, Z;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &U );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Z );

  r_val = FLA_Hess_UT_step_opt_var5( A, U, Z, T );

  FLA_Obj_free( &U );
  FLA_Obj_free( &Z );

  return r_val;
}

FLA_Error FLA_Hess_UT_step_opt_var5( FLA_Obj A, FLA_Obj U, FLA_Obj Z, FLA_Obj T )
{
  FLA_Datatype datatype;
  int          m_A, m_T;
  int          rs_A, cs_A;
  int          rs_U, cs_U;
  int          rs_Z, cs_Z;
  int          rs_T, cs_T;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  m_T      = FLA_Obj_length( T );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_U     = FLA_Obj_row_stride( U );
  cs_U     = FLA_Obj_col_stride( U );

  rs_Z     = FLA_Obj_row_stride( Z );
  cs_Z     = FLA_Obj_col_stride( Z );

  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_U = FLA_FLOAT_PTR( U );
      float* buff_Z = FLA_FLOAT_PTR( Z );
      float* buff_T = FLA_FLOAT_PTR( T );

      FLA_Hess_UT_step_ops_var5( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_U, rs_U, cs_U,
                                 buff_Z, rs_Z, cs_Z,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_U = FLA_DOUBLE_PTR( U );
      double* buff_Z = FLA_DOUBLE_PTR( Z );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_Hess_UT_step_opd_var5( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_U, rs_U, cs_U,
                                 buff_Z, rs_Z, cs_Z,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_U = FLA_COMPLEX_PTR( U );
      scomplex* buff_Z = FLA_COMPLEX_PTR( Z );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_Hess_UT_step_opc_var5( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_U, rs_U, cs_U,
                                 buff_Z, rs_Z, cs_Z,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_U = FLA_DOUBLE_COMPLEX_PTR( U );
      dcomplex* buff_Z = FLA_DOUBLE_COMPLEX_PTR( Z );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_Hess_UT_step_opz_var5( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_U, rs_U, cs_U,
                                 buff_Z, rs_Z, cs_Z,
                                 buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_ops_var5( int m_A,
                                     int m_T,
                                     float* buff_A, int rs_A, int cs_A, 
                                     float* buff_U, int rs_U, int cs_U,
                                     float* buff_Z, int rs_Z, int cs_Z,
                                     float* buff_T, int rs_T, int cs_T )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  float*    buff_w = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_w  = 1;

  // FLA_Set( FLA_ZERO, U );
  // FLA_Set( FLA_ZERO, Z );
  bl1_ssetm( m_A,
             b_alg,
             buff_0,
             buff_U, rs_U, cs_U );
  bl1_ssetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    float*    a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    U00      = buff_U + (0  )*cs_U + (0  )*rs_U;
    float*    u10t     = buff_U + (0  )*cs_U + (i  )*rs_U;
    float*    U20      = buff_U + (0  )*cs_U + (i+1)*rs_U;
    float*    u21      = buff_U + (i  )*cs_U + (i+1)*rs_U;

    float*    Z00      = buff_Z + (0  )*cs_Z + (0  )*rs_Z;
    float*    z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    float*    Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    float*    z01      = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    float*    zeta11   = buff_Z + (i  )*cs_Z + (i  )*rs_Z;
    float*    z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    float*    T00      = buff_T + (0  )*cs_T + (0  )*rs_T;
    float*    t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    float*    tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    float*    w0       = buff_w + (0  )*inc_w;

    float*    a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    float*    a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    float*    u21_t    = u21    + (0  )*cs_U + (0  )*rs_U;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, u10t, w0 );
      // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
      //           T00, w0 );
      bl1_scopyv( BLIS1_CONJUGATE,
                  m_behind,
                  u10t, cs_U,
                  w0,   inc_w );
      bl1_strsv( BLIS1_UPPER_TRIANGULAR,
                 BLIS1_NO_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 T00, rs_T, cs_T,
                 w0,  inc_w );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z00, w0, FLA_ONE, a01 );
      // FLA_Dots( FLA_MINUS_ONE, z10t, w0, FLA_ONE, alpha11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, w0, FLA_ONE, a21 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_behind,
                 buff_m1,
                 Z00,  rs_Z, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 a01,  rs_A );
      bl1_sdots( BLIS1_NO_CONJUGATE,
                 m_behind,
                 buff_m1,
                 z10t, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 alpha11 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20,  rs_Z, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 a21,  rs_A );

      // FLA_Trmvsx( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
      //             FLA_ONE, U00, a01, FLA_ZERO, w0 );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, alpha11, u10t, w0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, U20, a21, FLA_ONE, w0 );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  m_behind,
                  a01, rs_A,
                  w0,  inc_w );
      bl1_strmv( BLIS1_LOWER_TRIANGULAR,
                 BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 U00, rs_U, cs_U,
                 w0,  inc_w );
      bl1_saxpyv( BLIS1_CONJUGATE,
                  m_behind,
                  alpha11,
                  u10t, cs_U,
                  w0,   inc_w );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 U20, rs_U, cs_U,
                 a21, rs_A,
                 buff_1,
                 w0,  inc_w );

      // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
      //           T00, w0 );
      bl1_strsv( BLIS1_UPPER_TRIANGULAR,
                 BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 T00, rs_T, cs_T,
                 w0,  inc_w );

      // FLA_Trmvsx( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
      //             FLA_MINUS_ONE, U00, w0, FLA_ONE, a01 );
      // FLA_Dots( FLA_MINUS_ONE, u10t, w0, FLA_ONE, alpha11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, U20, w0, FLA_ONE, a21 );
      bl1_strmvsx( BLIS1_LOWER_TRIANGULAR,
                   BLIS1_NO_TRANSPOSE,
                   BLIS1_NONUNIT_DIAG,
                   m_behind,
                   buff_m1,
                   U00, rs_U, cs_U,
                   w0,  inc_w,
                   buff_1,
                   a01, rs_A );
      bl1_sdots( BLIS1_NO_CONJUGATE,
                 m_behind,
                 buff_m1,
                 u10t, cs_U,
                 w0,   inc_w,
                 buff_1,
                 alpha11 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 U20, rs_U, cs_U,
                 w0,  inc_w,
                 buff_1,
                 a21, rs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Househ2_UT( FLA_LEFT,
      //                 a21_t,
      //                 a21_b, tau11 );
      FLA_Househ2_UT_l_ops( m_ahead - 1,
                            a21_t,
                            a21_b, rs_A,
                            tau11 );

      // FLA_Copy( a21, u21 );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u21, rs_U );

      // FLA_Set( FLA_ONE, u21_t );
      *u21_t = *buff_1;

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, u21, FLA_ZERO, z01 );
      // FLA_Dot( a12t, u21, zeta11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, u21, FLA_ZERO, z21 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 u21, rs_U,
                 buff_0,
                 z01, rs_Z );
      bl1_sdot( BLIS1_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                u21,  rs_U,
                zeta11 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 u21, rs_U,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, U20, u21, FLA_ZERO, t01 );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 U20, rs_U, cs_U,
                 u21, rs_U,
                 buff_0,
                 t01, rs_T );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opd_var5( int m_A,
                                     int m_T,
                                     double* buff_A, int rs_A, int cs_A, 
                                     double* buff_U, int rs_U, int cs_U,
                                     double* buff_Z, int rs_Z, int cs_Z,
                                     double* buff_T, int rs_T, int cs_T )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  double*   buff_w = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_w  = 1;

  // FLA_Set( FLA_ZERO, U );
  // FLA_Set( FLA_ZERO, Z );
  bl1_dsetm( m_A,
             b_alg,
             buff_0,
             buff_U, rs_U, cs_U );
  bl1_dsetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    double*   a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   U00      = buff_U + (0  )*cs_U + (0  )*rs_U;
    double*   u10t     = buff_U + (0  )*cs_U + (i  )*rs_U;
    double*   U20      = buff_U + (0  )*cs_U + (i+1)*rs_U;
    double*   u21      = buff_U + (i  )*cs_U + (i+1)*rs_U;

    double*   Z00      = buff_Z + (0  )*cs_Z + (0  )*rs_Z;
    double*   z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    double*   Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    double*   z01      = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    double*   zeta11   = buff_Z + (i  )*cs_Z + (i  )*rs_Z;
    double*   z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    double*   T00      = buff_T + (0  )*cs_T + (0  )*rs_T;
    double*   t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    double*   tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    double*   w0       = buff_w + (0  )*inc_w;

    double*   a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    double*   a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    double*   u21_t    = u21    + (0  )*cs_U + (0  )*rs_U;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, u10t, w0 );
      // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
      //           T00, w0 );
      bl1_dcopyv( BLIS1_CONJUGATE,
                  m_behind,
                  u10t, cs_U,
                  w0,   inc_w );
      bl1_dtrsv( BLIS1_UPPER_TRIANGULAR,
                 BLIS1_NO_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 T00, rs_T, cs_T,
                 w0,  inc_w );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z00, w0, FLA_ONE, a01 );
      // FLA_Dots( FLA_MINUS_ONE, z10t, w0, FLA_ONE, alpha11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, w0, FLA_ONE, a21 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_behind,
                 buff_m1,
                 Z00,  rs_Z, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 a01,  rs_A );
      bl1_ddots( BLIS1_NO_CONJUGATE,
                 m_behind,
                 buff_m1,
                 z10t, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 alpha11 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20,  rs_Z, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 a21,  rs_A );

      // FLA_Trmvsx( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
      //             FLA_ONE, U00, a01, FLA_ZERO, w0 );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, alpha11, u10t, w0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, U20, a21, FLA_ONE, w0 );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  m_behind,
                  a01, rs_A,
                  w0,  inc_w );
      bl1_dtrmv( BLIS1_LOWER_TRIANGULAR,
                 BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 U00, rs_U, cs_U,
                 w0,  inc_w );
      bl1_daxpyv( BLIS1_CONJUGATE,
                  m_behind,
                  alpha11,
                  u10t, cs_U,
                  w0,   inc_w );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 U20, rs_U, cs_U,
                 a21, rs_A,
                 buff_1,
                 w0,  inc_w );

      // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
      //           T00, w0 );
      bl1_dtrsv( BLIS1_UPPER_TRIANGULAR,
                 BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 T00, rs_T, cs_T,
                 w0,  inc_w );

      // FLA_Trmvsx( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
      //             FLA_MINUS_ONE, U00, w0, FLA_ONE, a01 );
      // FLA_Dots( FLA_MINUS_ONE, u10t, w0, FLA_ONE, alpha11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, U20, w0, FLA_ONE, a21 );
      bl1_dtrmvsx( BLIS1_LOWER_TRIANGULAR,
                   BLIS1_NO_TRANSPOSE,
                   BLIS1_NONUNIT_DIAG,
                   m_behind,
                   buff_m1,
                   U00, rs_U, cs_U,
                   w0,  inc_w,
                   buff_1,
                   a01, rs_A );
      bl1_ddots( BLIS1_NO_CONJUGATE,
                 m_behind,
                 buff_m1,
                 u10t, cs_U,
                 w0,   inc_w,
                 buff_1,
                 alpha11 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 U20, rs_U, cs_U,
                 w0,  inc_w,
                 buff_1,
                 a21, rs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Househ2_UT( FLA_LEFT,
      //                 a21_t,
      //                 a21_b, tau11 );
      FLA_Househ2_UT_l_opd( m_ahead - 1,
                            a21_t,
                            a21_b, rs_A,
                            tau11 );

      // FLA_Copy( a21, u21 );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u21, rs_U );

      // FLA_Set( FLA_ONE, u21_t );
      *u21_t = *buff_1;

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, u21, FLA_ZERO, z01 );
      // FLA_Dot( a12t, u21, zeta11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, u21, FLA_ZERO, z21 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 u21, rs_U,
                 buff_0,
                 z01, rs_Z );
      bl1_ddot( BLIS1_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                u21,  rs_U,
                zeta11 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 u21, rs_U,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, U20, u21, FLA_ZERO, t01 );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 U20, rs_U, cs_U,
                 u21, rs_U,
                 buff_0,
                 t01, rs_T );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opc_var5( int m_A,
                                     int m_T,
                                     scomplex* buff_A, int rs_A, int cs_A, 
                                     scomplex* buff_U, int rs_U, int cs_U,
                                     scomplex* buff_Z, int rs_Z, int cs_Z,
                                     scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  scomplex* buff_w = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_w  = 1;

  // FLA_Set( FLA_ZERO, U );
  // FLA_Set( FLA_ZERO, Z );
  bl1_csetm( m_A,
             b_alg,
             buff_0,
             buff_U, rs_U, cs_U );
  bl1_csetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* U00      = buff_U + (0  )*cs_U + (0  )*rs_U;
    scomplex* u10t     = buff_U + (0  )*cs_U + (i  )*rs_U;
    scomplex* U20      = buff_U + (0  )*cs_U + (i+1)*rs_U;
    scomplex* u21      = buff_U + (i  )*cs_U + (i+1)*rs_U;

    scomplex* Z00      = buff_Z + (0  )*cs_Z + (0  )*rs_Z;
    scomplex* z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    scomplex* Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    scomplex* z01      = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    scomplex* zeta11   = buff_Z + (i  )*cs_Z + (i  )*rs_Z;
    scomplex* z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    scomplex* T00      = buff_T + (0  )*cs_T + (0  )*rs_T;
    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    scomplex* w0       = buff_w + (0  )*inc_w;

    scomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    scomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    scomplex* u21_t    = u21    + (0  )*cs_U + (0  )*rs_U;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, u10t, w0 );
      // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
      //           T00, w0 );
      bl1_ccopyv( BLIS1_CONJUGATE,
                  m_behind,
                  u10t, cs_U,
                  w0,   inc_w );
      bl1_ctrsv( BLIS1_UPPER_TRIANGULAR,
                 BLIS1_NO_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 T00, rs_T, cs_T,
                 w0,  inc_w );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z00, w0, FLA_ONE, a01 );
      // FLA_Dots( FLA_MINUS_ONE, z10t, w0, FLA_ONE, alpha11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, w0, FLA_ONE, a21 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_behind,
                 buff_m1,
                 Z00,  rs_Z, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 a01,  rs_A );
      bl1_cdots( BLIS1_NO_CONJUGATE,
                 m_behind,
                 buff_m1,
                 z10t, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 alpha11 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20,  rs_Z, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 a21,  rs_A );

      // FLA_Trmvsx( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
      //             FLA_ONE, U00, a01, FLA_ZERO, w0 );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, alpha11, u10t, w0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, U20, a21, FLA_ONE, w0 );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  m_behind,
                  a01, rs_A,
                  w0,  inc_w );
      bl1_ctrmv( BLIS1_LOWER_TRIANGULAR,
                 BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 U00, rs_U, cs_U,
                 w0,  inc_w );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  m_behind,
                  alpha11,
                  u10t, cs_U,
                  w0,   inc_w );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 U20, rs_U, cs_U,
                 a21, rs_A,
                 buff_1,
                 w0,  inc_w );

      // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
      //           T00, w0 );
      bl1_ctrsv( BLIS1_UPPER_TRIANGULAR,
                 BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 T00, rs_T, cs_T,
                 w0,  inc_w );

      // FLA_Trmvsx( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
      //             FLA_MINUS_ONE, U00, w0, FLA_ONE, a01 );
      // FLA_Dots( FLA_MINUS_ONE, u10t, w0, FLA_ONE, alpha11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, U20, w0, FLA_ONE, a21 );
      bl1_ctrmvsx( BLIS1_LOWER_TRIANGULAR,
                   BLIS1_NO_TRANSPOSE,
                   BLIS1_NONUNIT_DIAG,
                   m_behind,
                   buff_m1,
                   U00, rs_U, cs_U,
                   w0,  inc_w,
                   buff_1,
                   a01, rs_A );
      bl1_cdots( BLIS1_NO_CONJUGATE,
                 m_behind,
                 buff_m1,
                 u10t, cs_U,
                 w0,   inc_w,
                 buff_1,
                 alpha11 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 U20, rs_U, cs_U,
                 w0,  inc_w,
                 buff_1,
                 a21, rs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Househ2_UT( FLA_LEFT,
      //                 a21_t,
      //                 a21_b, tau11 );
      FLA_Househ2_UT_l_opc( m_ahead - 1,
                            a21_t,
                            a21_b, rs_A,
                            tau11 );

      // FLA_Copy( a21, u21 );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u21, rs_U );

      // FLA_Set( FLA_ONE, u21_t );
      *u21_t = *buff_1;

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, u21, FLA_ZERO, z01 );
      // FLA_Dot( a12t, u21, zeta11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, u21, FLA_ZERO, z21 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 u21, rs_U,
                 buff_0,
                 z01, rs_Z );
      bl1_cdot( BLIS1_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                u21,  rs_U,
                zeta11 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 u21, rs_U,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, U20, u21, FLA_ZERO, t01 );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 U20, rs_U, cs_U,
                 u21, rs_U,
                 buff_0,
                 t01, rs_T );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opz_var5( int m_A,
                                     int m_T,
                                     dcomplex* buff_A, int rs_A, int cs_A, 
                                     dcomplex* buff_U, int rs_U, int cs_U,
                                     dcomplex* buff_Z, int rs_Z, int cs_Z,
                                     dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  dcomplex* buff_w = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_w  = 1;

  // FLA_Set( FLA_ZERO, U );
  // FLA_Set( FLA_ZERO, Z );
  bl1_zsetm( m_A,
             b_alg,
             buff_0,
             buff_U, rs_U, cs_U );
  bl1_zsetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* U00      = buff_U + (0  )*cs_U + (0  )*rs_U;
    dcomplex* u10t     = buff_U + (0  )*cs_U + (i  )*rs_U;
    dcomplex* U20      = buff_U + (0  )*cs_U + (i+1)*rs_U;
    dcomplex* u21      = buff_U + (i  )*cs_U + (i+1)*rs_U;

    dcomplex* Z00      = buff_Z + (0  )*cs_Z + (0  )*rs_Z;
    dcomplex* z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    dcomplex* Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    dcomplex* z01      = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    dcomplex* zeta11   = buff_Z + (i  )*cs_Z + (i  )*rs_Z;
    dcomplex* z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    dcomplex* T00      = buff_T + (0  )*cs_T + (0  )*rs_T;
    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    dcomplex* w0       = buff_w + (0  )*inc_w;

    dcomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    dcomplex* u21_t    = u21    + (0  )*cs_U + (0  )*rs_U;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, u10t, w0 );
      // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
      //           T00, w0 );
      bl1_zcopyv( BLIS1_CONJUGATE,
                  m_behind,
                  u10t, cs_U,
                  w0,   inc_w );
      bl1_ztrsv( BLIS1_UPPER_TRIANGULAR,
                 BLIS1_NO_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 T00, rs_T, cs_T,
                 w0,  inc_w );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z00, w0, FLA_ONE, a01 );
      // FLA_Dots( FLA_MINUS_ONE, z10t, w0, FLA_ONE, alpha11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, w0, FLA_ONE, a21 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_behind,
                 buff_m1,
                 Z00,  rs_Z, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 a01,  rs_A );
      bl1_zdots( BLIS1_NO_CONJUGATE,
                 m_behind,
                 buff_m1,
                 z10t, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 alpha11 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20,  rs_Z, cs_Z,
                 w0,   inc_w,
                 buff_1,
                 a21,  rs_A );

      // FLA_Trmvsx( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
      //             FLA_ONE, U00, a01, FLA_ZERO, w0 );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, alpha11, u10t, w0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, U20, a21, FLA_ONE, w0 );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  m_behind,
                  a01, rs_A,
                  w0,  inc_w );
      bl1_ztrmv( BLIS1_LOWER_TRIANGULAR,
                 BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 U00, rs_U, cs_U,
                 w0,  inc_w );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  m_behind,
                  alpha11,
                  u10t, cs_U,
                  w0,   inc_w );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 U20, rs_U, cs_U,
                 a21, rs_A,
                 buff_1,
                 w0,  inc_w );

      // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
      //           T00, w0 );
      bl1_ztrsv( BLIS1_UPPER_TRIANGULAR,
                 BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NONUNIT_DIAG,
                 m_behind,
                 T00, rs_T, cs_T,
                 w0,  inc_w );

      // FLA_Trmvsx( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
      //             FLA_MINUS_ONE, U00, w0, FLA_ONE, a01 );
      // FLA_Dots( FLA_MINUS_ONE, u10t, w0, FLA_ONE, alpha11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, U20, w0, FLA_ONE, a21 );
      bl1_ztrmvsx( BLIS1_LOWER_TRIANGULAR,
                   BLIS1_NO_TRANSPOSE,
                   BLIS1_NONUNIT_DIAG,
                   m_behind,
                   buff_m1,
                   U00, rs_U, cs_U,
                   w0,  inc_w,
                   buff_1,
                   a01, rs_A );
      bl1_zdots( BLIS1_NO_CONJUGATE,
                 m_behind,
                 buff_m1,
                 u10t, cs_U,
                 w0,   inc_w,
                 buff_1,
                 alpha11 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 U20, rs_U, cs_U,
                 w0,  inc_w,
                 buff_1,
                 a21, rs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Househ2_UT( FLA_LEFT,
      //                 a21_t,
      //                 a21_b, tau11 );
      FLA_Househ2_UT_l_opz( m_ahead - 1,
                            a21_t,
                            a21_b, rs_A,
                            tau11 );

      // FLA_Copy( a21, u21 );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u21, rs_U );

      // FLA_Set( FLA_ONE, u21_t );
      *u21_t = *buff_1;

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, u21, FLA_ZERO, z01 );
      // FLA_Dot( a12t, u21, zeta11 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, u21, FLA_ZERO, z21 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 u21, rs_U,
                 buff_0,
                 z01, rs_Z );
      bl1_zdot( BLIS1_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                u21,  rs_U,
                zeta11 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 u21, rs_U,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, U20, u21, FLA_ZERO, t01 );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 U20, rs_U, cs_U,
                 u21, rs_U,
                 buff_0,
                 t01, rs_T );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}

