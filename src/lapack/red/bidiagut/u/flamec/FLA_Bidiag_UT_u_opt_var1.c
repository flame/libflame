/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_opt_var1( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  return FLA_Bidiag_UT_u_step_opt_var1( A, TU, TV );
}

FLA_Error FLA_Bidiag_UT_u_step_opt_var1( FLA_Obj A, FLA_Obj T, FLA_Obj S )
{
  FLA_Datatype datatype;
  int          m_A, n_A, m_TS;
  int          rs_A, cs_A;
  int          rs_T, cs_T;
  int          rs_S, cs_S;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  m_TS     = FLA_Obj_length( T );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );
  
  rs_S     = FLA_Obj_row_stride( S );
  cs_S     = FLA_Obj_col_stride( S );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_T = FLA_FLOAT_PTR( T );
      float* buff_S = FLA_FLOAT_PTR( S );

      FLA_Bidiag_UT_u_step_ops_var1( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_T = FLA_DOUBLE_PTR( T );
      double* buff_S = FLA_DOUBLE_PTR( S );

      FLA_Bidiag_UT_u_step_opd_var1( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );
      scomplex* buff_S = FLA_COMPLEX_PTR( S );

      FLA_Bidiag_UT_u_step_opc_var1( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );
      dcomplex* buff_S = FLA_DOUBLE_COMPLEX_PTR( S );

      FLA_Bidiag_UT_u_step_opz_var1( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_ops_var1( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );

  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  float*    buff_v = ( float* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  int       inc_v  = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    float*    a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    float*    tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    float*    s01      = buff_S + (i  )*cs_S + (0  )*rs_S;
    float*    sigma11  = buff_S + (i  )*cs_S + (i  )*rs_S;

    float*    v21      = buff_v + (i+1)*inc_v;

    float*    a12t_l   = a12t   + (0  )*cs_A + (0  )*rs_A;
    float*    a12t_r   = a12t   + (1  )*cs_A + (0  )*rs_A;

    float*    A22_l    = A22    + (0  )*cs_A + (0  )*rs_A;
    float*    A22_r    = A22    + (1  )*cs_A + (0  )*rs_A;

    float*    v21_t    = v21    + (0  )*inc_v;
    float*    v21_b    = v21    + (1  )*inc_v;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    FLA_Househ2_UT_l_ops( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );

    if ( n_ahead > 0 )
    {
      // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21, a12t, A22 );
      FLA_Apply_H2_UT_l_ops_var1( m_ahead,
                                  n_ahead,
                                  tau11,
                                  a21,  rs_A,
                                  a12t, cs_A,
                                  A22,  rs_A, cs_A );

      // FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );
      FLA_Househ2_UT_r_ops( n_ahead - 1,
                            a12t_l,
                            a12t_r, cs_A,
                            sigma11 );

      // FLA_Set( FLA_ONE, v21_t );
      // FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );
      *v21_t = *buff_1;
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  n_ahead - 1,
                  a12t_r, cs_A,
                  v21_b,  inc_v );

      // FLA_Apply_H2_UT( FLA_RIGHT, sigma11, v21_b, A22_l, A22_r );
      FLA_Apply_H2_UT_r_ops_var1( m_ahead,
                                  n_ahead - 1,
                                  sigma11,
                                  v21_b, inc_v,
                                  A22_l, rs_A,
                                  A22_r, rs_A, cs_A );

      // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, s01 );
      bl1_sgemv( BLIS1_CONJ_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 v21, inc_v,
                 buff_0,
                 s01, rs_S );
    }

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_scopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01,  rs_T );
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

  // FLA_Obj_free( &v );
  FLA_free( buff_v );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_opd_var1( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );

  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  double*   buff_v = ( double* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  int       inc_v  = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    double*   a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    double*   tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    double*   s01      = buff_S + (i  )*cs_S + (0  )*rs_S;
    double*   sigma11  = buff_S + (i  )*cs_S + (i  )*rs_S;

    double*   v21      = buff_v + (i+1)*inc_v;

    double*   a12t_l   = a12t   + (0  )*cs_A + (0  )*rs_A;
    double*   a12t_r   = a12t   + (1  )*cs_A + (0  )*rs_A;

    double*   A22_l    = A22    + (0  )*cs_A + (0  )*rs_A;
    double*   A22_r    = A22    + (1  )*cs_A + (0  )*rs_A;

    double*   v21_t    = v21    + (0  )*inc_v;
    double*   v21_b    = v21    + (1  )*inc_v;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    FLA_Househ2_UT_l_opd( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );

    if ( n_ahead > 0 )
    {
      // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21, a12t, A22 );
      FLA_Apply_H2_UT_l_opd_var1( m_ahead,
                                  n_ahead,
                                  tau11,
                                  a21,  rs_A,
                                  a12t, cs_A,
                                  A22,  rs_A, cs_A );

      // FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );
      FLA_Househ2_UT_r_opd( n_ahead - 1,
                            a12t_l,
                            a12t_r, cs_A,
                            sigma11 );

      // FLA_Set( FLA_ONE, v21_t );
      // FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );
      *v21_t = *buff_1;
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  n_ahead - 1,
                  a12t_r, cs_A,
                  v21_b,  inc_v );

      // FLA_Apply_H2_UT( FLA_RIGHT, sigma11, v21_b, A22_l, A22_r );
      FLA_Apply_H2_UT_r_opd_var1( m_ahead,
                                  n_ahead - 1,
                                  sigma11,
                                  v21_b, inc_v,
                                  A22_l, rs_A,
                                  A22_r, rs_A, cs_A );

      // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, s01 );
      bl1_dgemv( BLIS1_CONJ_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 v21, inc_v,
                 buff_0,
                 s01, rs_S );
    }

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_dcopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01,  rs_T );
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

  // FLA_Obj_free( &v );
  FLA_free( buff_v );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_opc_var1( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );

  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  scomplex* buff_v = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  int       inc_v  = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    scomplex* s01      = buff_S + (i  )*cs_S + (0  )*rs_S;
    scomplex* sigma11  = buff_S + (i  )*cs_S + (i  )*rs_S;

    scomplex* v21      = buff_v + (i+1)*inc_v;

    scomplex* a12t_l   = a12t   + (0  )*cs_A + (0  )*rs_A;
    scomplex* a12t_r   = a12t   + (1  )*cs_A + (0  )*rs_A;

    scomplex* A22_l    = A22    + (0  )*cs_A + (0  )*rs_A;
    scomplex* A22_r    = A22    + (1  )*cs_A + (0  )*rs_A;

    scomplex* v21_t    = v21    + (0  )*inc_v;
    scomplex* v21_b    = v21    + (1  )*inc_v;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    FLA_Househ2_UT_l_opc( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );

    if ( n_ahead > 0 )
    {
      // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21, a12t, A22 );
      FLA_Apply_H2_UT_l_opc_var1( m_ahead,
                                  n_ahead,
                                  tau11,
                                  a21,  rs_A,
                                  a12t, cs_A,
                                  A22,  rs_A, cs_A );

      // FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );
      FLA_Househ2_UT_r_opc( n_ahead - 1,
                            a12t_l,
                            a12t_r, cs_A,
                            sigma11 );

      // FLA_Set( FLA_ONE, v21_t );
      // FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );
      *v21_t = *buff_1;
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  n_ahead - 1,
                  a12t_r, cs_A,
                  v21_b,  inc_v );

      // FLA_Apply_H2_UT( FLA_RIGHT, sigma11, v21_b, A22_l, A22_r );
      FLA_Apply_H2_UT_r_opc_var1( m_ahead,
                                  n_ahead - 1,
                                  sigma11,
                                  v21_b, inc_v,
                                  A22_l, rs_A,
                                  A22_r, rs_A, cs_A );

      // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, s01 );
      bl1_cgemv( BLIS1_CONJ_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 v21, inc_v,
                 buff_0,
                 s01, rs_S );
    }

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_ccopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01,  rs_T );
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

  // FLA_Obj_free( &v );
  FLA_free( buff_v );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_opz_var1( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );

  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  dcomplex* buff_v = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  int       inc_v  = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    dcomplex* s01      = buff_S + (i  )*cs_S + (0  )*rs_S;
    dcomplex* sigma11  = buff_S + (i  )*cs_S + (i  )*rs_S;

    dcomplex* v21      = buff_v + (i+1)*inc_v;

    dcomplex* a12t_l   = a12t   + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a12t_r   = a12t   + (1  )*cs_A + (0  )*rs_A;

    dcomplex* A22_l    = A22    + (0  )*cs_A + (0  )*rs_A;
    dcomplex* A22_r    = A22    + (1  )*cs_A + (0  )*rs_A;

    dcomplex* v21_t    = v21    + (0  )*inc_v;
    dcomplex* v21_b    = v21    + (1  )*inc_v;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    FLA_Househ2_UT_l_opz( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );

    if ( n_ahead > 0 )
    {
      // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21, a12t, A22 );
      FLA_Apply_H2_UT_l_opz_var1( m_ahead,
                                  n_ahead,
                                  tau11,
                                  a21,  rs_A,
                                  a12t, cs_A,
                                  A22,  rs_A, cs_A );

      // FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );
      FLA_Househ2_UT_r_opz( n_ahead - 1,
                            a12t_l,
                            a12t_r, cs_A,
                            sigma11 );

      // FLA_Set( FLA_ONE, v21_t );
      // FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );
      *v21_t = *buff_1;
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  n_ahead - 1,
                  a12t_r, cs_A,
                  v21_b,  inc_v );

      // FLA_Apply_H2_UT( FLA_RIGHT, sigma11, v21_b, A22_l, A22_r );
      FLA_Apply_H2_UT_r_opz_var1( m_ahead,
                                  n_ahead - 1,
                                  sigma11,
                                  v21_b, inc_v,
                                  A22_l, rs_A,
                                  A22_r, rs_A, cs_A );

      // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, s01 );
      bl1_zgemv( BLIS1_CONJ_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 v21, inc_v,
                 buff_0,
                 s01, rs_S );
    }

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_zcopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01,  rs_T );
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

  // FLA_Obj_free( &v );
  FLA_free( buff_v );

  return FLA_SUCCESS;
}

