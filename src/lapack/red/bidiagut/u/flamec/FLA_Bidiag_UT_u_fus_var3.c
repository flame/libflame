/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_ofu_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  return FLA_Bidiag_UT_u_step_ofu_var3( A, TU, TV );
}

FLA_Error FLA_Bidiag_UT_u_step_ofu_var3( FLA_Obj A, FLA_Obj T, FLA_Obj S )
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

      FLA_Bidiag_UT_u_step_ofs_var3( m_A,
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

      FLA_Bidiag_UT_u_step_ofd_var3( m_A,
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

      FLA_Bidiag_UT_u_step_ofc_var3( m_A,
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

      FLA_Bidiag_UT_u_step_ofz_var3( m_A,
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



FLA_Error FLA_Bidiag_UT_u_step_ofs_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );

  float     alpha12;
  float     minus_conj_alpha12;
  float     psi11_minus_alpha12;
  float     minus_inv_tau11;
  float     minus_upsilon11;
  float     minus_conj_nu11;
  float     minus_conj_psi11;
  float     minus_zeta11;
  float     beta;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &ap );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &up );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  float*    buff_w  = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_ap = ( float* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_u  = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_up = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_v  = ( float* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_y  = ( float* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_z  = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_w   = 1;
  int       inc_ap  = 1;
  int       inc_u   = 1;
  int       inc_up  = 1;
  int       inc_v   = 1;
  int       inc_y   = 1;
  int       inc_z   = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    float*    a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    float*    A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    float*    alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    float*    a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    float*    A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    float*    A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    float*    t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    float*    tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    float*    s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    float*    sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    float*    w21       = buff_w  + (i+1)*inc_w;

    float*    a12p      = buff_ap + (i+1)*inc_ap;

    float*    upsilon11 = buff_u  + (i  )*inc_u;
    float*    u21       = buff_u  + (i+1)*inc_u;

    float*    u21p      = buff_up + (i+1)*inc_up;

    float*    nu11      = buff_v  + (i  )*inc_v;
    float*    v21       = buff_v  + (i+1)*inc_v;

    float*    psi11     = buff_y  + (i  )*inc_y;
    float*    y21       = buff_y  + (i+1)*inc_y;

    float*    zeta11    = buff_z  + (i  )*inc_z;
    float*    z21       = buff_z  + (i+1)*inc_z;

    float*    a12p_t    = a12p    + (0  )*inc_ap;
    float*    a12p_b    = a12p    + (1  )*inc_ap;

    float*    v21_t     = v21     + (0  )*inc_v;
    float*    v21_b     = v21     + (1  )*inc_v;

    float*    a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    float*    a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    float*    A22_l     = A22     + (0  )*cs_A + (0  )*rs_A;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       m_behind  = i;
    int       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon11, minus_upsilon11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon11 );
      bl1_smult3( buff_m1, upsilon11, &minus_upsilon11 );

      // FLA_Copy( zeta11, minus_zeta11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta11 );
      bl1_smult3( buff_m1, zeta11, &minus_zeta11 );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, psi11, minus_conj_psi11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_psi11 );
      bl1_scopyconj( psi11, &minus_conj_psi11 );
      bl1_sscals( buff_m1, &minus_conj_psi11 );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, nu11, minus_conj_nu11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_nu11 );
      bl1_scopyconj( nu11, &minus_conj_nu11 );
      bl1_sscals( buff_m1, &minus_conj_nu11 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_upsilon11, psi11, alpha11 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_zeta11,    nu11,  alpha11 );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  1,
                  &minus_upsilon11,
                  psi11,   1,
                  alpha11, 1 );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  1,
                  &minus_zeta11,
                  nu11,    1,
                  alpha11, 1 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi11, u21, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_nu11,  z21, a21 );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_psi11,
                  u21, inc_u,
                  a21, rs_A );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_nu11,
                  z21, inc_z,
                  a21, rs_A );

      // FLA_Axpyt( FLA_TRANSPOSE, minus_upsilon11, y21, a12t );
      // FLA_Axpyt( FLA_TRANSPOSE, minus_zeta11,    v21, a12t );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  &minus_upsilon11,
                  y21,  inc_y,
                  a12t, cs_A );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  &minus_zeta11,
                  v21,  inc_v,
                  a12t, cs_A );
    }

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    // FLA_Copy( a21, u21p );
    FLA_Househ2_UT_l_ops( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );
    bl1_scopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                a21,  rs_A,
                u21p, inc_up );

    if ( n_ahead > 0 )
    {
      // FLA_Copy( FLA_MINUS_ONE, minus_inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, minus_inv_tau11 );
      bl1_sdiv3( buff_m1, tau11, &minus_inv_tau11 );

      // FLA_Copyt( FLA_TRANSPOSE, a12t, a12p );
      // FLA_Axpyt( FLA_TRANSPOSE, minus_inv_tau11, a12t, a12p );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  a12p, inc_ap );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  &minus_inv_tau11,
                  a12t, cs_A,
                  a12p, inc_ap );
    }

    if ( m_behind > 0 && n_ahead > 0 )
    {
      // FLA_Ger( FLA_MINUS_ONE, u21, y21, A22 );
      // FLA_Ger( FLA_MINUS_ONE, z21, v21, A22 );
      // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, u21p, FLA_ZERO, y21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      FLA_Fused_Gerc2_Ahx_Axpy_Ax_ops_var1( m_ahead,
                                            n_ahead,
                                            tau11,
                                            buff_m1,
                                            u21, inc_u,
                                            y21, inc_y,
                                            z21, inc_z,
                                            v21, inc_v,
                                            A22, rs_A, cs_A,
                                            u21p, inc_up,
                                            a12p, inc_ap,
                                            w21,  inc_w );
                                           
                                           
    }
    else if ( n_ahead > 0 )
    {
      // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, u21p, FLA_ZERO, y21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      FLA_Fused_Ahx_Axpy_Ax_ops_var1( m_ahead,
                                      n_ahead,
                                      tau11,
                                      buff_0,
                                      A22,  rs_A, cs_A,
                                      u21p, inc_up,
                                      a12p, inc_ap,
                                      y21,  inc_y,
                                      w21,  inc_w );
    }

    if ( n_ahead > 0 )
    {
      // FLA_Axpyt( FLA_TRANSPOSE, FLA_ONE, a12t, y21 );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  buff_1,
                  a12t, cs_A,
                  y21,  inc_y );

      // FLA_Househ2s_UT( FLA_RIGHT,
      //                  a12p_t,
      //                  a12p_b,
      //                  alpha12, psi11_minus_alpha12, sigma11 );
      FLA_Househ2s_UT_r_ops( n_ahead - 1,
                             a12p_t,
                             a12p_b, inc_ap,
                             &alpha12,
                             &psi11_minus_alpha12,
                             sigma11 );

      // FLA_Copy( a12p, v21 );
      // FLA_Mult_add( FLA_MINUS_ONE, alpha12, v21_t );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, psi11_minus_alpha12, v21 );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12p, inc_ap,
                  v21,  inc_v );
      bl1_smult4( buff_m1, &alpha12, v21_t, v21_t );
      bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                     n_ahead,
                     &psi11_minus_alpha12,
                     v21, inc_v );

      // FLA_Copy( alpha12, a12t_l );
      // FLA_Copyt( FLA_TRANSPOSE, v21_b, a12t_r );
      *a12t_l = alpha12;
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  n_ahead - 1,
                  v21_b,  inc_v,
                  a12t_r, cs_A );
    }

    // FLA_Copy( u21p, u21 );
    bl1_scopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                u21p, inc_up,
                u21,  inc_u );

    if ( n_ahead > 0 )
    {
      // FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, beta );
      bl1_sdot( BLIS1_CONJUGATE,
                n_ahead,
                y21, inc_y,
                v21, inc_v,
                &beta );
      bl1_sscals( &minus_inv_tau11, &beta );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha12, minus_conj_alpha12 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_alpha12 );
      bl1_scopyconj( &alpha12, &minus_conj_alpha12 );
      bl1_sneg1( &minus_conj_alpha12 );

      // FLA_Copy( w21, z21 );
      // FLA_Axpy( minus_conj_alpha12, A22_l, z21 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, z21 );
      // FLA_Axpy( beta, u21, z21 );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, inc_z );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_alpha12,
                  A22_l, rs_A,
                  z21,   inc_z );
      bl1_sinvscalv( BLIS1_CONJUGATE,
                     m_ahead,
                     &psi11_minus_alpha12,
                     z21, inc_z );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  u21, inc_u,
                  z21, inc_z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11,   y21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_sinvscalv( BLIS1_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, inc_y );
      bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, inc_z );

      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A02, v21, FLA_ZERO, s01 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 v21, inc_v,
                 buff_0,
                 s01, rs_S );
    }

    // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ONE, t01 );
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
               u21, inc_u,
               buff_1,
               t01, rs_T );

    if ( m_behind + 1 == b_alg && n_ahead > 0 )
    {
      // FLA_Ger( FLA_MINUS_ONE, u21, y21, A22 );
      // FLA_Ger( FLA_MINUS_ONE, z21, v21, A22 );
      bl1_sger( BLIS1_NO_CONJUGATE,
                BLIS1_NO_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u21, inc_u,
                y21, inc_y,
                A22, rs_A, cs_A );
      bl1_sger( BLIS1_NO_CONJUGATE,
                BLIS1_NO_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z21, inc_z,
                v21, inc_v,
                A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  // FLA_Obj_free( &ap );
  // FLA_Obj_free( &u );
  // FLA_Obj_free( &up );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_w );
  FLA_free( buff_ap );
  FLA_free( buff_u );
  FLA_free( buff_up );
  FLA_free( buff_v );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_ofd_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );

  double    alpha12;
  double    minus_conj_alpha12;
  double    psi11_minus_alpha12;
  double    minus_inv_tau11;
  double    minus_upsilon11;
  double    minus_conj_nu11;
  double    minus_conj_psi11;
  double    minus_zeta11;
  double    beta;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &ap );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &up );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  double*   buff_w  = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_ap = ( double* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_u  = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_up = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_v  = ( double* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_y  = ( double* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_z  = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_w   = 1;
  int       inc_ap  = 1;
  int       inc_u   = 1;
  int       inc_up  = 1;
  int       inc_v   = 1;
  int       inc_y   = 1;
  int       inc_z   = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    double*   a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    double*   A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    double*   alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    double*   a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    double*   A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    double*   A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    double*   t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    double*   tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    double*   s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    double*   sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    double*   w21       = buff_w  + (i+1)*inc_w;

    double*   a12p      = buff_ap + (i+1)*inc_ap;

    double*   upsilon11 = buff_u  + (i  )*inc_u;
    double*   u21       = buff_u  + (i+1)*inc_u;

    double*   u21p      = buff_up + (i+1)*inc_up;

    double*   nu11      = buff_v  + (i  )*inc_v;
    double*   v21       = buff_v  + (i+1)*inc_v;

    double*   psi11     = buff_y  + (i  )*inc_y;
    double*   y21       = buff_y  + (i+1)*inc_y;

    double*   zeta11    = buff_z  + (i  )*inc_z;
    double*   z21       = buff_z  + (i+1)*inc_z;

    double*   a12p_t    = a12p    + (0  )*inc_ap;
    double*   a12p_b    = a12p    + (1  )*inc_ap;

    double*   v21_t     = v21     + (0  )*inc_v;
    double*   v21_b     = v21     + (1  )*inc_v;

    double*   a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    double*   a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    double*   A22_l     = A22     + (0  )*cs_A + (0  )*rs_A;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       m_behind  = i;
    int       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon11, minus_upsilon11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon11 );
      bl1_dmult3( buff_m1, upsilon11, &minus_upsilon11 );

      // FLA_Copy( zeta11, minus_zeta11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta11 );
      bl1_dmult3( buff_m1, zeta11, &minus_zeta11 );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, psi11, minus_conj_psi11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_psi11 );
      bl1_dcopyconj( psi11, &minus_conj_psi11 );
      bl1_dscals( buff_m1, &minus_conj_psi11 );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, nu11, minus_conj_nu11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_nu11 );
      bl1_dcopyconj( nu11, &minus_conj_nu11 );
      bl1_dscals( buff_m1, &minus_conj_nu11 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi11, upsilon11, alpha11 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_nu11,  zeta11,    alpha11 );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  1,
                  &minus_conj_psi11,
                  upsilon11, 1,
                  alpha11,   1 );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  1,
                  &minus_conj_nu11,
                  zeta11,  1,
                  alpha11, 1 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi11, u21, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_nu11,  z21, a21 );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_psi11,
                  u21, inc_u,
                  a21, rs_A );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_nu11,
                  z21, inc_z,
                  a21, rs_A );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_upsilon11, y21, a12t );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_zeta11,    v21, a12t );
      bl1_daxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  &minus_upsilon11,
                  y21,  inc_y,
                  a12t, cs_A );
      bl1_daxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  &minus_zeta11,
                  v21,  inc_v,
                  a12t, cs_A );
    }

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    // FLA_Copy( a21, u21p );
    FLA_Househ2_UT_l_opd( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );
    bl1_dcopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                a21,  rs_A,
                u21p, inc_up );

    if ( n_ahead > 0 )
    {
      // FLA_Copy( FLA_MINUS_ONE, minus_inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, minus_inv_tau11 );
      bl1_ddiv3( buff_m1, tau11, &minus_inv_tau11 );

      // FLA_Copyt( FLA_TRANSPOSE, a12t, a12p );
      // FLA_Axpyt( FLA_TRANSPOSE, minus_inv_tau11, a12t, a12p );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  a12p, inc_ap );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  &minus_inv_tau11,
                  a12t, cs_A,
                  a12p, inc_ap );
    }

    if ( m_behind > 0 && n_ahead > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u21, y21, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );
      // FLA_Gemvc( FLA_CONJ_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, u21p, FLA_ZERO, y21 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      FLA_Fused_Gerc2_Ahx_Axpy_Ax_opd_var1( m_ahead,
                                            n_ahead,
                                            tau11,
                                            buff_m1,
                                            u21, inc_u,
                                            y21, inc_y,
                                            z21, inc_z,
                                            v21, inc_v,
                                            A22, rs_A, cs_A,
                                            u21p, inc_up,
                                            a12p, inc_ap,
                                            w21,  inc_w );
                                           
                                           
    }
    else if ( n_ahead > 0 )
    {
      // FLA_Gemvc( FLA_CONJ_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, u21p, FLA_ZERO, y21 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      FLA_Fused_Ahx_Axpy_Ax_opd_var1( m_ahead,
                                      n_ahead,
                                      tau11,
                                      buff_0,
                                      A22,  rs_A, cs_A,
                                      u21p, inc_up,
                                      a12p, inc_ap,
                                      y21,  inc_y,
                                      w21,  inc_w );
    }

    if ( n_ahead > 0 )
    {
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_ONE, a12t, y21 );
      bl1_daxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  buff_1,
                  a12t, cs_A,
                  y21,  inc_y );

      // FLA_Househ2s_UT( FLA_RIGHT,
      //                  a12p_t,
      //                  a12p_b,
      //                  alpha12, psi11_minus_alpha12, sigma11 );
      FLA_Househ2s_UT_r_opd( n_ahead - 1,
                             a12p_t,
                             a12p_b, inc_ap,
                             &alpha12,
                             &psi11_minus_alpha12,
                             sigma11 );

      // FLA_Copy( a12p, v21 );
      // FLA_Mult_add( FLA_MINUS_ONE, alpha12, v21_t );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, psi11_minus_alpha12, v21 );
      // FLA_Conjugate( v21_b );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12p, inc_ap,
                  v21,  inc_v );
      bl1_dmult4( buff_m1, &alpha12, v21_t, v21_t );
      bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                     n_ahead,
                     &psi11_minus_alpha12,
                     v21, inc_v );
      bl1_dconjv( n_ahead - 1,
                  v21_b, inc_v );

      // FLA_Copyt( FLA_NO_TRANSPOSE, alpha12, a12t_l );
      // FLA_Copyt( FLA_TRANSPOSE, v21_b, a12t_r );
      *a12t_l = alpha12;
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  n_ahead - 1,
                  v21_b,  inc_v,
                  a12t_r, cs_A );
    }

    // FLA_Copy( u21p, u21 );
    bl1_dcopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                u21p, inc_up,
                u21,  inc_u );

    if ( n_ahead > 0 )
    {
      // FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, beta );
      bl1_ddot( BLIS1_CONJUGATE,
                n_ahead,
                y21, inc_y,
                v21, inc_v,
                &beta );
      bl1_dscals( &minus_inv_tau11, &beta );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha12, minus_conj_alpha12 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_alpha12 );
      bl1_dcopyconj( &alpha12, &minus_conj_alpha12 );
      bl1_dneg1( &minus_conj_alpha12 );

      // FLA_Copy( w21, z21 );
      // FLA_Axpy( minus_conj_alpha12, A22_l, z21 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, z21 );
      // FLA_Axpy( beta, u21, z21 );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, inc_z );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_alpha12,
                  A22_l, rs_A,
                  z21,   inc_z );
      bl1_dinvscalv( BLIS1_CONJUGATE,
                     m_ahead,
                     &psi11_minus_alpha12,
                     z21, inc_z );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  u21, inc_u,
                  z21, inc_z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11,   y21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_dinvscalv( BLIS1_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, inc_y );
      bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, inc_z );

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

    // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ONE, t01 );
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
               u21, inc_u,
               buff_1,
               t01, rs_T );

    if ( m_behind + 1 == b_alg && n_ahead > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u21, y21, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );
      bl1_dger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u21, inc_u,
                y21, inc_y,
                A22, rs_A, cs_A );
      bl1_dger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z21, inc_z,
                v21, inc_v,
                A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  // FLA_Obj_free( &ap );
  // FLA_Obj_free( &u );
  // FLA_Obj_free( &up );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_w );
  FLA_free( buff_ap );
  FLA_free( buff_u );
  FLA_free( buff_up );
  FLA_free( buff_v );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_ofc_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );

  scomplex  alpha12;
  scomplex  minus_conj_alpha12;
  scomplex  psi11_minus_alpha12;
  scomplex  minus_inv_tau11;
  scomplex  minus_upsilon11;
  scomplex  minus_conj_nu11;
  scomplex  minus_conj_psi11;
  scomplex  minus_zeta11;
  scomplex  beta;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &ap );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &up );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  scomplex* buff_w  = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_ap = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_u  = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_up = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_v  = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_y  = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_z  = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_w   = 1;
  int       inc_ap  = 1;
  int       inc_u   = 1;
  int       inc_up  = 1;
  int       inc_v   = 1;
  int       inc_y   = 1;
  int       inc_z   = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    scomplex* alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    scomplex* s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    scomplex* sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    scomplex* w21       = buff_w  + (i+1)*inc_w;

    scomplex* a12p      = buff_ap + (i+1)*inc_ap;

    scomplex* upsilon11 = buff_u  + (i  )*inc_u;
    scomplex* u21       = buff_u  + (i+1)*inc_u;

    scomplex* u21p      = buff_up + (i+1)*inc_up;

    scomplex* nu11      = buff_v  + (i  )*inc_v;
    scomplex* v21       = buff_v  + (i+1)*inc_v;

    scomplex* psi11     = buff_y  + (i  )*inc_y;
    scomplex* y21       = buff_y  + (i+1)*inc_y;

    scomplex* zeta11    = buff_z  + (i  )*inc_z;
    scomplex* z21       = buff_z  + (i+1)*inc_z;

    scomplex* a12p_t    = a12p    + (0  )*inc_ap;
    scomplex* a12p_b    = a12p    + (1  )*inc_ap;

    scomplex* v21_t     = v21     + (0  )*inc_v;
    scomplex* v21_b     = v21     + (1  )*inc_v;

    scomplex* a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    scomplex* a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    scomplex* A22_l     = A22     + (0  )*cs_A + (0  )*rs_A;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       m_behind  = i;
    int       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon11, minus_upsilon11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon11 );
      bl1_cmult3( buff_m1, upsilon11, &minus_upsilon11 );

      // FLA_Copy( zeta11, minus_zeta11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta11 );
      bl1_cmult3( buff_m1, zeta11, &minus_zeta11 );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, psi11, minus_conj_psi11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_psi11 );
      bl1_ccopyconj( psi11, &minus_conj_psi11 );
      bl1_cscals( buff_m1, &minus_conj_psi11 );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, nu11, minus_conj_nu11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_nu11 );
      bl1_ccopyconj( nu11, &minus_conj_nu11 );
      bl1_cscals( buff_m1, &minus_conj_nu11 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi11, upsilon11, alpha11 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_nu11,  zeta11,    alpha11 );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  1,
                  &minus_conj_psi11,
                  upsilon11, 1,
                  alpha11,   1 );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  1,
                  &minus_conj_nu11,
                  zeta11,  1,
                  alpha11, 1 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi11, u21, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_nu11,  z21, a21 );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_psi11,
                  u21, inc_u,
                  a21, rs_A );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_nu11,
                  z21, inc_z,
                  a21, rs_A );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_upsilon11, y21, a12t );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_zeta11,    v21, a12t );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  &minus_upsilon11,
                  y21,  inc_y,
                  a12t, cs_A );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  &minus_zeta11,
                  v21,  inc_v,
                  a12t, cs_A );
    }

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    // FLA_Copy( a21, u21p );
    FLA_Househ2_UT_l_opc( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );
    bl1_ccopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                a21,  rs_A,
                u21p, inc_up );

    if ( n_ahead > 0 )
    {
      // FLA_Copy( FLA_MINUS_ONE, minus_inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, minus_inv_tau11 );
      bl1_cdiv3( buff_m1, tau11, &minus_inv_tau11 );

      // FLA_Copyt( FLA_TRANSPOSE, a12t, a12p );
      // FLA_Axpyt( FLA_TRANSPOSE, minus_inv_tau11, a12t, a12p );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  a12p, inc_ap );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  &minus_inv_tau11,
                  a12t, cs_A,
                  a12p, inc_ap );
    }

    if ( m_behind > 0 && n_ahead > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u21, y21, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );
      // FLA_Gemvc( FLA_CONJ_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, u21p, FLA_ZERO, y21 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      FLA_Fused_Gerc2_Ahx_Axpy_Ax_opc_var1( m_ahead,
                                            n_ahead,
                                            tau11,
                                            buff_m1,
                                            u21, inc_u,
                                            y21, inc_y,
                                            z21, inc_z,
                                            v21, inc_v,
                                            A22, rs_A, cs_A,
                                            u21p, inc_up,
                                            a12p, inc_ap,
                                            w21,  inc_w );
                                           
                                           
    }
    else if ( n_ahead > 0 )
    {
      // FLA_Gemvc( FLA_CONJ_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, u21p, FLA_ZERO, y21 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      FLA_Fused_Ahx_Axpy_Ax_opc_var1( m_ahead,
                                      n_ahead,
                                      tau11,
                                      buff_0,
                                      A22,  rs_A, cs_A,
                                      u21p, inc_up,
                                      a12p, inc_ap,
                                      y21,  inc_y,
                                      w21,  inc_w );
    }

    if ( n_ahead > 0 )
    {
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_ONE, a12t, y21 );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  buff_1,
                  a12t, cs_A,
                  y21,  inc_y );

      // FLA_Househ2s_UT( FLA_RIGHT,
      //                  a12p_t,
      //                  a12p_b,
      //                  alpha12, psi11_minus_alpha12, sigma11 );
      FLA_Househ2s_UT_r_opc( n_ahead - 1,
                             a12p_t,
                             a12p_b, inc_ap,
                             &alpha12,
                             &psi11_minus_alpha12,
                             sigma11 );

      // FLA_Copy( a12p, v21 );
      // FLA_Mult_add( FLA_MINUS_ONE, alpha12, v21_t );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, psi11_minus_alpha12, v21 );
      // FLA_Conjugate( v21_b );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12p, inc_ap,
                  v21,  inc_v );
      bl1_cmult4( buff_m1, &alpha12, v21_t, v21_t );
      bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                     n_ahead,
                     &psi11_minus_alpha12,
                     v21, inc_v );
      bl1_cconjv( n_ahead - 1,
                  v21_b, inc_v );

      // FLA_Copyt( FLA_NO_TRANSPOSE, alpha12, a12t_l );
      // FLA_Copyt( FLA_TRANSPOSE, v21_b, a12t_r );
      *a12t_l = alpha12;
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  n_ahead - 1,
                  v21_b,  inc_v,
                  a12t_r, cs_A );
    }

    // FLA_Copy( u21p, u21 );
    bl1_ccopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                u21p, inc_up,
                u21,  inc_u );

    if ( n_ahead > 0 )
    {
      // FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, beta );
      bl1_cdot( BLIS1_CONJUGATE,
                n_ahead,
                y21, inc_y,
                v21, inc_v,
                &beta );
      bl1_cscals( &minus_inv_tau11, &beta );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha12, minus_conj_alpha12 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_alpha12 );
      bl1_ccopyconj( &alpha12, &minus_conj_alpha12 );
      bl1_cneg1( &minus_conj_alpha12 );

      // FLA_Copy( w21, z21 );
      // FLA_Axpy( minus_conj_alpha12, A22_l, z21 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, z21 );
      // FLA_Axpy( beta, u21, z21 );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, inc_z );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_alpha12,
                  A22_l, rs_A,
                  z21,   inc_z );
      bl1_cinvscalv( BLIS1_CONJUGATE,
                     m_ahead,
                     &psi11_minus_alpha12,
                     z21, inc_z );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  u21, inc_u,
                  z21, inc_z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11,   y21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_cinvscalv( BLIS1_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, inc_y );
      bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, inc_z );

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

    // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ONE, t01 );
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
               u21, inc_u,
               buff_1,
               t01, rs_T );

    if ( m_behind + 1 == b_alg && n_ahead > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u21, y21, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );
      bl1_cger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u21, inc_u,
                y21, inc_y,
                A22, rs_A, cs_A );
      bl1_cger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z21, inc_z,
                v21, inc_v,
                A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  // FLA_Obj_free( &ap );
  // FLA_Obj_free( &u );
  // FLA_Obj_free( &up );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_w );
  FLA_free( buff_ap );
  FLA_free( buff_u );
  FLA_free( buff_up );
  FLA_free( buff_v );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_ofz_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );

  dcomplex  alpha12;
  dcomplex  minus_conj_alpha12;
  dcomplex  psi11_minus_alpha12;
  dcomplex  minus_inv_tau11;
  dcomplex  minus_upsilon11;
  dcomplex  minus_conj_nu11;
  dcomplex  minus_conj_psi11;
  dcomplex  minus_zeta11;
  dcomplex  beta;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &ap );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &up );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  dcomplex* buff_w  = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_ap = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_u  = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_up = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_v  = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_y  = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_z  = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_w   = 1;
  int       inc_ap  = 1;
  int       inc_u   = 1;
  int       inc_up  = 1;
  int       inc_v   = 1;
  int       inc_y   = 1;
  int       inc_z   = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    dcomplex* s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    dcomplex* sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    dcomplex* w21       = buff_w  + (i+1)*inc_w;

    dcomplex* a12p      = buff_ap + (i+1)*inc_ap;

    dcomplex* upsilon11 = buff_u  + (i  )*inc_u;
    dcomplex* u21       = buff_u  + (i+1)*inc_u;

    dcomplex* u21p      = buff_up + (i+1)*inc_up;

    dcomplex* nu11      = buff_v  + (i  )*inc_v;
    dcomplex* v21       = buff_v  + (i+1)*inc_v;

    dcomplex* psi11     = buff_y  + (i  )*inc_y;
    dcomplex* y21       = buff_y  + (i+1)*inc_y;

    dcomplex* zeta11    = buff_z  + (i  )*inc_z;
    dcomplex* z21       = buff_z  + (i+1)*inc_z;

    dcomplex* a12p_t    = a12p    + (0  )*inc_ap;
    dcomplex* a12p_b    = a12p    + (1  )*inc_ap;

    dcomplex* v21_t     = v21     + (0  )*inc_v;
    dcomplex* v21_b     = v21     + (1  )*inc_v;

    dcomplex* a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    dcomplex* A22_l     = A22     + (0  )*cs_A + (0  )*rs_A;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       m_behind  = i;
    int       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon11, minus_upsilon11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon11 );
      bl1_zmult3( buff_m1, upsilon11, &minus_upsilon11 );

      // FLA_Copy( zeta11, minus_zeta11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta11 );
      bl1_zmult3( buff_m1, zeta11, &minus_zeta11 );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, psi11, minus_conj_psi11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_psi11 );
      bl1_zcopyconj( psi11, &minus_conj_psi11 );
      bl1_zscals( buff_m1, &minus_conj_psi11 );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, nu11, minus_conj_nu11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_nu11 );
      bl1_zcopyconj( nu11, &minus_conj_nu11 );
      bl1_zscals( buff_m1, &minus_conj_nu11 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi11, upsilon11, alpha11 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_nu11,  zeta11,    alpha11 );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  1,
                  &minus_conj_psi11,
                  upsilon11, 1,
                  alpha11,   1 );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  1,
                  &minus_conj_nu11,
                  zeta11,  1,
                  alpha11, 1 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi11, u21, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_nu11,  z21, a21 );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_psi11,
                  u21, inc_u,
                  a21, rs_A );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_nu11,
                  z21, inc_z,
                  a21, rs_A );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_upsilon11, y21, a12t );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_zeta11,    v21, a12t );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  &minus_upsilon11,
                  y21,  inc_y,
                  a12t, cs_A );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  &minus_zeta11,
                  v21,  inc_v,
                  a12t, cs_A );
    }

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    // FLA_Copy( a21, u21p );
    FLA_Househ2_UT_l_opz( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );
    bl1_zcopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                a21,  rs_A,
                u21p, inc_up );

    if ( n_ahead > 0 )
    {
      // FLA_Copy( FLA_MINUS_ONE, minus_inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, minus_inv_tau11 );
      bl1_zdiv3( buff_m1, tau11, &minus_inv_tau11 );

      // FLA_Copyt( FLA_TRANSPOSE, a12t, a12p );
      // FLA_Axpyt( FLA_TRANSPOSE, minus_inv_tau11, a12t, a12p );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  a12p, inc_ap );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  &minus_inv_tau11,
                  a12t, cs_A,
                  a12p, inc_ap );
    }

    if ( m_behind > 0 && n_ahead > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u21, y21, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );
      // FLA_Gemvc( FLA_CONJ_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, u21p, FLA_ZERO, y21 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      FLA_Fused_Gerc2_Ahx_Axpy_Ax_opz_var1( m_ahead,
                                            n_ahead,
                                            tau11,
                                            buff_m1,
                                            u21, inc_u,
                                            y21, inc_y,
                                            z21, inc_z,
                                            v21, inc_v,
                                            A22, rs_A, cs_A,
                                            u21p, inc_up,
                                            a12p, inc_ap,
                                            w21,  inc_w );
                                           
                                           
    }
    else if ( n_ahead > 0 )
    {
      // FLA_Gemvc( FLA_CONJ_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, u21p, FLA_ZERO, y21 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      FLA_Fused_Ahx_Axpy_Ax_opz_var1( m_ahead,
                                      n_ahead,
                                      tau11,
                                      buff_0,
                                      A22,  rs_A, cs_A,
                                      u21p, inc_up,
                                      a12p, inc_ap,
                                      y21,  inc_y,
                                      w21,  inc_w );
    }

    if ( n_ahead > 0 )
    {
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_ONE, a12t, y21 );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  buff_1,
                  a12t, cs_A,
                  y21,  inc_y );

      // FLA_Househ2s_UT( FLA_RIGHT,
      //                  a12p_t,
      //                  a12p_b,
      //                  alpha12, psi11_minus_alpha12, sigma11 );
      FLA_Househ2s_UT_r_opz( n_ahead - 1,
                             a12p_t,
                             a12p_b, inc_ap,
                             &alpha12,
                             &psi11_minus_alpha12,
                             sigma11 );

      // FLA_Copy( a12p, v21 );
      // FLA_Mult_add( FLA_MINUS_ONE, alpha12, v21_t );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, psi11_minus_alpha12, v21 );
      // FLA_Conjugate( v21_b );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12p, inc_ap,
                  v21,  inc_v );
      bl1_zmult4( buff_m1, &alpha12, v21_t, v21_t );
      bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                     n_ahead,
                     &psi11_minus_alpha12,
                     v21, inc_v );
      bl1_zconjv( n_ahead - 1,
                  v21_b, inc_v );

      // FLA_Copyt( FLA_NO_TRANSPOSE, alpha12, a12t_l );
      // FLA_Copyt( FLA_TRANSPOSE, v21_b, a12t_r );
      *a12t_l = alpha12;
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  n_ahead - 1,
                  v21_b,  inc_v,
                  a12t_r, cs_A );
    }

    // FLA_Copy( u21p, u21 );
    bl1_zcopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                u21p, inc_up,
                u21,  inc_u );

    if ( n_ahead > 0 )
    {
      // FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, beta );
      bl1_zdot( BLIS1_CONJUGATE,
                n_ahead,
                y21, inc_y,
                v21, inc_v,
                &beta );
      bl1_zscals( &minus_inv_tau11, &beta );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha12, minus_conj_alpha12 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_alpha12 );
      bl1_zcopyconj( &alpha12, &minus_conj_alpha12 );
      bl1_zneg1( &minus_conj_alpha12 );

      // FLA_Copy( w21, z21 );
      // FLA_Axpy( minus_conj_alpha12, A22_l, z21 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, z21 );
      // FLA_Axpy( beta, u21, z21 );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, inc_z );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_alpha12,
                  A22_l, rs_A,
                  z21,   inc_z );
      bl1_zinvscalv( BLIS1_CONJUGATE,
                     m_ahead,
                     &psi11_minus_alpha12,
                     z21, inc_z );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  u21, inc_u,
                  z21, inc_z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11,   y21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_zinvscalv( BLIS1_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, inc_y );
      bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, inc_z );

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

    // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ONE, t01 );
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
               u21, inc_u,
               buff_1,
               t01, rs_T );

    if ( m_behind + 1 == b_alg && n_ahead > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u21, y21, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );
      bl1_zger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u21, inc_u,
                y21, inc_y,
                A22, rs_A, cs_A );
      bl1_zger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z21, inc_z,
                v21, inc_v,
                A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  // FLA_Obj_free( &ap );
  // FLA_Obj_free( &u );
  // FLA_Obj_free( &up );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_w );
  FLA_free( buff_ap );
  FLA_free( buff_u );
  FLA_free( buff_up );
  FLA_free( buff_v );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}

