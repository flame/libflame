/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_opt_var4( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  FLA_Error    r_val;
  FLA_Obj      Y, Z;
  FLA_Datatype datatype_A;
  dim_t        m_A, n_A;

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );
  n_A        = FLA_Obj_width( A );
  
  FLA_Obj_create( datatype_A, n_A, n_A, 0, 0, &Y );
  FLA_Obj_create( datatype_A, m_A, n_A, 0, 0, &Z );

  r_val = FLA_Bidiag_UT_u_step_opt_var4( A, Y, Z, TU, TV );

  FLA_Obj_free( &Y );
  FLA_Obj_free( &Z );

  return r_val;
}

FLA_Error FLA_Bidiag_UT_u_step_opt_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj S )
{
  FLA_Datatype datatype;
  integer          m_A, n_A, m_TS;
  integer          rs_A, cs_A;
  integer          rs_Y, cs_Y;
  integer          rs_Z, cs_Z;
  integer          rs_T, cs_T;
  integer          rs_S, cs_S;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  m_TS     = FLA_Obj_length( T );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_Y     = FLA_Obj_row_stride( Y );
  cs_Y     = FLA_Obj_col_stride( Y );

  rs_Z     = FLA_Obj_row_stride( Z );
  cs_Z     = FLA_Obj_col_stride( Z );

  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );
  
  rs_S     = FLA_Obj_row_stride( S );
  cs_S     = FLA_Obj_col_stride( S );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_Y = FLA_FLOAT_PTR( Y );
      float* buff_Z = FLA_FLOAT_PTR( Z );
      float* buff_T = FLA_FLOAT_PTR( T );
      float* buff_S = FLA_FLOAT_PTR( S );

      FLA_Bidiag_UT_u_step_ops_var4( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_Y, rs_Y, cs_Y,
                                     buff_Z, rs_Z, cs_Z,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_Y = FLA_DOUBLE_PTR( Y );
      double* buff_Z = FLA_DOUBLE_PTR( Z );
      double* buff_T = FLA_DOUBLE_PTR( T );
      double* buff_S = FLA_DOUBLE_PTR( S );

      FLA_Bidiag_UT_u_step_opd_var4( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_Y, rs_Y, cs_Y,
                                     buff_Z, rs_Z, cs_Z,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_Y = FLA_COMPLEX_PTR( Y );
      scomplex* buff_Z = FLA_COMPLEX_PTR( Z );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );
      scomplex* buff_S = FLA_COMPLEX_PTR( S );

      FLA_Bidiag_UT_u_step_opc_var4( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_Y, rs_Y, cs_Y,
                                     buff_Z, rs_Z, cs_Z,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_Y = FLA_DOUBLE_COMPLEX_PTR( Y );
      dcomplex* buff_Z = FLA_DOUBLE_COMPLEX_PTR( Z );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );
      dcomplex* buff_S = FLA_DOUBLE_COMPLEX_PTR( S );

      FLA_Bidiag_UT_u_step_opz_var4( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_Y, rs_Y, cs_Y,
                                     buff_Z, rs_Z, cs_Z,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_ops_var4( integer m_A,
                                         integer n_A,
                                         integer m_TS,
                                         float* buff_A, integer rs_A, integer cs_A, 
                                         float* buff_Y, integer rs_Y, integer cs_Y, 
                                         float* buff_Z, integer rs_Z, integer cs_Z, 
                                         float* buff_T, integer rs_T, integer cs_T, 
                                         float* buff_S, integer rs_S, integer cs_S )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );

  float     alpha12;
  float     minus_conj_alpha12;
  float     psi11_minus_alpha12;
  float     minus_inv_tau11;
  float     beta;
  float     last_elem;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &al );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &ap );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &up );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &g );
  float*    buff_w  = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_al = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_ap = ( float* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_u  = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_up = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_v  = ( float* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_d  = ( float* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_e  = ( float* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_f  = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_g  = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_w   = 1;
  integer       inc_al  = 1;
  integer       inc_ap  = 1;
  integer       inc_u   = 1;
  integer       inc_up  = 1;
  integer       inc_v   = 1;
  integer       inc_d   = 1;
  integer       inc_e   = 1;
  integer       inc_f   = 1;
  integer       inc_g   = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bl1_ssetm( n_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bl1_ssetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    float*    a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    float*    A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    float*    a01       = buff_A  + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    float*    a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    float*    A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    float*    A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    float*    y10t      = buff_Y  + (0  )*cs_Y + (i  )*rs_Y;
    float*    Y20       = buff_Y  + (0  )*cs_Y + (i+1)*rs_Y;
    float*    y21       = buff_Y  + (i  )*cs_Y + (i+1)*rs_Y;

    float*    z10t      = buff_Z  + (0  )*cs_Z + (i  )*rs_Z;
    float*    Z20       = buff_Z  + (0  )*cs_Z + (i+1)*rs_Z;
    float*    z21       = buff_Z  + (i  )*cs_Z + (i+1)*rs_Z;

    float*    t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    float*    tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    float*    s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    float*    sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    float*    w21       = buff_w  + (i+1)*inc_w;

    float*    a22l      = buff_al + (i+1)*inc_al;

    float*    a12p      = buff_ap + (i+1)*inc_ap;

    float*    u21       = buff_u  + (i+1)*inc_u;

    float*    u21p      = buff_up + (i+1)*inc_up;

    float*    v21       = buff_v  + (i+1)*inc_v;

    float*    d0        = buff_d  + (0  )*inc_d;

    float*    e0        = buff_e  + (0  )*inc_e;

    float*    f0        = buff_f  + (0  )*inc_f;

    float*    g0        = buff_g  + (0  )*inc_g;

    float*    a12p_t    = a12p    + (0  )*inc_ap;
    float*    a12p_b    = a12p    + (1  )*inc_ap;

    float*    v21_t     = v21     + (0  )*inc_v;
    float*    v21_b     = v21     + (1  )*inc_v;

    float*    a01_b     = a01     + (0  )*cs_A + (i-1)*rs_A;

    float*    a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    float*    a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    float*    A02_l     = A02     + (0  )*cs_A + (0  )*rs_A;

    float*    A22_l     = A22     + (0  )*cs_A + (0  )*rs_A;

    float*    Y20_t     = Y20     + (0  )*cs_Y + (0  )*rs_Y;

    float*    ABL       = a10t;
    float*    ZBL       = z10t;

    float*    a2        = alpha11;

    integer       m_ahead   = m_A - i - 1;
    integer       n_ahead   = n_A - i - 1;
    integer       m_behind  = i;
    integer       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a01_b, last_elem );
      // FLA_Set( FLA_ONE, a01_b );
      last_elem = *a01_b;
      *a01_b = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a01,  FLA_ONE, a2 );
    bl1_sgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               y10t, cs_Y,
               buff_1,
               a2,   rs_A );
    bl1_sgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL, rs_Z, cs_Z,
               a01, rs_A,
               buff_1,
               a2,  rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE,    FLA_MINUS_ONE, A02, z10t, FLA_ONE, a12t );
    bl1_sgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               n_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a01_b );
      *a01_b = last_elem;
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

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21p, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, u21p, FLA_ZERO, e0 );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20,  rs_A, cs_A,
                 u21p, inc_up,
                 buff_0,
                 d0,   inc_d );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20,  rs_Z, cs_Z,
                 u21p, inc_up,
                 buff_0,
                 e0,   inc_e );

      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Axpy( FLA_ONE, d0, t01 );
      bl1_scopyv( BLIS1_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  buff_1,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Set( FLA_ZERO, y21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_TRANSPOSE,    FLA_MINUS_ONE, A02, e0, FLA_ONE, y21 );
      bl1_ssetv( n_ahead,
                 buff_0,
                 y21, rs_Y );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bl1_sgemv( BLIS1_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_m1,
                 A02, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, u21p, FLA_ONE, y21 );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22,  rs_A, cs_A,
                 u21p, inc_up,
                 buff_1,
                 y21,  rs_Y );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      bl1_saxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  &minus_inv_tau11,
                  y21,  rs_Y,
                  a12p, inc_ap );

      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22,  rs_A, cs_A,
                 a12p, inc_ap,
                 buff_0,
                 w21,  inc_w );

      // FLA_Gemvc( FLA_CONJ_TRANSPOSE,    FLA_CONJUGATE, FLA_ONE, Y20, a12p, FLA_ZERO, f0 );
      // FLA_Gemvc( FLA_CONJ_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A02, a12p, FLA_ZERO, g0 );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_1,
                 Y20,  rs_Y, cs_Y,
                 a12p, inc_ap,
                 buff_0,
                 f0,   inc_f );
      bl1_sgemv( BLIS1_CONJ_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02,  rs_A, cs_A,
                 a12p, inc_ap,
                 buff_0,
                 g0,   inc_g );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, w21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, g0, FLA_ONE, w21 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 w21, inc_w );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 g0,  inc_g,
                 buff_1,
                 w21, inc_w );

      // FLA_Copy( A22_l, a22l );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A20, Y20_t, FLA_ONE, a22l );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, Z20, A02_l, FLA_ONE, a22l );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  A22_l, rs_A,
                  a22l,  inc_al );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20,   rs_A, cs_A,
                 Y20_t, cs_Y,
                 buff_1,
                 a22l,  inc_al );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20,   rs_Z, cs_Z,
                 A02_l, rs_A,
                 buff_1,
                 a22l,  inc_al );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_ONE, a12t, y21 );
      bl1_saxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  buff_1, 
                  a12t, cs_A,
                  y21,  rs_Y );

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
      // FLA_Conjugate( v21_b );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12p, inc_ap,
                  v21,  inc_v );
      bl1_smult4( buff_m1, &alpha12, v21_t, v21_t );
      bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                     n_ahead,
                     &psi11_minus_alpha12,
                     v21, inc_v );
      bl1_sconjv( n_ahead - 1,
                  v21_b, inc_v );

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha12, minus_conj_alpha12 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_alpha12 );
      bl1_scopyconj( &alpha12, &minus_conj_alpha12 );
      bl1_sneg1( &minus_conj_alpha12 );

      // FLA_Copy( g0, s01 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_conj_alpha12, A02_l, s01 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, s01 );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  g0,  inc_g,
                  s01, rs_S );
      bl1_saxpyv( BLIS1_CONJUGATE,
                  n_behind,
                  &minus_conj_alpha12,
                  A02_l, rs_A,
                  s01,   rs_S );
      bl1_sinvscalv( BLIS1_CONJUGATE,
                     n_behind,
                     &psi11_minus_alpha12,
                     s01, rs_S );

      // FLA_Copyt( FLA_NO_TRANSPOSE, alpha12, a12t_l );
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
                y21, rs_Y,
                v21, inc_v,
                &beta );
      bl1_sscals( &minus_inv_tau11, &beta );

      // FLA_Copy( w21, z21 );
      // FLA_Axpy( minus_conj_alpha12, a22l, z21 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, z21 );
      // FLA_Axpy( beta, u21, z21 );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, rs_Z );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_alpha12,
                  a22l, inc_al,
                  z21,  rs_Z );
      bl1_sinvscalv( BLIS1_CONJUGATE,
                     m_ahead,
                     &psi11_minus_alpha12,
                     z21, rs_Z );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  u21, inc_u,
                  z21, rs_Z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11,   y21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_sinvscalv( BLIS1_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, rs_Y );
      bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, rs_Z );
    }
	else // if ( n_ahead == 0 )
    {
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
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  // FLA_Obj_free( &al );
  // FLA_Obj_free( &ap );
  // FLA_Obj_free( &u );
  // FLA_Obj_free( &up );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  // FLA_Obj_free( &g );
  FLA_free( buff_w );
  FLA_free( buff_al );
  FLA_free( buff_ap );
  FLA_free( buff_u );
  FLA_free( buff_up );
  FLA_free( buff_v );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );
  FLA_free( buff_g );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_opd_var4( integer m_A,
                                         integer n_A,
                                         integer m_TS,
                                         double* buff_A, integer rs_A, integer cs_A, 
                                         double* buff_Y, integer rs_Y, integer cs_Y, 
                                         double* buff_Z, integer rs_Z, integer cs_Z, 
                                         double* buff_T, integer rs_T, integer cs_T, 
                                         double* buff_S, integer rs_S, integer cs_S )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );

  double    alpha12;
  double    minus_conj_alpha12;
  double    psi11_minus_alpha12;
  double    minus_inv_tau11;
  double    beta;
  double    last_elem;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &al );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &ap );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &up );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &g );
  double*   buff_w  = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_al = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_ap = ( double* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_u  = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_up = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_v  = ( double* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_d  = ( double* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_e  = ( double* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_f  = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_g  = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_w   = 1;
  integer       inc_al  = 1;
  integer       inc_ap  = 1;
  integer       inc_u   = 1;
  integer       inc_up  = 1;
  integer       inc_v   = 1;
  integer       inc_d   = 1;
  integer       inc_e   = 1;
  integer       inc_f   = 1;
  integer       inc_g   = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bl1_dsetm( n_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bl1_dsetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    double*   a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    double*   A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    double*   a01       = buff_A  + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    double*   a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    double*   A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    double*   A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    double*   y10t      = buff_Y  + (0  )*cs_Y + (i  )*rs_Y;
    double*   Y20       = buff_Y  + (0  )*cs_Y + (i+1)*rs_Y;
    double*   y21       = buff_Y  + (i  )*cs_Y + (i+1)*rs_Y;

    double*   z10t      = buff_Z  + (0  )*cs_Z + (i  )*rs_Z;
    double*   Z20       = buff_Z  + (0  )*cs_Z + (i+1)*rs_Z;
    double*   z21       = buff_Z  + (i  )*cs_Z + (i+1)*rs_Z;

    double*   t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    double*   tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    double*   s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    double*   sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    double*   w21       = buff_w  + (i+1)*inc_w;

    double*   a22l      = buff_al + (i+1)*inc_al;

    double*   a12p      = buff_ap + (i+1)*inc_ap;

    double*   u21       = buff_u  + (i+1)*inc_u;

    double*   u21p      = buff_up + (i+1)*inc_up;

    double*   v21       = buff_v  + (i+1)*inc_v;

    double*   d0        = buff_d  + (0  )*inc_d;

    double*   e0        = buff_e  + (0  )*inc_e;

    double*   f0        = buff_f  + (0  )*inc_f;

    double*   g0        = buff_g  + (0  )*inc_g;

    double*   a12p_t    = a12p    + (0  )*inc_ap;
    double*   a12p_b    = a12p    + (1  )*inc_ap;

    double*   v21_t     = v21     + (0  )*inc_v;
    double*   v21_b     = v21     + (1  )*inc_v;

    double*   a01_b     = a01     + (0  )*cs_A + (i-1)*rs_A;

    double*   a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    double*   a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    double*   A02_l     = A02     + (0  )*cs_A + (0  )*rs_A;

    double*   A22_l     = A22     + (0  )*cs_A + (0  )*rs_A;

    double*   Y20_t     = Y20     + (0  )*cs_Y + (0  )*rs_Y;

    double*   ABL       = a10t;
    double*   ZBL       = z10t;

    double*   a2        = alpha11;

    integer       m_ahead   = m_A - i - 1;
    integer       n_ahead   = n_A - i - 1;
    integer       m_behind  = i;
    integer       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a01_b, last_elem );
      // FLA_Set( FLA_ONE, a01_b );
      last_elem = *a01_b;
      *a01_b = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a01,  FLA_ONE, a2 );
    bl1_dgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               y10t, cs_Y,
               buff_1,
               a2,   rs_A );
    bl1_dgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL, rs_Z, cs_Z,
               a01, rs_A,
               buff_1,
               a2,  rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE,    FLA_MINUS_ONE, A02, z10t, FLA_ONE, a12t );
    bl1_dgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               n_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a01_b );
      *a01_b = last_elem;
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

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21p, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, u21p, FLA_ZERO, e0 );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20,  rs_A, cs_A,
                 u21p, inc_up,
                 buff_0,
                 d0,   inc_d );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20,  rs_Z, cs_Z,
                 u21p, inc_up,
                 buff_0,
                 e0,   inc_e );

      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Axpy( FLA_ONE, d0, t01 );
      bl1_dcopyv( BLIS1_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  buff_1,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Set( FLA_ZERO, y21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_TRANSPOSE,    FLA_MINUS_ONE, A02, e0, FLA_ONE, y21 );
      bl1_dsetv( n_ahead,
                 buff_0,
                 y21, rs_Y );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bl1_dgemv( BLIS1_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_m1,
                 A02, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, u21p, FLA_ONE, y21 );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22,  rs_A, cs_A,
                 u21p, inc_up,
                 buff_1,
                 y21,  rs_Y );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      bl1_daxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  &minus_inv_tau11,
                  y21,  rs_Y,
                  a12p, inc_ap );

      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22,  rs_A, cs_A,
                 a12p, inc_ap,
                 buff_0,
                 w21,  inc_w );

      // FLA_Gemvc( FLA_CONJ_TRANSPOSE,    FLA_CONJUGATE, FLA_ONE, Y20, a12p, FLA_ZERO, f0 );
      // FLA_Gemvc( FLA_CONJ_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A02, a12p, FLA_ZERO, g0 );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_1,
                 Y20,  rs_Y, cs_Y,
                 a12p, inc_ap,
                 buff_0,
                 f0,   inc_f );
      bl1_dgemv( BLIS1_CONJ_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02,  rs_A, cs_A,
                 a12p, inc_ap,
                 buff_0,
                 g0,   inc_g );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, w21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, g0, FLA_ONE, w21 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 w21, inc_w );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 g0,  inc_g,
                 buff_1,
                 w21, inc_w );

      // FLA_Copy( A22_l, a22l );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A20, Y20_t, FLA_ONE, a22l );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, Z20, A02_l, FLA_ONE, a22l );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  A22_l, rs_A,
                  a22l,  inc_al );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20,   rs_A, cs_A,
                 Y20_t, cs_Y,
                 buff_1,
                 a22l,  inc_al );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20,   rs_Z, cs_Z,
                 A02_l, rs_A,
                 buff_1,
                 a22l,  inc_al );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_ONE, a12t, y21 );
      bl1_daxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  buff_1, 
                  a12t, cs_A,
                  y21,  rs_Y );

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

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha12, minus_conj_alpha12 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_alpha12 );
      bl1_dcopyconj( &alpha12, &minus_conj_alpha12 );
      bl1_dneg1( &minus_conj_alpha12 );

      // FLA_Copy( g0, s01 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_conj_alpha12, A02_l, s01 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, s01 );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  g0,  inc_g,
                  s01, rs_S );
      bl1_daxpyv( BLIS1_CONJUGATE,
                  n_behind,
                  &minus_conj_alpha12,
                  A02_l, rs_A,
                  s01,   rs_S );
      bl1_dinvscalv( BLIS1_CONJUGATE,
                     n_behind,
                     &psi11_minus_alpha12,
                     s01, rs_S );

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
                y21, rs_Y,
                v21, inc_v,
                &beta );
      bl1_dscals( &minus_inv_tau11, &beta );

      // FLA_Copy( w21, z21 );
      // FLA_Axpy( minus_conj_alpha12, a22l, z21 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, z21 );
      // FLA_Axpy( beta, u21, z21 );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, rs_Z );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_alpha12,
                  a22l, inc_al,
                  z21,  rs_Z );
      bl1_dinvscalv( BLIS1_CONJUGATE,
                     m_ahead,
                     &psi11_minus_alpha12,
                     z21, rs_Z );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  u21, inc_u,
                  z21, rs_Z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11,   y21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_dinvscalv( BLIS1_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, rs_Y );
      bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, rs_Z );
    }
	else // if ( n_ahead == 0 )
    {
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
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  // FLA_Obj_free( &al );
  // FLA_Obj_free( &ap );
  // FLA_Obj_free( &u );
  // FLA_Obj_free( &up );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  // FLA_Obj_free( &g );
  FLA_free( buff_w );
  FLA_free( buff_al );
  FLA_free( buff_ap );
  FLA_free( buff_u );
  FLA_free( buff_up );
  FLA_free( buff_v );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );
  FLA_free( buff_g );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_opc_var4( integer m_A,
                                         integer n_A,
                                         integer m_TS,
                                         scomplex* buff_A, integer rs_A, integer cs_A, 
                                         scomplex* buff_Y, integer rs_Y, integer cs_Y, 
                                         scomplex* buff_Z, integer rs_Z, integer cs_Z, 
                                         scomplex* buff_T, integer rs_T, integer cs_T, 
                                         scomplex* buff_S, integer rs_S, integer cs_S )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );

  scomplex  alpha12;
  scomplex  minus_conj_alpha12;
  scomplex  psi11_minus_alpha12;
  scomplex  minus_inv_tau11;
  scomplex  beta;
  scomplex  last_elem;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &al );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &ap );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &up );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &g );
  scomplex* buff_w  = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_al = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_ap = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_u  = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_up = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_v  = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_d  = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_e  = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_f  = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_g  = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_w   = 1;
  integer       inc_al  = 1;
  integer       inc_ap  = 1;
  integer       inc_u   = 1;
  integer       inc_up  = 1;
  integer       inc_v   = 1;
  integer       inc_d   = 1;
  integer       inc_e   = 1;
  integer       inc_f   = 1;
  integer       inc_g   = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bl1_csetm( n_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bl1_csetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    scomplex* a01       = buff_A  + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* y10t      = buff_Y  + (0  )*cs_Y + (i  )*rs_Y;
    scomplex* Y20       = buff_Y  + (0  )*cs_Y + (i+1)*rs_Y;
    scomplex* y21       = buff_Y  + (i  )*cs_Y + (i+1)*rs_Y;

    scomplex* z10t      = buff_Z  + (0  )*cs_Z + (i  )*rs_Z;
    scomplex* Z20       = buff_Z  + (0  )*cs_Z + (i+1)*rs_Z;
    scomplex* z21       = buff_Z  + (i  )*cs_Z + (i+1)*rs_Z;

    scomplex* t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    scomplex* s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    scomplex* sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    scomplex* w21       = buff_w  + (i+1)*inc_w;

    scomplex* a22l      = buff_al + (i+1)*inc_al;

    scomplex* a12p      = buff_ap + (i+1)*inc_ap;

    scomplex* u21       = buff_u  + (i+1)*inc_u;

    scomplex* u21p      = buff_up + (i+1)*inc_up;

    scomplex* v21       = buff_v  + (i+1)*inc_v;

    scomplex* d0        = buff_d  + (0  )*inc_d;

    scomplex* e0        = buff_e  + (0  )*inc_e;

    scomplex* f0        = buff_f  + (0  )*inc_f;

    scomplex* g0        = buff_g  + (0  )*inc_g;

    scomplex* a12p_t    = a12p    + (0  )*inc_ap;
    scomplex* a12p_b    = a12p    + (1  )*inc_ap;

    scomplex* v21_t     = v21     + (0  )*inc_v;
    scomplex* v21_b     = v21     + (1  )*inc_v;

    scomplex* a01_b     = a01     + (0  )*cs_A + (i-1)*rs_A;

    scomplex* a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    scomplex* a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    scomplex* A02_l     = A02     + (0  )*cs_A + (0  )*rs_A;

    scomplex* A22_l     = A22     + (0  )*cs_A + (0  )*rs_A;

    scomplex* Y20_t     = Y20     + (0  )*cs_Y + (0  )*rs_Y;

    scomplex* ABL       = a10t;
    scomplex* ZBL       = z10t;

    scomplex* a2        = alpha11;

    integer       m_ahead   = m_A - i - 1;
    integer       n_ahead   = n_A - i - 1;
    integer       m_behind  = i;
    integer       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a01_b, last_elem );
      // FLA_Set( FLA_ONE, a01_b );
      last_elem = *a01_b;
      *a01_b = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a01,  FLA_ONE, a2 );
    bl1_cgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               y10t, cs_Y,
               buff_1,
               a2,   rs_A );
    bl1_cgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL, rs_Z, cs_Z,
               a01, rs_A,
               buff_1,
               a2,  rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE,    FLA_MINUS_ONE, A02, z10t, FLA_ONE, a12t );
    bl1_cgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               n_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a01_b );
      *a01_b = last_elem;
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

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21p, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, u21p, FLA_ZERO, e0 );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20,  rs_A, cs_A,
                 u21p, inc_up,
                 buff_0,
                 d0,   inc_d );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20,  rs_Z, cs_Z,
                 u21p, inc_up,
                 buff_0,
                 e0,   inc_e );

      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Axpy( FLA_ONE, d0, t01 );
      bl1_ccopyv( BLIS1_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  buff_1,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Set( FLA_ZERO, y21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_TRANSPOSE,    FLA_MINUS_ONE, A02, e0, FLA_ONE, y21 );
      bl1_csetv( n_ahead,
                 buff_0,
                 y21, rs_Y );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bl1_cgemv( BLIS1_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_m1,
                 A02, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, u21p, FLA_ONE, y21 );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22,  rs_A, cs_A,
                 u21p, inc_up,
                 buff_1,
                 y21,  rs_Y );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  &minus_inv_tau11,
                  y21,  rs_Y,
                  a12p, inc_ap );

      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22,  rs_A, cs_A,
                 a12p, inc_ap,
                 buff_0,
                 w21,  inc_w );

      // FLA_Gemvc( FLA_CONJ_TRANSPOSE,    FLA_CONJUGATE, FLA_ONE, Y20, a12p, FLA_ZERO, f0 );
      // FLA_Gemvc( FLA_CONJ_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A02, a12p, FLA_ZERO, g0 );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_1,
                 Y20,  rs_Y, cs_Y,
                 a12p, inc_ap,
                 buff_0,
                 f0,   inc_f );
      bl1_cgemv( BLIS1_CONJ_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02,  rs_A, cs_A,
                 a12p, inc_ap,
                 buff_0,
                 g0,   inc_g );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, w21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, g0, FLA_ONE, w21 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 w21, inc_w );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 g0,  inc_g,
                 buff_1,
                 w21, inc_w );

      // FLA_Copy( A22_l, a22l );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A20, Y20_t, FLA_ONE, a22l );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, Z20, A02_l, FLA_ONE, a22l );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  A22_l, rs_A,
                  a22l,  inc_al );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20,   rs_A, cs_A,
                 Y20_t, cs_Y,
                 buff_1,
                 a22l,  inc_al );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20,   rs_Z, cs_Z,
                 A02_l, rs_A,
                 buff_1,
                 a22l,  inc_al );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_ONE, a12t, y21 );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  buff_1, 
                  a12t, cs_A,
                  y21,  rs_Y );

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

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha12, minus_conj_alpha12 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_alpha12 );
      bl1_ccopyconj( &alpha12, &minus_conj_alpha12 );
      bl1_cneg1( &minus_conj_alpha12 );

      // FLA_Copy( g0, s01 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_conj_alpha12, A02_l, s01 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, s01 );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  g0,  inc_g,
                  s01, rs_S );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  n_behind,
                  &minus_conj_alpha12,
                  A02_l, rs_A,
                  s01,   rs_S );
      bl1_cinvscalv( BLIS1_CONJUGATE,
                     n_behind,
                     &psi11_minus_alpha12,
                     s01, rs_S );

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
                y21, rs_Y,
                v21, inc_v,
                &beta );
      bl1_cscals( &minus_inv_tau11, &beta );

      // FLA_Copy( w21, z21 );
      // FLA_Axpy( minus_conj_alpha12, a22l, z21 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, z21 );
      // FLA_Axpy( beta, u21, z21 );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, rs_Z );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_alpha12,
                  a22l, inc_al,
                  z21,  rs_Z );
      bl1_cinvscalv( BLIS1_CONJUGATE,
                     m_ahead,
                     &psi11_minus_alpha12,
                     z21, rs_Z );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  u21, inc_u,
                  z21, rs_Z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11,   y21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_cinvscalv( BLIS1_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, rs_Y );
      bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, rs_Z );
    }
	else // if ( n_ahead == 0 )
    {
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
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  // FLA_Obj_free( &al );
  // FLA_Obj_free( &ap );
  // FLA_Obj_free( &u );
  // FLA_Obj_free( &up );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  // FLA_Obj_free( &g );
  FLA_free( buff_w );
  FLA_free( buff_al );
  FLA_free( buff_ap );
  FLA_free( buff_u );
  FLA_free( buff_up );
  FLA_free( buff_v );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );
  FLA_free( buff_g );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_opz_var4( integer m_A,
                                         integer n_A,
                                         integer m_TS,
                                         dcomplex* buff_A, integer rs_A, integer cs_A, 
                                         dcomplex* buff_Y, integer rs_Y, integer cs_Y, 
                                         dcomplex* buff_Z, integer rs_Z, integer cs_Z, 
                                         dcomplex* buff_T, integer rs_T, integer cs_T, 
                                         dcomplex* buff_S, integer rs_S, integer cs_S )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );

  dcomplex  alpha12;
  dcomplex  minus_conj_alpha12;
  dcomplex  psi11_minus_alpha12;
  dcomplex  minus_inv_tau11;
  dcomplex  beta;
  dcomplex  last_elem;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &al );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &ap );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &up );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &g );
  dcomplex* buff_w  = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_al = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_ap = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_u  = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_up = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_v  = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_d  = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_e  = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_f  = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_g  = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_w   = 1;
  integer       inc_al  = 1;
  integer       inc_ap  = 1;
  integer       inc_u   = 1;
  integer       inc_up  = 1;
  integer       inc_v   = 1;
  integer       inc_d   = 1;
  integer       inc_e   = 1;
  integer       inc_f   = 1;
  integer       inc_g   = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bl1_zsetm( n_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bl1_zsetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* a01       = buff_A  + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* y10t      = buff_Y  + (0  )*cs_Y + (i  )*rs_Y;
    dcomplex* Y20       = buff_Y  + (0  )*cs_Y + (i+1)*rs_Y;
    dcomplex* y21       = buff_Y  + (i  )*cs_Y + (i+1)*rs_Y;

    dcomplex* z10t      = buff_Z  + (0  )*cs_Z + (i  )*rs_Z;
    dcomplex* Z20       = buff_Z  + (0  )*cs_Z + (i+1)*rs_Z;
    dcomplex* z21       = buff_Z  + (i  )*cs_Z + (i+1)*rs_Z;

    dcomplex* t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    dcomplex* s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    dcomplex* sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    dcomplex* w21       = buff_w  + (i+1)*inc_w;

    dcomplex* a22l      = buff_al + (i+1)*inc_al;

    dcomplex* a12p      = buff_ap + (i+1)*inc_ap;

    dcomplex* u21       = buff_u  + (i+1)*inc_u;

    dcomplex* u21p      = buff_up + (i+1)*inc_up;

    dcomplex* v21       = buff_v  + (i+1)*inc_v;

    dcomplex* d0        = buff_d  + (0  )*inc_d;

    dcomplex* e0        = buff_e  + (0  )*inc_e;

    dcomplex* f0        = buff_f  + (0  )*inc_f;

    dcomplex* g0        = buff_g  + (0  )*inc_g;

    dcomplex* a12p_t    = a12p    + (0  )*inc_ap;
    dcomplex* a12p_b    = a12p    + (1  )*inc_ap;

    dcomplex* v21_t     = v21     + (0  )*inc_v;
    dcomplex* v21_b     = v21     + (1  )*inc_v;

    dcomplex* a01_b     = a01     + (0  )*cs_A + (i-1)*rs_A;

    dcomplex* a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    dcomplex* A02_l     = A02     + (0  )*cs_A + (0  )*rs_A;

    dcomplex* A22_l     = A22     + (0  )*cs_A + (0  )*rs_A;

    dcomplex* Y20_t     = Y20     + (0  )*cs_Y + (0  )*rs_Y;

    dcomplex* ABL       = a10t;
    dcomplex* ZBL       = z10t;

    dcomplex* a2        = alpha11;

    integer       m_ahead   = m_A - i - 1;
    integer       n_ahead   = n_A - i - 1;
    integer       m_behind  = i;
    integer       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a01_b, last_elem );
      // FLA_Set( FLA_ONE, a01_b );
      last_elem = *a01_b;
      *a01_b = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a01,  FLA_ONE, a2 );
    bl1_zgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               y10t, cs_Y,
               buff_1,
               a2,   rs_A );
    bl1_zgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL, rs_Z, cs_Z,
               a01, rs_A,
               buff_1,
               a2,  rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE,    FLA_MINUS_ONE, A02, z10t, FLA_ONE, a12t );
    bl1_zgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               n_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a01_b );
      *a01_b = last_elem;
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

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21p, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, u21p, FLA_ZERO, e0 );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20,  rs_A, cs_A,
                 u21p, inc_up,
                 buff_0,
                 d0,   inc_d );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20,  rs_Z, cs_Z,
                 u21p, inc_up,
                 buff_0,
                 e0,   inc_e );

      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Axpy( FLA_ONE, d0, t01 );
      bl1_zcopyv( BLIS1_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  buff_1,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Set( FLA_ZERO, y21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_TRANSPOSE,    FLA_MINUS_ONE, A02, e0, FLA_ONE, y21 );
      bl1_zsetv( n_ahead,
                 buff_0,
                 y21, rs_Y );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bl1_zgemv( BLIS1_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_m1,
                 A02, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, u21p, FLA_ONE, y21 );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22,  rs_A, cs_A,
                 u21p, inc_up,
                 buff_1,
                 y21,  rs_Y );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  &minus_inv_tau11,
                  y21,  rs_Y,
                  a12p, inc_ap );

      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22,  rs_A, cs_A,
                 a12p, inc_ap,
                 buff_0,
                 w21,  inc_w );

      // FLA_Gemvc( FLA_CONJ_TRANSPOSE,    FLA_CONJUGATE, FLA_ONE, Y20, a12p, FLA_ZERO, f0 );
      // FLA_Gemvc( FLA_CONJ_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A02, a12p, FLA_ZERO, g0 );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_1,
                 Y20,  rs_Y, cs_Y,
                 a12p, inc_ap,
                 buff_0,
                 f0,   inc_f );
      bl1_zgemv( BLIS1_CONJ_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02,  rs_A, cs_A,
                 a12p, inc_ap,
                 buff_0,
                 g0,   inc_g );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, w21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, g0, FLA_ONE, w21 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 w21, inc_w );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 g0,  inc_g,
                 buff_1,
                 w21, inc_w );

      // FLA_Copy( A22_l, a22l );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A20, Y20_t, FLA_ONE, a22l );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, Z20, A02_l, FLA_ONE, a22l );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  A22_l, rs_A,
                  a22l,  inc_al );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20,   rs_A, cs_A,
                 Y20_t, cs_Y,
                 buff_1,
                 a22l,  inc_al );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20,   rs_Z, cs_Z,
                 A02_l, rs_A,
                 buff_1,
                 a22l,  inc_al );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_ONE, a12t, y21 );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  buff_1, 
                  a12t, cs_A,
                  y21,  rs_Y );

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

      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha12, minus_conj_alpha12 );
      // FLA_Scal( FLA_MINUS_ONE, minus_conj_alpha12 );
      bl1_zcopyconj( &alpha12, &minus_conj_alpha12 );
      bl1_zneg1( &minus_conj_alpha12 );

      // FLA_Copy( g0, s01 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_conj_alpha12, A02_l, s01 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, s01 );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  g0,  inc_g,
                  s01, rs_S );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  n_behind,
                  &minus_conj_alpha12,
                  A02_l, rs_A,
                  s01,   rs_S );
      bl1_zinvscalv( BLIS1_CONJUGATE,
                     n_behind,
                     &psi11_minus_alpha12,
                     s01, rs_S );

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
                y21, rs_Y,
                v21, inc_v,
                &beta );
      bl1_zscals( &minus_inv_tau11, &beta );

      // FLA_Copy( w21, z21 );
      // FLA_Axpy( minus_conj_alpha12, a22l, z21 );
      // FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, z21 );
      // FLA_Axpy( beta, u21, z21 );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, rs_Z );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_alpha12,
                  a22l, inc_al,
                  z21,  rs_Z );
      bl1_zinvscalv( BLIS1_CONJUGATE,
                     m_ahead,
                     &psi11_minus_alpha12,
                     z21, rs_Z );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  u21, inc_u,
                  z21, rs_Z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11,   y21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_zinvscalv( BLIS1_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, rs_Y );
      bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, rs_Z );
    }
	else // if ( n_ahead == 0 )
    {
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
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &w );
  // FLA_Obj_free( &al );
  // FLA_Obj_free( &ap );
  // FLA_Obj_free( &u );
  // FLA_Obj_free( &up );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  // FLA_Obj_free( &g );
  FLA_free( buff_w );
  FLA_free( buff_al );
  FLA_free( buff_ap );
  FLA_free( buff_u );
  FLA_free( buff_up );
  FLA_free( buff_v );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );
  FLA_free( buff_g );

  return FLA_SUCCESS;
}

