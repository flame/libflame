/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Hess_UT_opt_var4( FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;
  FLA_Obj   Y, Z;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Y );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Z );

  r_val = FLA_Hess_UT_step_opt_var4( A, Y, Z, T );

  FLA_Obj_free( &Y );
  FLA_Obj_free( &Z );

  return r_val;
}

FLA_Error FLA_Hess_UT_step_opt_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T )
{
  FLA_Datatype datatype;
  int          m_A, m_T;
  int          rs_A, cs_A;
  int          rs_Y, cs_Y;
  int          rs_Z, cs_Z;
  int          rs_T, cs_T;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  m_T      = FLA_Obj_length( T );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_Y     = FLA_Obj_row_stride( Y );
  cs_Y     = FLA_Obj_col_stride( Y );

  rs_Z     = FLA_Obj_row_stride( Z );
  cs_Z     = FLA_Obj_col_stride( Z );

  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_Y = FLA_FLOAT_PTR( Y );
      float* buff_Z = FLA_FLOAT_PTR( Z );
      float* buff_T = FLA_FLOAT_PTR( T );

      FLA_Hess_UT_step_ops_var4( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_Y, rs_Y, cs_Y,
                                 buff_Z, rs_Z, cs_Z,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_Y = FLA_DOUBLE_PTR( Y );
      double* buff_Z = FLA_DOUBLE_PTR( Z );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_Hess_UT_step_opd_var4( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_Y, rs_Y, cs_Y,
                                 buff_Z, rs_Z, cs_Z,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_Y = FLA_COMPLEX_PTR( Y );
      scomplex* buff_Z = FLA_COMPLEX_PTR( Z );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_Hess_UT_step_opc_var4( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_Y, rs_Y, cs_Y,
                                 buff_Z, rs_Z, cs_Z,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_Y = FLA_DOUBLE_COMPLEX_PTR( Y );
      dcomplex* buff_Z = FLA_DOUBLE_COMPLEX_PTR( Z );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_Hess_UT_step_opz_var4( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_Y, rs_Y, cs_Y,
                                 buff_Z, rs_Z, cs_Z,
                                 buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_ops_var4( int m_A,
                                     int m_T,
                                     float* buff_A, int rs_A, int cs_A, 
                                     float* buff_Y, int rs_Y, int cs_Y,
                                     float* buff_Z, int rs_Z, int cs_Z,
                                     float* buff_T, int rs_T, int cs_T )
{
  float*    buff_2  = FLA_FLOAT_PTR( FLA_TWO );
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );

  float     first_elem, last_elem;
  float     dot_product;
  float     beta, conj_beta;
  float     inv_tau11;
  float     minus_inv_tau11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  float*    buff_d = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_e = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_f = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_d  = 1;
  int       inc_e  = 1;
  int       inc_f  = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bl1_ssetm( m_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bl1_ssetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    float*    a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    y10t     = buff_Y + (0  )*cs_Y + (i  )*rs_Y;
    float*    Y20      = buff_Y + (0  )*cs_Y + (i+1)*rs_Y;
    float*    y21      = buff_Y + (i  )*cs_Y + (i+1)*rs_Y;

    float*    z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    float*    Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    float*    z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    float*    t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    float*    tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    float*    d0       = buff_d + (0  )*inc_d;

    float*    e0       = buff_e + (0  )*inc_e;

    float*    f0       = buff_f + (0  )*inc_f;

    float*    a10t_r   = a10t   + (i-1)*cs_A + (0  )*rs_A;

    float*    a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    float*    a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    float*    ABL      = a10t;
    float*    ZBL      = z10t;

    float*    a2       = alpha11;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a10t_r, last_elem );
      // FLA_Set( FLA_ONE, a10t_r );
      last_elem = *a10t_r;
      *a10t_r = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a10t, FLA_ONE, a2 );
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
               ZBL,  rs_Z, cs_Z,
               a10t, cs_A,
               buff_1,
               a2,   rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, A20, z10t, FLA_ONE, a12t );
    bl1_sgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bl1_sgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_m1,
               A20,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a10t_r );
      *a10t_r = last_elem;
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

      // FLA_Set( FLA_ONE, inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      // FLA_Copy( inv_tau11, minus_inv_tau11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );
      bl1_sdiv3( buff_1, tau11, &inv_tau11 );
      bl1_sneg2( &inv_tau11, &minus_inv_tau11 );

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, y21 );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y21, rs_Y );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, z21 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Y20, a21, FLA_ZERO, e0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, a21, FLA_ZERO, f0 );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 d0,  inc_d );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Y20, rs_Y, cs_Y,
                 a21, rs_A,
                 buff_0,
                 e0,  inc_e );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20, rs_Z, cs_Z,
                 a21, rs_A,
                 buff_0,
                 f0,  inc_f );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, y21 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, e0, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, d0, FLA_ONE, z21 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 z21, rs_Z );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 d0,  inc_d,
                 buff_1,
                 z21, rs_Z );

      // FLA_Copy( d0, t01 );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bl1_sdot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, rs_Z,
                &beta );
      bl1_sinvscals( buff_2, &beta );
      bl1_scopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y21 );
      // FLA_Scal( inv_tau11, y21 );
      bl1_sscals( &minus_inv_tau11, &conj_beta );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y21, rs_Y );
      bl1_sscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y21, rs_Y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bl1_sscals( &minus_inv_tau11, &beta );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, rs_Z );
      bl1_sscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, rs_Z );

      // FLA_Dot( a12t, a21, dot_product );
      // FLA_Scal( minus_inv_tau11, dot_product );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, dot_product, a21, a12t );
      bl1_sdot( BLIS1_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                a21,  rs_A,
                &dot_product );
      bl1_sscals( &minus_inv_tau11, &dot_product );
      bl1_saxpyv( BLIS1_CONJUGATE,
                  m_ahead,
                  &dot_product,
                  a21,  rs_A,
                  a12t, cs_A );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, e0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, e0, a21, A02 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 e0,  inc_e );
      bl1_sger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                e0,  inc_e,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opd_var4( int m_A,
                                     int m_T,
                                     double* buff_A, int rs_A, int cs_A, 
                                     double* buff_Y, int rs_Y, int cs_Y,
                                     double* buff_Z, int rs_Z, int cs_Z,
                                     double* buff_T, int rs_T, int cs_T )
{
  double*   buff_2  = FLA_DOUBLE_PTR( FLA_TWO );
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );

  double    first_elem, last_elem;
  double    dot_product;
  double    beta, conj_beta;
  double    inv_tau11;
  double    minus_inv_tau11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  double*   buff_d = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_e = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_f = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_d  = 1;
  int       inc_e  = 1;
  int       inc_f  = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bl1_dsetm( m_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bl1_dsetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    double*   a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   y10t     = buff_Y + (0  )*cs_Y + (i  )*rs_Y;
    double*   Y20      = buff_Y + (0  )*cs_Y + (i+1)*rs_Y;
    double*   y21      = buff_Y + (i  )*cs_Y + (i+1)*rs_Y;

    double*   z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    double*   Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    double*   z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    double*   t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    double*   tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    double*   d0       = buff_d + (0  )*inc_d;

    double*   e0       = buff_e + (0  )*inc_e;

    double*   f0       = buff_f + (0  )*inc_f;

    double*   a10t_r   = a10t   + (i-1)*cs_A + (0  )*rs_A;

    double*   a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    double*   a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    double*   ABL      = a10t;
    double*   ZBL      = z10t;

    double*   a2       = alpha11;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a10t_r, last_elem );
      // FLA_Set( FLA_ONE, a10t_r );
      last_elem = *a10t_r;
      *a10t_r = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a10t, FLA_ONE, a2 );
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
               ZBL,  rs_Z, cs_Z,
               a10t, cs_A,
               buff_1,
               a2,   rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, A20, z10t, FLA_ONE, a12t );
    bl1_dgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bl1_dgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_m1,
               A20,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a10t_r );
      *a10t_r = last_elem;
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

      // FLA_Set( FLA_ONE, inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      // FLA_Copy( inv_tau11, minus_inv_tau11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );
      bl1_ddiv3( buff_1, tau11, &inv_tau11 );
      bl1_dneg2( &inv_tau11, &minus_inv_tau11 );

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, y21 );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y21, rs_Y );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, z21 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Y20, a21, FLA_ZERO, e0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, a21, FLA_ZERO, f0 );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 d0,  inc_d );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Y20, rs_Y, cs_Y,
                 a21, rs_A,
                 buff_0,
                 e0,  inc_e );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20, rs_Z, cs_Z,
                 a21, rs_A,
                 buff_0,
                 f0,  inc_f );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, y21 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, e0, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, d0, FLA_ONE, z21 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 z21, rs_Z );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 d0,  inc_d,
                 buff_1,
                 z21, rs_Z );

      // FLA_Copy( d0, t01 );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bl1_ddot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, rs_Z,
                &beta );
      bl1_dinvscals( buff_2, &beta );
      bl1_dcopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y21 );
      // FLA_Scal( inv_tau11, y21 );
      bl1_dscals( &minus_inv_tau11, &conj_beta );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y21, rs_Y );
      bl1_dscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y21, rs_Y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bl1_dscals( &minus_inv_tau11, &beta );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, rs_Z );
      bl1_dscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, rs_Z );

      // FLA_Dot( a12t, a21, dot_product );
      // FLA_Scal( minus_inv_tau11, dot_product );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, dot_product, a21, a12t );
      bl1_ddot( BLIS1_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                a21,  rs_A,
                &dot_product );
      bl1_dscals( &minus_inv_tau11, &dot_product );
      bl1_daxpyv( BLIS1_CONJUGATE,
                  m_ahead,
                  &dot_product,
                  a21,  rs_A,
                  a12t, cs_A );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, e0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, e0, a21, A02 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 e0,  inc_e );
      bl1_dger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                e0,  inc_e,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opc_var4( int m_A,
                                     int m_T,
                                     scomplex* buff_A, int rs_A, int cs_A, 
                                     scomplex* buff_Y, int rs_Y, int cs_Y,
                                     scomplex* buff_Z, int rs_Z, int cs_Z,
                                     scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_2  = FLA_COMPLEX_PTR( FLA_TWO );
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );

  scomplex  first_elem, last_elem;
  scomplex  dot_product;
  scomplex  beta, conj_beta;
  scomplex  inv_tau11;
  scomplex  minus_inv_tau11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  scomplex* buff_d = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_e = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_f = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_d  = 1;
  int       inc_e  = 1;
  int       inc_f  = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bl1_csetm( m_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bl1_csetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* y10t     = buff_Y + (0  )*cs_Y + (i  )*rs_Y;
    scomplex* Y20      = buff_Y + (0  )*cs_Y + (i+1)*rs_Y;
    scomplex* y21      = buff_Y + (i  )*cs_Y + (i+1)*rs_Y;

    scomplex* z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    scomplex* Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    scomplex* z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    scomplex* d0       = buff_d + (0  )*inc_d;

    scomplex* e0       = buff_e + (0  )*inc_e;

    scomplex* f0       = buff_f + (0  )*inc_f;

    scomplex* a10t_r   = a10t   + (i-1)*cs_A + (0  )*rs_A;

    scomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    scomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    scomplex* ABL      = a10t;
    scomplex* ZBL      = z10t;

    scomplex* a2       = alpha11;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a10t_r, last_elem );
      // FLA_Set( FLA_ONE, a10t_r );
      last_elem = *a10t_r;
      *a10t_r = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a10t, FLA_ONE, a2 );
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
               ZBL,  rs_Z, cs_Z,
               a10t, cs_A,
               buff_1,
               a2,   rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, A20, z10t, FLA_ONE, a12t );
    bl1_cgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bl1_cgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_m1,
               A20,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a10t_r );
      *a10t_r = last_elem;
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

      // FLA_Set( FLA_ONE, inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      // FLA_Copy( inv_tau11, minus_inv_tau11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );
      bl1_cdiv3( buff_1, tau11, &inv_tau11 );
      bl1_cneg2( &inv_tau11, &minus_inv_tau11 );

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, y21 );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y21, rs_Y );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, z21 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Y20, a21, FLA_ZERO, e0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, a21, FLA_ZERO, f0 );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 d0,  inc_d );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Y20, rs_Y, cs_Y,
                 a21, rs_A,
                 buff_0,
                 e0,  inc_e );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20, rs_Z, cs_Z,
                 a21, rs_A,
                 buff_0,
                 f0,  inc_f );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, y21 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, e0, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, d0, FLA_ONE, z21 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 z21, rs_Z );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 d0,  inc_d,
                 buff_1,
                 z21, rs_Z );

      // FLA_Copy( d0, t01 );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bl1_cdot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, rs_Z,
                &beta );
      bl1_cinvscals( buff_2, &beta );
      bl1_ccopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y21 );
      // FLA_Scal( inv_tau11, y21 );
      bl1_cscals( &minus_inv_tau11, &conj_beta );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y21, rs_Y );
      bl1_cscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y21, rs_Y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bl1_cscals( &minus_inv_tau11, &beta );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, rs_Z );
      bl1_cscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, rs_Z );

      // FLA_Dot( a12t, a21, dot_product );
      // FLA_Scal( minus_inv_tau11, dot_product );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, dot_product, a21, a12t );
      bl1_cdot( BLIS1_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                a21,  rs_A,
                &dot_product );
      bl1_cscals( &minus_inv_tau11, &dot_product );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  m_ahead,
                  &dot_product,
                  a21,  rs_A,
                  a12t, cs_A );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, e0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, e0, a21, A02 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 e0,  inc_e );
      bl1_cger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                e0,  inc_e,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opz_var4( int m_A,
                                     int m_T,
                                     dcomplex* buff_A, int rs_A, int cs_A, 
                                     dcomplex* buff_Y, int rs_Y, int cs_Y,
                                     dcomplex* buff_Z, int rs_Z, int cs_Z,
                                     dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_2  = FLA_DOUBLE_COMPLEX_PTR( FLA_TWO );
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );

  dcomplex  first_elem, last_elem;
  dcomplex  dot_product;
  dcomplex  beta, conj_beta;
  dcomplex  inv_tau11;
  dcomplex  minus_inv_tau11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  dcomplex* buff_d = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_e = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_f = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_d  = 1;
  int       inc_e  = 1;
  int       inc_f  = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bl1_zsetm( m_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bl1_zsetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* y10t     = buff_Y + (0  )*cs_Y + (i  )*rs_Y;
    dcomplex* Y20      = buff_Y + (0  )*cs_Y + (i+1)*rs_Y;
    dcomplex* y21      = buff_Y + (i  )*cs_Y + (i+1)*rs_Y;

    dcomplex* z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    dcomplex* Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    dcomplex* z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    dcomplex* d0       = buff_d + (0  )*inc_d;

    dcomplex* e0       = buff_e + (0  )*inc_e;

    dcomplex* f0       = buff_f + (0  )*inc_f;

    dcomplex* a10t_r   = a10t   + (i-1)*cs_A + (0  )*rs_A;

    dcomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    dcomplex* ABL      = a10t;
    dcomplex* ZBL      = z10t;

    dcomplex* a2       = alpha11;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a10t_r, last_elem );
      // FLA_Set( FLA_ONE, a10t_r );
      last_elem = *a10t_r;
      *a10t_r = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a10t, FLA_ONE, a2 );
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
               ZBL,  rs_Z, cs_Z,
               a10t, cs_A,
               buff_1,
               a2,   rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, A20, z10t, FLA_ONE, a12t );
    bl1_zgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bl1_zgemv( BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_m1,
               A20,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a10t_r );
      *a10t_r = last_elem;
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

      // FLA_Set( FLA_ONE, inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      // FLA_Copy( inv_tau11, minus_inv_tau11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );
      bl1_zdiv3( buff_1, tau11, &inv_tau11 );
      bl1_zneg2( &inv_tau11, &minus_inv_tau11 );

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, y21 );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y21, rs_Y );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, z21 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Y20, a21, FLA_ZERO, e0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, a21, FLA_ZERO, f0 );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 d0,  inc_d );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Y20, rs_Y, cs_Y,
                 a21, rs_A,
                 buff_0,
                 e0,  inc_e );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20, rs_Z, cs_Z,
                 a21, rs_A,
                 buff_0,
                 f0,  inc_f );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, y21 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, e0, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, d0, FLA_ONE, z21 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 z21, rs_Z );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 d0,  inc_d,
                 buff_1,
                 z21, rs_Z );

      // FLA_Copy( d0, t01 );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  n_behind,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bl1_zdot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, rs_Z,
                &beta );
      bl1_zinvscals( buff_2, &beta );
      bl1_zcopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y21 );
      // FLA_Scal( inv_tau11, y21 );
      bl1_zscals( &minus_inv_tau11, &conj_beta );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y21, rs_Y );
      bl1_zscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y21, rs_Y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bl1_zscals( &minus_inv_tau11, &beta );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, rs_Z );
      bl1_zscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, rs_Z );

      // FLA_Dot( a12t, a21, dot_product );
      // FLA_Scal( minus_inv_tau11, dot_product );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, dot_product, a21, a12t );
      bl1_zdot( BLIS1_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                a21,  rs_A,
                &dot_product );
      bl1_zscals( &minus_inv_tau11, &dot_product );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  m_ahead,
                  &dot_product,
                  a21,  rs_A,
                  a12t, cs_A );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, e0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, e0, a21, A02 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 e0,  inc_e );
      bl1_zger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                e0,  inc_e,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );

  return FLA_SUCCESS;
}

