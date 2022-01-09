/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_ofu_var2( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  return FLA_Bidiag_UT_u_step_ofu_var2( A, TU, TV );
}

FLA_Error FLA_Bidiag_UT_u_step_ofu_var2( FLA_Obj A, FLA_Obj T, FLA_Obj S )
{
  FLA_Datatype datatype;
  integer          m_A, n_A, m_TS;
  integer          rs_A, cs_A;
  integer          rs_T, cs_T;
  integer          rs_S, cs_S;

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

      FLA_Bidiag_UT_u_step_ofs_var2( m_A,
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

      FLA_Bidiag_UT_u_step_ofd_var2( m_A,
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

      FLA_Bidiag_UT_u_step_ofc_var2( m_A,
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

      FLA_Bidiag_UT_u_step_ofz_var2( m_A,
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



FLA_Error FLA_Bidiag_UT_u_step_ofs_var2( integer m_A,
                                         integer n_A,
                                         integer m_TS,
                                         float* buff_A, integer rs_A, integer cs_A, 
                                         float* buff_T, integer rs_T, integer cs_T, 
                                         float* buff_S, integer rs_S, integer cs_S )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );

  float     beta;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  float*    buff_v = ( float* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_y = ( float* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_z = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_v  = 1;
  integer       inc_y  = 1;
  integer       inc_z  = 1;

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

    float*    y21      = buff_y + (i+1)*inc_y;

    float*    z21      = buff_z + (i+1)*inc_z;

    float*    a12t_l   = a12t   + (0  )*cs_A + (0  )*rs_A;
    float*    a12t_r   = a12t   + (1  )*cs_A + (0  )*rs_A;

    float*    v21_t    = v21    + (0  )*inc_v;
    float*    v21_b    = v21    + (1  )*inc_v;

    integer       m_ahead  = m_A - i - 1;
    integer       n_ahead  = n_A - i - 1;
    integer       m_behind = i;
    integer       n_behind = i;

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
      // FLA_Copyt( FLA_TRANSPOSE, a12t, y21 );
      // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a21, FLA_ONE, y21 );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  y21,  inc_y );
      bl1_sgemv( BLIS1_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_1,
                 y21, inc_y );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, y21 );
      bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, inc_y );

      // FLA_Axpyt( FLA_TRANSPOSE, FLA_MINUS_ONE, y21, a12t );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  buff_m1,
                  y21,  inc_y,
                  a12t, cs_A );

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
                  v21_b,  inc_y );

      // FLA_Dotc( FLA_CONJUGATE, v21, y21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      bl1_sdot( BLIS1_CONJUGATE,
                n_ahead,
                v21, inc_v,
                y21, inc_y,
                &beta );
      bl1_sneg1( &beta );

      // FLA_Copy( a21, z21 );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, v21, beta, z21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  z21, inc_z );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 v21, inc_v,
                 &beta,
                 z21, inc_z );
      bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, inc_z );

      // FLA_Ger( FLA_MINUS_ONE, a21, y21, A22 );
      // FLA_Ger( FLA_MINUS_ONE, z21, v21, A22 );
      FLA_Fused_Gerc2_ops_var1( m_ahead,
                                n_ahead,
                                buff_m1,
                                a21, rs_A,
                                y21, inc_y,
                                z21, inc_z,
                                v21, inc_v,
                                A22, rs_A, cs_A );

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
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_v );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_ofd_var2( integer m_A,
                                         integer n_A,
                                         integer m_TS,
                                         double* buff_A, integer rs_A, integer cs_A, 
                                         double* buff_T, integer rs_T, integer cs_T, 
                                         double* buff_S, integer rs_S, integer cs_S )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );

  double    beta;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  double*   buff_v = ( double* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_y = ( double* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_z = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_v  = 1;
  integer       inc_y  = 1;
  integer       inc_z  = 1;

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

    double*   y21      = buff_y + (i+1)*inc_y;

    double*   z21      = buff_z + (i+1)*inc_z;

    double*   a12t_l   = a12t   + (0  )*cs_A + (0  )*rs_A;
    double*   a12t_r   = a12t   + (1  )*cs_A + (0  )*rs_A;

    double*   v21_t    = v21    + (0  )*inc_v;
    double*   v21_b    = v21    + (1  )*inc_v;

    integer       m_ahead  = m_A - i - 1;
    integer       n_ahead  = n_A - i - 1;
    integer       m_behind = i;
    integer       n_behind = i;

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
      // FLA_Copyt( FLA_TRANSPOSE, a12t, y21 );
      // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a21, FLA_ONE, y21 );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  y21,  inc_y );
      bl1_dgemv( BLIS1_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_1,
                 y21, inc_y );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, y21 );
      bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, inc_y );

      // FLA_Axpyt( FLA_TRANSPOSE, FLA_MINUS_ONE, y21, a12t );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  n_ahead,
                  buff_m1,
                  y21,  inc_y,
                  a12t, cs_A );

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
                  v21_b,  inc_y );

      // FLA_Dotc( FLA_CONJUGATE, v21, y21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      bl1_ddot( BLIS1_CONJUGATE,
                n_ahead,
                v21, inc_v,
                y21, inc_y,
                &beta );
      bl1_dneg1( &beta );

      // FLA_Copy( a21, z21 );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, v21, beta, z21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  z21, inc_z );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 v21, inc_v,
                 &beta,
                 z21, inc_z );
      bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, inc_z );

      // FLA_Ger( FLA_MINUS_ONE, a21, y21, A22 );
      // FLA_Ger( FLA_MINUS_ONE, z21, v21, A22 );
      FLA_Fused_Gerc2_opd_var1( m_ahead,
                                n_ahead,
                                buff_m1,
                                a21, rs_A,
                                y21, inc_y,
                                z21, inc_z,
                                v21, inc_v,
                                A22, rs_A, cs_A );

      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A02, v21, FLA_ZERO, s01 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_CONJUGATE,
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
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_v );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_ofc_var2( integer m_A,
                                         integer n_A,
                                         integer m_TS,
                                         scomplex* buff_A, integer rs_A, integer cs_A, 
                                         scomplex* buff_T, integer rs_T, integer cs_T, 
                                         scomplex* buff_S, integer rs_S, integer cs_S )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );

  scomplex  beta;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  scomplex* buff_v = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_y = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_z = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_v  = 1;
  integer       inc_y  = 1;
  integer       inc_z  = 1;

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

    scomplex* y21      = buff_y + (i+1)*inc_y;

    scomplex* z21      = buff_z + (i+1)*inc_z;

    scomplex* a12t_l   = a12t   + (0  )*cs_A + (0  )*rs_A;
    scomplex* a12t_r   = a12t   + (1  )*cs_A + (0  )*rs_A;

    scomplex* v21_t    = v21    + (0  )*inc_v;
    scomplex* v21_b    = v21    + (1  )*inc_v;

    integer       m_ahead  = m_A - i - 1;
    integer       n_ahead  = n_A - i - 1;
    integer       m_behind = i;
    integer       n_behind = i;

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
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a12t, y21 );
      // FLA_Gemvc( FLA_CONJ_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, a21, FLA_ONE, y21 );
      bl1_ccopyv( BLIS1_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  y21,  inc_y );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_1,
                 y21, inc_y );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, y21 );
      bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, inc_y );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_MINUS_ONE, y21, a12t );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  buff_m1,
                  y21,  inc_y,
                  a12t, cs_A );

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
                  v21_b,  inc_y );

      // FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      bl1_cdot( BLIS1_CONJUGATE,
                n_ahead,
                y21, inc_y,
                v21, inc_v,
                &beta );
      bl1_cneg1( &beta );

      // FLA_Copy( a21, z21 );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, v21, beta, z21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  z21, inc_z );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 v21, inc_v,
                 &beta,
                 z21, inc_z );
      bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, inc_z );

      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, a21, y21, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );
      FLA_Fused_Gerc2_opc_var1( m_ahead,
                                n_ahead,
                                buff_m1,
                                a21, rs_A,
                                y21, inc_y,
                                z21, inc_z,
                                v21, inc_v,
                                A22, rs_A, cs_A );

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
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_v );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_ofz_var2( integer m_A,
                                         integer n_A,
                                         integer m_TS,
                                         dcomplex* buff_A, integer rs_A, integer cs_A, 
                                         dcomplex* buff_T, integer rs_T, integer cs_T, 
                                         dcomplex* buff_S, integer rs_S, integer cs_S )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );

  dcomplex  beta;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  dcomplex* buff_v = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_y = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_z = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_v  = 1;
  integer       inc_y  = 1;
  integer       inc_z  = 1;

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

    dcomplex* y21      = buff_y + (i+1)*inc_y;

    dcomplex* z21      = buff_z + (i+1)*inc_z;

    dcomplex* a12t_l   = a12t   + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a12t_r   = a12t   + (1  )*cs_A + (0  )*rs_A;

    dcomplex* v21_t    = v21    + (0  )*inc_v;
    dcomplex* v21_b    = v21    + (1  )*inc_v;

    integer       m_ahead  = m_A - i - 1;
    integer       n_ahead  = n_A - i - 1;
    integer       m_behind = i;
    integer       n_behind = i;

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
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a12t, y21 );
      // FLA_Gemvc( FLA_CONJ_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, a21, FLA_ONE, y21 );
      bl1_zcopyv( BLIS1_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  y21,  inc_y );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_1,
                 y21, inc_y );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, y21 );
      bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, inc_y );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_MINUS_ONE, y21, a12t );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  n_ahead,
                  buff_m1,
                  y21,  inc_y,
                  a12t, cs_A );

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
                  v21_b,  inc_y );

      // FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      bl1_zdot( BLIS1_CONJUGATE,
                n_ahead,
                y21, inc_y,
                v21, inc_v,
                &beta );
      bl1_zneg1( &beta );

      // FLA_Copy( a21, z21 );
      // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, v21, beta, z21 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  z21, inc_z );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 v21, inc_v,
                 &beta,
                 z21, inc_z );
      bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, inc_z );

      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, a21, y21, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );
      FLA_Fused_Gerc2_opz_var1( m_ahead,
                                n_ahead,
                                buff_m1,
                                a21, rs_A,
                                y21, inc_y,
                                z21, inc_z,
                                v21, inc_v,
                                A22, rs_A, cs_A );

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
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_v );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}

