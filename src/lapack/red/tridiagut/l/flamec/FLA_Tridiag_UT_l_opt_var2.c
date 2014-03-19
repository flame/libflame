/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_l_opt_var2( FLA_Obj A, FLA_Obj T )
{
  return FLA_Tridiag_UT_l_step_opt_var2( A, T );
}

FLA_Error FLA_Tridiag_UT_l_step_opt_var2( FLA_Obj A, FLA_Obj T )
{
  FLA_Datatype datatype;
  int          m_A, m_T;
  int          rs_A, cs_A;
  int          rs_T, cs_T;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  m_T      = FLA_Obj_length( T );

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

      FLA_Tridiag_UT_l_step_ops_var2( m_A,
                                      m_T,
                                      buff_A, rs_A, cs_A,
                                      buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_Tridiag_UT_l_step_opd_var2( m_A,
                                      m_T,
                                      buff_A, rs_A, cs_A,
                                      buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_Tridiag_UT_l_step_opc_var2( m_A,
                                      m_T,
                                      buff_A, rs_A, cs_A,
                                      buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_Tridiag_UT_l_step_opz_var2( m_A,
                                      m_T,
                                      buff_A, rs_A, cs_A,
                                      buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Tridiag_UT_l_step_ops_var2( int m_A,
                                          int m_T,
                                          float* buff_A, int rs_A, int cs_A, 
                                          float* buff_T, int rs_T, int cs_T )
{
  float*    buff_2  = FLA_FLOAT_PTR( FLA_TWO );
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );

  float     first_elem;
  float     beta;
  float     inv_tau11;
  float     minus_inv_tau11;
  float     minus_upsilon11, minus_conj_upsilon11;
  float     minus_zeta11, minus_conj_zeta11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  float*    buff_u = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_z = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_w = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u  = 1;
  int       inc_z  = 1;
  int       inc_w  = 1;

  // Initialize some variables (only to prevent compiler warnings).
  first_elem      = *buff_0;
  minus_inv_tau11 = *buff_0;

  for ( i = 0; i < b_alg; ++i )
  {
    float*    A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    float*    tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    float*    upsilon11= buff_u + (i  )*inc_u;
    float*    u21      = buff_u + (i+1)*inc_u;

    float*    zeta11   = buff_z + (i  )*inc_z;
    float*    z21      = buff_z + (i+1)*inc_z;

    float*    w21      = buff_w + (i+1)*inc_w;

    float*    a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    float*    a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    int       m_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon11, minus_upsilon11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon11 );
      // FLA_Copy( minus_upsilon11, minus_conj_upsilon11 );
      bl1_smult3( buff_m1, upsilon11, &minus_upsilon11 );
      bl1_scopyconj( &minus_upsilon11, &minus_conj_upsilon11 );

      // FLA_Copy( zeta11, minus_zeta11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta11 );
      // FLA_Copy( minus_zeta11, minus_conj_zeta11 );
      bl1_smult3( buff_m1, zeta11, &minus_zeta11 );
      bl1_scopyconj( &minus_zeta11, &minus_conj_zeta11 );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon11, zeta11,    alpha11 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_zeta11,    upsilon11, alpha11 );
      bl1_saxpyv( BLIS1_CONJUGATE,
                  1,
                  &minus_upsilon11,
                  zeta11,  1,
                  alpha11, 1 );
      bl1_saxpyv( BLIS1_CONJUGATE,
                  1,
                  &minus_zeta11,
                  upsilon11, 1,
                  alpha11,  1 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_zeta11,    u21, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_upsilon11, z21, a21 );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_zeta11,
                  u21, inc_u,
                  a21, rs_A );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_upsilon11,
                  z21, inc_z,
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
    }

    if ( m_behind > 0 )
    {
      // FLA_Her2( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, u21, z21, A22 );
      bl1_ssyr2( BLIS1_LOWER_TRIANGULAR,
                 m_ahead,
                 buff_m1,
                 u21, inc_u,
                 z21, inc_z,
                 A22, rs_A, cs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Hemv( FLA_LOWER_TRIANGULAR, FLA_ONE, A22, a21, FLA_ZERO, w21 );
      bl1_ssymv( BLIS1_LOWER_TRIANGULAR,
                 m_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 w21, inc_w );

      // FLA_Copy( a21, u21 );
      // FLA_Copy( w21, z21 );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u21, inc_u );
      bl1_scopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      bl1_sdot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, inc_z,
                &beta );
      bl1_sinvscals( buff_2, &beta );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bl1_sscals( &minus_inv_tau11, &beta );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, inc_z );
      bl1_sscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, inc_z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );
      bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 t01, rs_T );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    if ( m_behind + 1 == b_alg && m_ahead > 0 )
    {
      // FLA_Her2( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, u21, z21, A22 );
      bl1_ssyr2( BLIS1_LOWER_TRIANGULAR,
                 m_ahead,
                 buff_m1,
                 u21, inc_u,
                 z21, inc_z,
                 A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &z );
  // FLA_Obj_free( &w );
  FLA_free( buff_u );
  FLA_free( buff_z );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}



FLA_Error FLA_Tridiag_UT_l_step_opd_var2( int m_A,
                                          int m_T,
                                          double* buff_A, int rs_A, int cs_A, 
                                          double* buff_T, int rs_T, int cs_T )
{
  double*   buff_2  = FLA_DOUBLE_PTR( FLA_TWO );
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );

  double    first_elem;
  double    beta;
  double    inv_tau11;
  double    minus_inv_tau11;
  double    minus_upsilon11, minus_conj_upsilon11;
  double    minus_zeta11, minus_conj_zeta11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  double*   buff_u = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_z = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_w = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u  = 1;
  int       inc_z  = 1;
  int       inc_w  = 1;

  // Initialize some variables (only to prevent compiler warnings).
  first_elem      = *buff_0;
  minus_inv_tau11 = *buff_0;

  for ( i = 0; i < b_alg; ++i )
  {
    double*   A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    double*   tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    double*   upsilon11= buff_u + (i  )*inc_u;
    double*   u21      = buff_u + (i+1)*inc_u;

    double*   zeta11   = buff_z + (i  )*inc_z;
    double*   z21      = buff_z + (i+1)*inc_z;

    double*   w21      = buff_w + (i+1)*inc_w;

    double*   a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    double*   a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    int       m_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon11, minus_upsilon11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon11 );
      // FLA_Copy( minus_upsilon11, minus_conj_upsilon11 );
      bl1_dmult3( buff_m1, upsilon11, &minus_upsilon11 );
      bl1_dcopyconj( &minus_upsilon11, &minus_conj_upsilon11 );

      // FLA_Copy( zeta11, minus_zeta11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta11 );
      // FLA_Copy( minus_zeta11, minus_conj_zeta11 );
      bl1_dmult3( buff_m1, zeta11, &minus_zeta11 );
      bl1_dcopyconj( &minus_zeta11, &minus_conj_zeta11 );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon11, zeta11,    alpha11 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_zeta11,    upsilon11, alpha11 );
      bl1_daxpyv( BLIS1_CONJUGATE,
                  1,
                  &minus_upsilon11,
                  zeta11,  1,
                  alpha11, 1 );
      bl1_daxpyv( BLIS1_CONJUGATE,
                  1,
                  &minus_zeta11,
                  upsilon11, 1,
                  alpha11,  1 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_zeta11,    u21, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_upsilon11, z21, a21 );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_zeta11,
                  u21, inc_u,
                  a21, rs_A );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_upsilon11,
                  z21, inc_z,
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
    }

    if ( m_behind > 0 )
    {
      // FLA_Her2( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, u21, z21, A22 );
      bl1_dsyr2( BLIS1_LOWER_TRIANGULAR,
                 m_ahead,
                 buff_m1,
                 u21, inc_u,
                 z21, inc_z,
                 A22, rs_A, cs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Hemv( FLA_LOWER_TRIANGULAR, FLA_ONE, A22, a21, FLA_ZERO, w21 );
      bl1_dsymv( BLIS1_LOWER_TRIANGULAR,
                 m_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 w21, inc_w );

      // FLA_Copy( a21, u21 );
      // FLA_Copy( w21, z21 );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u21, inc_u );
      bl1_dcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      bl1_ddot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, inc_z,
                &beta );
      bl1_dinvscals( buff_2, &beta );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bl1_dscals( &minus_inv_tau11, &beta );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, inc_z );
      bl1_dscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, inc_z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );
      bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 t01, rs_T );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    if ( m_behind + 1 == b_alg && m_ahead > 0 )
    {
      // FLA_Her2( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, u21, z21, A22 );
      bl1_dsyr2( BLIS1_LOWER_TRIANGULAR,
                 m_ahead,
                 buff_m1,
                 u21, inc_u,
                 z21, inc_z,
                 A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &z );
  // FLA_Obj_free( &w );
  FLA_free( buff_u );
  FLA_free( buff_z );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}



FLA_Error FLA_Tridiag_UT_l_step_opc_var2( int m_A,
                                          int m_T,
                                          scomplex* buff_A, int rs_A, int cs_A, 
                                          scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_2  = FLA_COMPLEX_PTR( FLA_TWO );
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );

  scomplex  first_elem;
  scomplex  beta;
  scomplex  inv_tau11;
  scomplex  minus_inv_tau11;
  scomplex  minus_upsilon11, minus_conj_upsilon11;
  scomplex  minus_zeta11, minus_conj_zeta11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  scomplex* buff_u = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_z = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_w = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u  = 1;
  int       inc_z  = 1;
  int       inc_w  = 1;

  // Initialize some variables (only to prevent compiler warnings).
  first_elem      = *buff_0;
  minus_inv_tau11 = *buff_0;

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    scomplex* upsilon11= buff_u + (i  )*inc_u;
    scomplex* u21      = buff_u + (i+1)*inc_u;

    scomplex* zeta11   = buff_z + (i  )*inc_z;
    scomplex* z21      = buff_z + (i+1)*inc_z;

    scomplex* w21      = buff_w + (i+1)*inc_w;

    scomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    scomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    int       m_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon11, minus_upsilon11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon11 );
      // FLA_Copy( minus_upsilon11, minus_conj_upsilon11 );
      bl1_cmult3( buff_m1, upsilon11, &minus_upsilon11 );
      bl1_ccopyconj( &minus_upsilon11, &minus_conj_upsilon11 );

      // FLA_Copy( zeta11, minus_zeta11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta11 );
      // FLA_Copy( minus_zeta11, minus_conj_zeta11 );
      bl1_cmult3( buff_m1, zeta11, &minus_zeta11 );
      bl1_ccopyconj( &minus_zeta11, &minus_conj_zeta11 );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon11, zeta11,    alpha11 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_zeta11,    upsilon11, alpha11 );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  1,
                  &minus_upsilon11,
                  zeta11,  1,
                  alpha11, 1 );
      bl1_caxpyv( BLIS1_CONJUGATE,
                  1,
                  &minus_zeta11,
                  upsilon11, 1,
                  alpha11,  1 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_zeta11,    u21, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_upsilon11, z21, a21 );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_zeta11,
                  u21, inc_u,
                  a21, rs_A );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_upsilon11,
                  z21, inc_z,
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
    }

    if ( m_behind > 0 )
    {
      // FLA_Her2( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, u21, z21, A22 );
      bl1_cher2( BLIS1_LOWER_TRIANGULAR,
                 BLIS1_NO_TRANSPOSE,
                 m_ahead,
                 buff_m1,
                 u21, inc_u,
                 z21, inc_z,
                 A22, rs_A, cs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Hemv( FLA_LOWER_TRIANGULAR, FLA_ONE, A22, a21, FLA_ZERO, w21 );
      bl1_chemv( BLIS1_LOWER_TRIANGULAR,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 w21, inc_w );

      // FLA_Copy( a21, u21 );
      // FLA_Copy( w21, z21 );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u21, inc_u );
      bl1_ccopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      bl1_cdot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, inc_z,
                &beta );
      bl1_cinvscals( buff_2, &beta );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bl1_cscals( &minus_inv_tau11, &beta );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, inc_z );
      bl1_cscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, inc_z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );
      bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 t01, rs_T );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    if ( m_behind + 1 == b_alg && m_ahead > 0 )
    {
      // FLA_Her2( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, u21, z21, A22 );
      bl1_cher2( BLIS1_LOWER_TRIANGULAR,
                 BLIS1_NO_TRANSPOSE,
                 m_ahead,
                 buff_m1,
                 u21, inc_u,
                 z21, inc_z,
                 A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &z );
  // FLA_Obj_free( &w );
  FLA_free( buff_u );
  FLA_free( buff_z );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}



FLA_Error FLA_Tridiag_UT_l_step_opz_var2( int m_A,
                                          int m_T,
                                          dcomplex* buff_A, int rs_A, int cs_A, 
                                          dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_2  = FLA_DOUBLE_COMPLEX_PTR( FLA_TWO );
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );

  dcomplex  first_elem;
  dcomplex  beta;
  dcomplex  inv_tau11;
  dcomplex  minus_inv_tau11;
  dcomplex  minus_upsilon11, minus_conj_upsilon11;
  dcomplex  minus_zeta11, minus_conj_zeta11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  dcomplex* buff_u = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_z = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_w = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u  = 1;
  int       inc_z  = 1;
  int       inc_w  = 1;

  // Initialize some variables (only to prevent compiler warnings).
  first_elem      = *buff_0;
  minus_inv_tau11 = *buff_0;

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    dcomplex* upsilon11= buff_u + (i  )*inc_u;
    dcomplex* u21      = buff_u + (i+1)*inc_u;

    dcomplex* zeta11   = buff_z + (i  )*inc_z;
    dcomplex* z21      = buff_z + (i+1)*inc_z;

    dcomplex* w21      = buff_w + (i+1)*inc_w;

    dcomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    int       m_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon11, minus_upsilon11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon11 );
      // FLA_Copy( minus_upsilon11, minus_conj_upsilon11 );
      bl1_zmult3( buff_m1, upsilon11, &minus_upsilon11 );
      bl1_zcopyconj( &minus_upsilon11, &minus_conj_upsilon11 );

      // FLA_Copy( zeta11, minus_zeta11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta11 );
      // FLA_Copy( minus_zeta11, minus_conj_zeta11 );
      bl1_zmult3( buff_m1, zeta11, &minus_zeta11 );
      bl1_zcopyconj( &minus_zeta11, &minus_conj_zeta11 );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon11, zeta11,    alpha11 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_zeta11,    upsilon11, alpha11 );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  1,
                  &minus_upsilon11,
                  zeta11,  1,
                  alpha11, 1 );
      bl1_zaxpyv( BLIS1_CONJUGATE,
                  1,
                  &minus_zeta11,
                  upsilon11, 1,
                  alpha11,  1 );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_zeta11,    u21, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_upsilon11, z21, a21 );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_zeta11,
                  u21, inc_u,
                  a21, rs_A );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_upsilon11,
                  z21, inc_z,
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
    }

    if ( m_behind > 0 )
    {
      // FLA_Her2( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, u21, z21, A22 );
      bl1_zher2( BLIS1_LOWER_TRIANGULAR,
                 BLIS1_NO_TRANSPOSE,
                 m_ahead,
                 buff_m1,
                 u21, inc_u,
                 z21, inc_z,
                 A22, rs_A, cs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Hemv( FLA_LOWER_TRIANGULAR, FLA_ONE, A22, a21, FLA_ZERO, w21 );
      bl1_zhemv( BLIS1_LOWER_TRIANGULAR,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 w21, inc_w );

      // FLA_Copy( a21, u21 );
      // FLA_Copy( w21, z21 );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u21, inc_u );
      bl1_zcopyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  w21, inc_w,
                  z21, inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      bl1_zdot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, inc_z,
                &beta );
      bl1_zinvscals( buff_2, &beta );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bl1_zscals( &minus_inv_tau11, &beta );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, inc_z );
      bl1_zscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, inc_z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );
      bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 t01, rs_T );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    if ( m_behind + 1 == b_alg && m_ahead > 0 )
    {
      // FLA_Her2( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, u21, z21, A22 );
      bl1_zher2( BLIS1_LOWER_TRIANGULAR,
                 BLIS1_NO_TRANSPOSE,
                 m_ahead,
                 buff_m1,
                 u21, inc_u,
                 z21, inc_z,
                 A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &z );
  // FLA_Obj_free( &w );
  FLA_free( buff_u );
  FLA_free( buff_z );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}

