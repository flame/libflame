/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Hess_UT_ofu_var2( FLA_Obj A, FLA_Obj T )
{
  return FLA_Hess_UT_step_ofu_var2( A, T );
}

FLA_Error FLA_Hess_UT_step_ofu_var2( FLA_Obj A, FLA_Obj T )
{
  FLA_Datatype datatype;
  integer          m_A, m_T;
  integer          rs_A, cs_A;
  integer          rs_T, cs_T;

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

      FLA_Hess_UT_step_ofs_var2( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_Hess_UT_step_ofd_var2( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_Hess_UT_step_ofc_var2( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_Hess_UT_step_ofz_var2( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_ofs_var2( integer m_A,
                                     integer m_T,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_T, integer rs_T, integer cs_T )
{
  float*    buff_2  = FLA_FLOAT_PTR( FLA_TWO );
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );

  float     first_elem;
  float     dot_product;
  float     beta, conj_beta;
  float     inv_tau11;
  float     minus_inv_tau11;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  float*    buff_y = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_z = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_y  = 1;
  integer       inc_z  = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    float*    A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    float*    tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    float*    y0       = buff_y + (0  )*inc_y;
    float*    y2       = buff_y + (i+1)*inc_y;

    float*    z2       = buff_z + (i+1)*inc_z;

    float*    a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    float*    a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    integer       m_ahead  = m_A - i - 1;
    integer       n_ahead  = m_A - i - 1;
    integer       m_behind = i;
    integer       n_behind = i;

    /*------------------------------------------------------------*/

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

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, y2 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, z2 );
      FLA_Fused_Ahx_Ax_ops_var1( m_ahead,
                                 n_ahead,
                                 A22, rs_A, cs_A,
                                 a21, rs_A,
                                 y2,  inc_y,
                                 z2,  inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z2, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bl1_sdot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z2,  inc_z,
                &beta );
      bl1_sinvscals( buff_2, &beta );
      bl1_scopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y2 );
      // FLA_Scal( inv_tau11, y2 );
      bl1_sscals( &minus_inv_tau11, &conj_beta );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y2, inc_y );
      bl1_sscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y2, inc_y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z2 );
      // FLA_Scal( inv_tau11, z2 );
      bl1_sscals( &minus_inv_tau11, &beta );
      bl1_saxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z2, inc_z );
      bl1_sscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z2, inc_z );

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

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, y0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, y0, a21, A02 );
      bl1_sgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y0,  inc_y );
      bl1_sger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                y0,  inc_y,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, a21, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, a21, A22 );
      FLA_Fused_Gerc2_ops_var1( m_ahead,
                                n_ahead,
                                buff_m1,
                                a21, rs_A,
                                y2,  inc_y,
                                z2,  inc_z,
                                a21, rs_A,
                                A22, rs_A, cs_A );

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

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_ofd_var2( integer m_A,
                                     integer m_T,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_T, integer rs_T, integer cs_T )
{
  double*   buff_2  = FLA_DOUBLE_PTR( FLA_TWO );
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );

  double    first_elem;
  double    dot_product;
  double    beta, conj_beta;
  double    inv_tau11;
  double    minus_inv_tau11;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  double*   buff_y = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_z = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_y  = 1;
  integer       inc_z  = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    double*   A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    double*   tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    double*   y0       = buff_y + (0  )*inc_y;
    double*   y2       = buff_y + (i+1)*inc_y;

    double*   z2       = buff_z + (i+1)*inc_z;

    double*   a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    double*   a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    integer       m_ahead  = m_A - i - 1;
    integer       n_ahead  = m_A - i - 1;
    integer       m_behind = i;
    integer       n_behind = i;

    /*------------------------------------------------------------*/

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

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, y2 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, z2 );
      FLA_Fused_Ahx_Ax_opd_var1( m_ahead,
                                 n_ahead,
                                 A22, rs_A, cs_A,
                                 a21, rs_A,
                                 y2,  inc_y,
                                 z2,  inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z2, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bl1_ddot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z2,  inc_z,
                &beta );
      bl1_dinvscals( buff_2, &beta );
      bl1_dcopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y2 );
      // FLA_Scal( inv_tau11, y2 );
      bl1_dscals( &minus_inv_tau11, &conj_beta );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y2, inc_y );
      bl1_dscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y2, inc_y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z2 );
      // FLA_Scal( inv_tau11, z2 );
      bl1_dscals( &minus_inv_tau11, &beta );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z2, inc_z );
      bl1_dscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z2, inc_z );

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

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, y0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, y0, a21, A02 );
      bl1_dgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y0,  inc_y );
      bl1_dger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                y0,  inc_y,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, a21, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, a21, A22 );
      FLA_Fused_Gerc2_opd_var1( m_ahead,
                                n_ahead,
                                buff_m1,
                                a21, rs_A,
                                y2,  inc_y,
                                z2,  inc_z,
                                a21, rs_A,
                                A22, rs_A, cs_A );

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

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_ofc_var2( integer m_A,
                                     integer m_T,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_T, integer rs_T, integer cs_T )
{
  scomplex* buff_2  = FLA_COMPLEX_PTR( FLA_TWO );
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );

  scomplex  first_elem;
  scomplex  dot_product;
  scomplex  beta, conj_beta;
  scomplex  inv_tau11;
  scomplex  minus_inv_tau11;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  scomplex* buff_y = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_z = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_y  = 1;
  integer       inc_z  = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    scomplex* y0       = buff_y + (0  )*inc_y;
    scomplex* y2       = buff_y + (i+1)*inc_y;

    scomplex* z2       = buff_z + (i+1)*inc_z;

    scomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    scomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    integer       m_ahead  = m_A - i - 1;
    integer       n_ahead  = m_A - i - 1;
    integer       m_behind = i;
    integer       n_behind = i;

    /*------------------------------------------------------------*/

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

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, y2 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, z2 );
      FLA_Fused_Ahx_Ax_opc_var1( m_ahead,
                                 n_ahead,
                                 A22, rs_A, cs_A,
                                 a21, rs_A,
                                 y2,  inc_y,
                                 z2,  inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z2, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bl1_cdot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z2,  inc_z,
                &beta );
      bl1_cinvscals( buff_2, &beta );
      bl1_ccopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y2 );
      // FLA_Scal( inv_tau11, y2 );
      bl1_cscals( &minus_inv_tau11, &conj_beta );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y2, inc_y );
      bl1_cscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y2, inc_y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z2 );
      // FLA_Scal( inv_tau11, z2 );
      bl1_cscals( &minus_inv_tau11, &beta );
      bl1_caxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z2, inc_z );
      bl1_cscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z2, inc_z );

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

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, y0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, y0, a21, A02 );
      bl1_cgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y0,  inc_y );
      bl1_cger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                y0,  inc_y,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, a21, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, a21, A22 );
      FLA_Fused_Gerc2_opc_var1( m_ahead,
                                n_ahead,
                                buff_m1,
                                a21, rs_A,
                                y2,  inc_y,
                                z2,  inc_z,
                                a21, rs_A,
                                A22, rs_A, cs_A );

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

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_ofz_var2( integer m_A,
                                     integer m_T,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_T, integer rs_T, integer cs_T )
{
  dcomplex* buff_2  = FLA_DOUBLE_COMPLEX_PTR( FLA_TWO );
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );

  dcomplex  first_elem;
  dcomplex  dot_product;
  dcomplex  beta, conj_beta;
  dcomplex  inv_tau11;
  dcomplex  minus_inv_tau11;
  integer       i;

  // b_alg = FLA_Obj_length( T );
  integer       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  dcomplex* buff_y = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_z = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  integer       inc_y  = 1;
  integer       inc_z  = 1;

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    dcomplex* y0       = buff_y + (0  )*inc_y;
    dcomplex* y2       = buff_y + (i+1)*inc_y;

    dcomplex* z2       = buff_z + (i+1)*inc_z;

    dcomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    integer       m_ahead  = m_A - i - 1;
    integer       n_ahead  = m_A - i - 1;
    integer       m_behind = i;
    integer       n_behind = i;

    /*------------------------------------------------------------*/

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

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, y2 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, z2 );
      FLA_Fused_Ahx_Ax_opz_var1( m_ahead,
                                 n_ahead,
                                 A22, rs_A, cs_A,
                                 a21, rs_A,
                                 y2,  inc_y,
                                 z2,  inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z2, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bl1_zdot( BLIS1_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z2,  inc_z,
                &beta );
      bl1_zinvscals( buff_2, &beta );
      bl1_zcopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y2 );
      // FLA_Scal( inv_tau11, y2 );
      bl1_zscals( &minus_inv_tau11, &conj_beta );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y2, inc_y );
      bl1_zscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y2, inc_y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z2 );
      // FLA_Scal( inv_tau11, z2 );
      bl1_zscals( &minus_inv_tau11, &beta );
      bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z2, inc_z );
      bl1_zscalv( BLIS1_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z2, inc_z );

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

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, y0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, y0, a21, A02 );
      bl1_zgemv( BLIS1_NO_TRANSPOSE,
                 BLIS1_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y0,  inc_y );
      bl1_zger( BLIS1_NO_CONJUGATE,
                BLIS1_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                y0,  inc_y,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, a21, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, a21, A22 );
      FLA_Fused_Gerc2_opz_var1( m_ahead,
                                n_ahead,
                                buff_m1,
                                a21, rs_A,
                                y2,  inc_y,
                                z2,  inc_z,
                                a21, rs_A,
                                A22, rs_A, cs_A );

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

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  FLA_free( buff_y );
  FLA_free( buff_z );

  return FLA_SUCCESS;
}

