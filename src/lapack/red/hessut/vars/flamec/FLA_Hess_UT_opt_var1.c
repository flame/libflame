/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Hess_UT_opt_var1( FLA_Obj A, FLA_Obj T )
{
  return FLA_Hess_UT_step_opt_var1( A, T );
}

FLA_Error FLA_Hess_UT_step_opt_var1( FLA_Obj A, FLA_Obj T )
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

      FLA_Hess_UT_step_ops_var1( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_Hess_UT_step_opd_var1( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_Hess_UT_step_opc_var1( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_Hess_UT_step_opz_var1( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_ops_var1( int m_A,
                                     int m_T,
                                     float* buff_A, int rs_A, int cs_A, 
                                     float* buff_T, int rs_T, int cs_T )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );

  float     first_elem;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  for ( i = 0; i < b_alg; ++i )
  {
    float*    A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    float*    a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    a21_b    = buff_A + (i  )*cs_A + (i+2)*rs_A;

    float*    A22_t    = buff_A + (i+1)*cs_A + (i+1)*rs_A;
    float*    A22_b    = buff_A + (i+1)*cs_A + (i+2)*rs_A;

    float*    A2_l     = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    A2_r     = buff_A + (i+2)*cs_A + (0  )*rs_A;

    float*    t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    float*    tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       n_behind = i;

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

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21_b, A22_t,
      //                                          A22_b );
      FLA_Apply_H2_UT_l_ops_var1( m_ahead - 1,
                                  n_ahead,
                                  tau11,
                                  a21_b, rs_A,
                                  A22_t, cs_A,
                                  A22_b, rs_A, cs_A );

      // FLA_Apply_H2_UT( FLA_RIGHT, tau11, a21_b, A2_l, A2_r );
      FLA_Apply_H2_UT_r_ops_var1( m_A,
                                  n_ahead - 1,
                                  tau11,
                                  a21_b, rs_A,
                                  A2_l,  rs_A,
                                  A2_r,  rs_A, cs_A );

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

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opd_var1( int m_A,
                                     int m_T,
                                     double* buff_A, int rs_A, int cs_A, 
                                     double* buff_T, int rs_T, int cs_T )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );

  double    first_elem;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  for ( i = 0; i < b_alg; ++i )
  {
    double*   A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    double*   a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   a21_b    = buff_A + (i  )*cs_A + (i+2)*rs_A;

    double*   A22_t    = buff_A + (i+1)*cs_A + (i+1)*rs_A;
    double*   A22_b    = buff_A + (i+1)*cs_A + (i+2)*rs_A;

    double*   A2_l     = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   A2_r     = buff_A + (i+2)*cs_A + (0  )*rs_A;

    double*   t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    double*   tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       n_behind = i;

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

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21_b, A22_t,
      //                                          A22_b );
      FLA_Apply_H2_UT_l_opd_var1( m_ahead - 1,
                                  n_ahead,
                                  tau11,
                                  a21_b, rs_A,
                                  A22_t, cs_A,
                                  A22_b, rs_A, cs_A );

      // FLA_Apply_H2_UT( FLA_RIGHT, tau11, a21_b, A2_l, A2_r );
      FLA_Apply_H2_UT_r_opd_var1( m_A,
                                  n_ahead - 1,
                                  tau11,
                                  a21_b, rs_A,
                                  A2_l,  rs_A,
                                  A2_r,  rs_A, cs_A );

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

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opc_var1( int m_A,
                                     int m_T,
                                     scomplex* buff_A, int rs_A, int cs_A, 
                                     scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );

  scomplex  first_elem;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    scomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* a21_b    = buff_A + (i  )*cs_A + (i+2)*rs_A;

    scomplex* A22_t    = buff_A + (i+1)*cs_A + (i+1)*rs_A;
    scomplex* A22_b    = buff_A + (i+1)*cs_A + (i+2)*rs_A;

    scomplex* A2_l     = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* A2_r     = buff_A + (i+2)*cs_A + (0  )*rs_A;

    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       n_behind = i;

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

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21_b, A22_t,
      //                                          A22_b );
      FLA_Apply_H2_UT_l_opc_var1( m_ahead - 1,
                                  n_ahead,
                                  tau11,
                                  a21_b, rs_A,
                                  A22_t, cs_A,
                                  A22_b, rs_A, cs_A );

      // FLA_Apply_H2_UT( FLA_RIGHT, tau11, a21_b, A2_l, A2_r );
      FLA_Apply_H2_UT_r_opc_var1( m_A,
                                  n_ahead - 1,
                                  tau11,
                                  a21_b, rs_A,
                                  A2_l,  rs_A,
                                  A2_r,  rs_A, cs_A );

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

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opz_var1( int m_A,
                                     int m_T,
                                     dcomplex* buff_A, int rs_A, int cs_A, 
                                     dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );

  dcomplex  first_elem;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    dcomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* a21_b    = buff_A + (i  )*cs_A + (i+2)*rs_A;

    dcomplex* A22_t    = buff_A + (i+1)*cs_A + (i+1)*rs_A;
    dcomplex* A22_b    = buff_A + (i+1)*cs_A + (i+2)*rs_A;

    dcomplex* A2_l     = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* A2_r     = buff_A + (i+2)*cs_A + (0  )*rs_A;

    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       n_behind = i;

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

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Apply_H2_UT( FLA_LEFT, tau11, a21_b, A22_t,
      //                                          A22_b );
      FLA_Apply_H2_UT_l_opz_var1( m_ahead - 1,
                                  n_ahead,
                                  tau11,
                                  a21_b, rs_A,
                                  A22_t, cs_A,
                                  A22_b, rs_A, cs_A );

      // FLA_Apply_H2_UT( FLA_RIGHT, tau11, a21_b, A2_l, A2_r );
      FLA_Apply_H2_UT_r_opz_var1( m_A,
                                  n_ahead - 1,
                                  tau11,
                                  a21_b, rs_A,
                                  A2_l,  rs_A,
                                  A2_r,  rs_A, cs_A );

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

  return FLA_SUCCESS;
}

