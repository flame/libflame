/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_QR2_UT_opt_var1( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T )
{
  FLA_Datatype datatype;
  integer          m_UT, m_D;
  integer          rs_U, cs_U;
  integer          rs_D, cs_D;
  integer          rs_T, cs_T;

  datatype = FLA_Obj_datatype( U );

  m_UT     = FLA_Obj_length( U );
  m_D      = FLA_Obj_length( D );

  rs_U     = FLA_Obj_row_stride( U );
  cs_U     = FLA_Obj_col_stride( U );
  rs_D     = FLA_Obj_row_stride( D );
  cs_D     = FLA_Obj_col_stride( D );
  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_U = FLA_FLOAT_PTR( U );
      float* buff_D = FLA_FLOAT_PTR( D );
      float* buff_T = FLA_FLOAT_PTR( T );

      FLA_QR2_UT_ops_var1( m_UT,
                           m_D,
                           buff_U, rs_U, cs_U,
                           buff_D, rs_D, cs_D,
                           buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_U = FLA_DOUBLE_PTR( U );
      double* buff_D = FLA_DOUBLE_PTR( D );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_QR2_UT_opd_var1( m_UT,
                           m_D,
                           buff_U, rs_U, cs_U,
                           buff_D, rs_D, cs_D,
                           buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_U = FLA_COMPLEX_PTR( U );
      scomplex* buff_D = FLA_COMPLEX_PTR( D );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_QR2_UT_opc_var1( m_UT,
                           m_D,
                           buff_U, rs_U, cs_U,
                           buff_D, rs_D, cs_D,
                           buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_U = FLA_DOUBLE_COMPLEX_PTR( U );
      dcomplex* buff_D = FLA_DOUBLE_COMPLEX_PTR( D );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_QR2_UT_opz_var1( m_UT,
                           m_D,
                           buff_U, rs_U, cs_U,
                           buff_D, rs_D, cs_D,
                           buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR2_UT_ops_var1( integer m_UT,
                               integer m_D,
                               float* buff_U, integer rs_U, integer cs_U,
                               float* buff_D, integer rs_D, integer cs_D,
                               float* buff_T, integer rs_T, integer cs_T )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  integer       i;

  for ( i = 0; i < m_UT; ++i )
  {
    float*    upsilon11 = buff_U + (i  )*cs_U + (i  )*rs_U;
    float*    u12t      = buff_U + (i+1)*cs_U + (i  )*rs_U;

    float*    D0        = buff_D + (0  )*cs_D + (0  )*rs_D;
    float*    d1        = buff_D + (i  )*cs_D + (0  )*rs_D;
    float*    D2        = buff_D + (i+1)*cs_D + (0  )*rs_D;

    float*    tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    float*    t01       = buff_T + (i  )*cs_T + (0  )*rs_T;

    integer       mn_ahead  = m_UT - i - 1;
    integer       mn_behind = i;

    //------------------------------------------------------------//

    // FLA_Househ2_UT( FLA_LEFT,
    //                 upsilon11,
    //                 d1, tau11 );
    FLA_Househ2_UT_l_ops( m_D,
                          upsilon11,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, d1, u12t,
    //                                       D2 );
    FLA_Apply_H2_UT_l_ops_var1( m_D,
                                mn_ahead,
                                tau11,
                                d1, rs_D,
                                u12t, cs_U,
                                D2, rs_D, cs_D );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, D0, d1, FLA_ZERO, t01 );
    bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_D,
               mn_behind,
               buff_1,
               D0, rs_D, cs_D,
               d1, rs_D,
               buff_0,
               t01, rs_T );

    //------------------------------------------------------------//

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR2_UT_opd_var1( integer m_UT,
                               integer m_D,
                               double* buff_U, integer rs_U, integer cs_U,
                               double* buff_D, integer rs_D, integer cs_D,
                               double* buff_T, integer rs_T, integer cs_T )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  integer       i;

  for ( i = 0; i < m_UT; ++i )
  {
    double*   upsilon11 = buff_U + (i  )*cs_U + (i  )*rs_U;
    double*   u12t      = buff_U + (i+1)*cs_U + (i  )*rs_U;

    double*   D0        = buff_D + (0  )*cs_D + (0  )*rs_D;
    double*   d1        = buff_D + (i  )*cs_D + (0  )*rs_D;
    double*   D2        = buff_D + (i+1)*cs_D + (0  )*rs_D;

    double*   tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    double*   t01       = buff_T + (i  )*cs_T + (0  )*rs_T;

    integer       mn_ahead  = m_UT - i - 1;
    integer       mn_behind = i;

    //------------------------------------------------------------//

    // FLA_Househ2_UT( FLA_LEFT,
    //                 upsilon11,
    //                 d1, tau11 );
    FLA_Househ2_UT_l_opd( m_D,
                          upsilon11,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, d1, u12t,
    //                                       D2 );
    FLA_Apply_H2_UT_l_opd_var1( m_D,
                                mn_ahead,
                                tau11,
                                d1, rs_D,
                                u12t, cs_U,
                                D2, rs_D, cs_D );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, D0, d1, FLA_ZERO, t01 );
    bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_D,
               mn_behind,
               buff_1,
               D0, rs_D, cs_D,
               d1, rs_D,
               buff_0,
               t01, rs_T );

    //------------------------------------------------------------//

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR2_UT_opc_var1( integer m_UT,
                               integer m_D,
                               scomplex* buff_U, integer rs_U, integer cs_U,
                               scomplex* buff_D, integer rs_D, integer cs_D,
                               scomplex* buff_T, integer rs_T, integer cs_T )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  integer       i;

  for ( i = 0; i < m_UT; ++i )
  {
    scomplex* upsilon11 = buff_U + (i  )*cs_U + (i  )*rs_U;
    scomplex* u12t      = buff_U + (i+1)*cs_U + (i  )*rs_U;

    scomplex* D0        = buff_D + (0  )*cs_D + (0  )*rs_D;
    scomplex* d1        = buff_D + (i  )*cs_D + (0  )*rs_D;
    scomplex* D2        = buff_D + (i+1)*cs_D + (0  )*rs_D;

    scomplex* tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    scomplex* t01       = buff_T + (i  )*cs_T + (0  )*rs_T;

    integer       mn_ahead  = m_UT - i - 1;
    integer       mn_behind = i;

    //------------------------------------------------------------//

    // FLA_Househ2_UT( FLA_LEFT,
    //                 upsilon11,
    //                 d1, tau11 );
    FLA_Househ2_UT_l_opc( m_D,
                          upsilon11,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, d1, u12t,
    //                                       D2 );
    FLA_Apply_H2_UT_l_opc_var1( m_D,
                                mn_ahead,
                                tau11,
                                d1, rs_D,
                                u12t, cs_U,
                                D2, rs_D, cs_D );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, D0, d1, FLA_ZERO, t01 );
    bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_D,
               mn_behind,
               buff_1,
               D0, rs_D, cs_D,
               d1, rs_D,
               buff_0,
               t01, rs_T );

    //------------------------------------------------------------//

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR2_UT_opz_var1( integer m_UT,
                               integer m_D,
                               dcomplex* buff_U, integer rs_U, integer cs_U,
                               dcomplex* buff_D, integer rs_D, integer cs_D,
                               dcomplex* buff_T, integer rs_T, integer cs_T )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  integer       i;

  for ( i = 0; i < m_UT; ++i )
  {
    dcomplex* upsilon11 = buff_U + (i  )*cs_U + (i  )*rs_U;
    dcomplex* u12t      = buff_U + (i+1)*cs_U + (i  )*rs_U;

    dcomplex* D0        = buff_D + (0  )*cs_D + (0  )*rs_D;
    dcomplex* d1        = buff_D + (i  )*cs_D + (0  )*rs_D;
    dcomplex* D2        = buff_D + (i+1)*cs_D + (0  )*rs_D;

    dcomplex* tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    dcomplex* t01       = buff_T + (i  )*cs_T + (0  )*rs_T;

    integer       mn_ahead  = m_UT - i - 1;
    integer       mn_behind = i;

    //------------------------------------------------------------//

    // FLA_Househ2_UT( FLA_LEFT,
    //                 upsilon11,
    //                 d1, tau11 );
    FLA_Househ2_UT_l_opz( m_D,
                          upsilon11,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, d1, u12t,
    //                                       D2 );
    FLA_Apply_H2_UT_l_opz_var1( m_D,
                                mn_ahead,
                                tau11,
                                d1, rs_D,
                                u12t, cs_U,
                                D2, rs_D, cs_D );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, D0, d1, FLA_ZERO, t01 );
    bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_D,
               mn_behind,
               buff_1,
               D0, rs_D, cs_D,
               d1, rs_D,
               buff_0,
               t01, rs_T );

    //------------------------------------------------------------//

  }

  return FLA_SUCCESS;
}

