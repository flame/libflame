/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Eig_gest_il_opt_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Datatype datatype;
  integer          m_AB;
  integer          rs_A, cs_A;
  integer          rs_Y, cs_Y;
  integer          rs_B, cs_B;

  datatype = FLA_Obj_datatype( A );

  m_AB     = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_Y     = FLA_Obj_row_stride( Y );
  cs_Y     = FLA_Obj_col_stride( Y );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );
 
  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_Y = FLA_FLOAT_PTR( Y );
      float* buff_B = FLA_FLOAT_PTR( B );

      FLA_Eig_gest_il_ops_var3( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_Y, rs_Y, cs_Y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_Y = FLA_DOUBLE_PTR( Y );
      double* buff_B = FLA_DOUBLE_PTR( B );

      FLA_Eig_gest_il_opd_var3( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_Y, rs_Y, cs_Y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_Y = FLA_COMPLEX_PTR( Y );
      scomplex* buff_B = FLA_COMPLEX_PTR( B );

      FLA_Eig_gest_il_opc_var3( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_Y, rs_Y, cs_Y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_Y = FLA_DOUBLE_COMPLEX_PTR( Y );
      dcomplex* buff_B = FLA_DOUBLE_COMPLEX_PTR( B );

      FLA_Eig_gest_il_opz_var3( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_Y, rs_Y, cs_Y,
                                buff_B, rs_B, cs_B );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_il_ops_var3( integer m_AB,
                                    float* buff_A, integer rs_A, integer cs_A, 
                                    float* buff_Y, integer rs_Y, integer cs_Y,
                                    float* buff_B, integer rs_B, integer cs_B )
{
  float*    buff_1   = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1  = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  float*    buff_m1h = FLA_FLOAT_PTR( FLA_MINUS_ONE_HALF );
  integer       i;

  for ( i = 0; i < m_AB; ++i )
  {
    float*    a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    float*    y10t     = buff_Y + (0  )*cs_Y + (i  )*rs_Y;
    float*    Y20      = buff_Y + (0  )*cs_Y + (i+1)*rs_Y;
    float*    y21      = buff_Y + (i  )*cs_Y + (i+1)*rs_Y;

    float*    b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    float*    B20      = buff_B + (0  )*cs_B + (i+1)*rs_B;
    float*    beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    float*    b21      = buff_B + (i  )*cs_B + (i+1)*rs_B;

    integer       m_ahead  = m_AB - i - 1;
    integer       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t, a10t );
    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, cs_Y,
                a10t, cs_A );

    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a10t, b10t, FLA_ONE, alpha11 );
    bl1_sdot2s( BLIS1_CONJUGATE,
                m_behind,
                buff_m1,
                a10t, cs_A,
                b10t, cs_B,
                buff_1,
                alpha11 );

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bl1_sinvscals( beta11, alpha11 );
    bl1_sinvscals( beta11, alpha11 );

    // FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, A20, b10t, FLA_ONE, a21 ); 
    bl1_sgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead,
               m_behind,
               buff_m1,
               A20,  rs_A, cs_A,
               b10t, cs_B,
               buff_1,
               a21,  rs_A );

    // FLA_Inv_scal_external( beta11, a21 );
    bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a21, rs_A );
  
    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t, a10t );
    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, cs_Y,
                a10t, cs_A );

    // FLA_Inv_scal_external( beta11, a10t );
    bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a10t, cs_A );

    // FLA_Ger_external( FLA_ONE, b21, a10t, Y20 );
    bl1_sger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_ahead,
              m_behind,
              buff_1,
              b21,  rs_B,
              a10t, cs_A,
              Y20,  rs_Y, cs_Y );

    // FLA_Copy_external( b21, y21 );
    // FLA_Scal_external( alpha11, y21 );
    bl1_scopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                b21, rs_B,
                y21, rs_Y );
    bl1_sscalv( BLIS1_NO_CONJUGATE,
                m_ahead,
                alpha11,
                y21, rs_Y );

    // FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_ONE, B20, a10t, FLA_ONE, y21 );
    bl1_sgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead,
               m_behind,
               buff_1,
               B20,  rs_B, cs_B,
               a10t, cs_A,
               buff_1,
               y21,  rs_Y );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_il_opd_var3( integer m_AB,
                                    double* buff_A, integer rs_A, integer cs_A, 
                                    double* buff_Y, integer rs_Y, integer cs_Y,
                                    double* buff_B, integer rs_B, integer cs_B )
{
  double*   buff_1   = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1  = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  double*   buff_m1h = FLA_DOUBLE_PTR( FLA_MINUS_ONE_HALF );
  integer       i;

  for ( i = 0; i < m_AB; ++i )
  {
    double*   a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    double*   y10t     = buff_Y + (0  )*cs_Y + (i  )*rs_Y;
    double*   Y20      = buff_Y + (0  )*cs_Y + (i+1)*rs_Y;
    double*   y21      = buff_Y + (i  )*cs_Y + (i+1)*rs_Y;

    double*   b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    double*   B20      = buff_B + (0  )*cs_B + (i+1)*rs_B;
    double*   beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    double*   b21      = buff_B + (i  )*cs_B + (i+1)*rs_B;

    integer       m_ahead  = m_AB - i - 1;
    integer       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t, a10t );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, cs_Y,
                a10t, cs_A );

    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a10t, b10t, FLA_ONE, alpha11 );
    bl1_ddot2s( BLIS1_CONJUGATE,
                m_behind,
                buff_m1,
                a10t, cs_A,
                b10t, cs_B,
                buff_1,
                alpha11 );

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bl1_dinvscals( beta11, alpha11 );
    bl1_dinvscals( beta11, alpha11 );

    // FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, A20, b10t, FLA_ONE, a21 ); 
    bl1_dgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead,
               m_behind,
               buff_m1,
               A20,  rs_A, cs_A,
               b10t, cs_B,
               buff_1,
               a21,  rs_A );

    // FLA_Inv_scal_external( beta11, a21 );
    bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a21, rs_A );
  
    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t, a10t );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, cs_Y,
                a10t, cs_A );

    // FLA_Inv_scal_external( beta11, a10t );
    bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a10t, cs_A );

    // FLA_Ger_external( FLA_ONE, b21, a10t, Y20 );
    bl1_dger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_ahead,
              m_behind,
              buff_1,
              b21,  rs_B,
              a10t, cs_A,
              Y20,  rs_Y, cs_Y );

    // FLA_Copy_external( b21, y21 );
    // FLA_Scal_external( alpha11, y21 );
    bl1_dcopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                b21, rs_B,
                y21, rs_Y );
    bl1_dscalv( BLIS1_NO_CONJUGATE,
                m_ahead,
                alpha11,
                y21, rs_Y );

    // FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_ONE, B20, a10t, FLA_ONE, y21 );
    bl1_dgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead,
               m_behind,
               buff_1,
               B20,  rs_B, cs_B,
               a10t, cs_A,
               buff_1,
               y21,  rs_Y );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_il_opc_var3( integer m_AB,
                                    scomplex* buff_A, integer rs_A, integer cs_A, 
                                    scomplex* buff_Y, integer rs_Y, integer cs_Y,
                                    scomplex* buff_B, integer rs_B, integer cs_B )
{
  scomplex* buff_1   = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1  = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  scomplex* buff_m1h = FLA_COMPLEX_PTR( FLA_MINUS_ONE_HALF );
  integer       i;

  for ( i = 0; i < m_AB; ++i )
  {
    scomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    scomplex* y10t     = buff_Y + (0  )*cs_Y + (i  )*rs_Y;
    scomplex* Y20      = buff_Y + (0  )*cs_Y + (i+1)*rs_Y;
    scomplex* y21      = buff_Y + (i  )*cs_Y + (i+1)*rs_Y;

    scomplex* b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    scomplex* B20      = buff_B + (0  )*cs_B + (i+1)*rs_B;
    scomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    scomplex* b21      = buff_B + (i  )*cs_B + (i+1)*rs_B;

    integer       m_ahead  = m_AB - i - 1;
    integer       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t, a10t );
    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, cs_Y,
                a10t, cs_A );

    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a10t, b10t, FLA_ONE, alpha11 );
    bl1_cdot2s( BLIS1_CONJUGATE,
                m_behind,
                buff_m1,
                a10t, cs_A,
                b10t, cs_B,
                buff_1,
                alpha11 );

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bl1_cinvscals( beta11, alpha11 );
    bl1_cinvscals( beta11, alpha11 );

    // FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, A20, b10t, FLA_ONE, a21 ); 
    bl1_cgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead,
               m_behind,
               buff_m1,
               A20,  rs_A, cs_A,
               b10t, cs_B,
               buff_1,
               a21,  rs_A );

    // FLA_Inv_scal_external( beta11, a21 );
    bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a21, rs_A );
  
    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t, a10t );
    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, cs_Y,
                a10t, cs_A );

    // FLA_Inv_scal_external( beta11, a10t );
    bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a10t, cs_A );

    // FLA_Ger_external( FLA_ONE, b21, a10t, Y20 );
    bl1_cger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_ahead,
              m_behind,
              buff_1,
              b21,  rs_B,
              a10t, cs_A,
              Y20,  rs_Y, cs_Y );

    // FLA_Copy_external( b21, y21 );
    // FLA_Scal_external( alpha11, y21 );
    bl1_ccopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                b21, rs_B,
                y21, rs_Y );
    bl1_cscalv( BLIS1_NO_CONJUGATE,
                m_ahead,
                alpha11,
                y21, rs_Y );

    // FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_ONE, B20, a10t, FLA_ONE, y21 );
    bl1_cgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead,
               m_behind,
               buff_1,
               B20,  rs_B, cs_B,
               a10t, cs_A,
               buff_1,
               y21,  rs_Y );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_il_opz_var3( integer m_AB,
                                    dcomplex* buff_A, integer rs_A, integer cs_A, 
                                    dcomplex* buff_Y, integer rs_Y, integer cs_Y,
                                    dcomplex* buff_B, integer rs_B, integer cs_B )
{
  dcomplex* buff_1   = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1  = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  dcomplex* buff_m1h = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE_HALF );
  integer       i;

  for ( i = 0; i < m_AB; ++i )
  {
    dcomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    dcomplex* y10t     = buff_Y + (0  )*cs_Y + (i  )*rs_Y;
    dcomplex* Y20      = buff_Y + (0  )*cs_Y + (i+1)*rs_Y;
    dcomplex* y21      = buff_Y + (i  )*cs_Y + (i+1)*rs_Y;

    dcomplex* b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    dcomplex* B20      = buff_B + (0  )*cs_B + (i+1)*rs_B;
    dcomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    dcomplex* b21      = buff_B + (i  )*cs_B + (i+1)*rs_B;

    integer       m_ahead  = m_AB - i - 1;
    integer       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t, a10t );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, cs_Y,
                a10t, cs_A );

    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a10t, b10t, FLA_ONE, alpha11 );
    bl1_zdot2s( BLIS1_CONJUGATE,
                m_behind,
                buff_m1,
                a10t, cs_A,
                b10t, cs_B,
                buff_1,
                alpha11 );

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bl1_zinvscals( beta11, alpha11 );
    bl1_zinvscals( beta11, alpha11 );

    // FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, A20, b10t, FLA_ONE, a21 ); 
    bl1_zgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead,
               m_behind,
               buff_m1,
               A20,  rs_A, cs_A,
               b10t, cs_B,
               buff_1,
               a21,  rs_A );

    // FLA_Inv_scal_external( beta11, a21 );
    bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a21, rs_A );
  
    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t, a10t );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, cs_Y,
                a10t, cs_A );

    // FLA_Inv_scal_external( beta11, a10t );
    bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a10t, cs_A );

    // FLA_Ger_external( FLA_ONE, b21, a10t, Y20 );
    bl1_zger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_ahead,
              m_behind,
              buff_1,
              b21,  rs_B,
              a10t, cs_A,
              Y20,  rs_Y, cs_Y );

    // FLA_Copy_external( b21, y21 );
    // FLA_Scal_external( alpha11, y21 );
    bl1_zcopyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                b21, rs_B,
                y21, rs_Y );
    bl1_zscalv( BLIS1_NO_CONJUGATE,
                m_ahead,
                alpha11,
                y21, rs_Y );

    // FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_ONE, B20, a10t, FLA_ONE, y21 );
    bl1_zgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_CONJUGATE,
               m_ahead,
               m_behind,
               buff_1,
               B20,  rs_B, cs_B,
               a10t, cs_A,
               buff_1,
               y21,  rs_Y );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

