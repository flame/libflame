/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Eig_gest_iu_opt_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Datatype datatype;
  integer          m_AB;
  integer          rs_A, cs_A;
  integer          rs_B, cs_B;
  integer          inc_y;
  FLA_Obj      yL, yR;

  datatype = FLA_Obj_datatype( A );

  m_AB     = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );
 
  FLA_Part_1x2( Y,    &yL, &yR,     1, FLA_LEFT );

  inc_y    = FLA_Obj_vector_inc( yL );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_y = FLA_FLOAT_PTR( yL );
      float* buff_B = FLA_FLOAT_PTR( B );

      FLA_Eig_gest_iu_ops_var5( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_y = FLA_DOUBLE_PTR( yL );
      double* buff_B = FLA_DOUBLE_PTR( B );

      FLA_Eig_gest_iu_opd_var5( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_y = FLA_COMPLEX_PTR( yL );
      scomplex* buff_B = FLA_COMPLEX_PTR( B );

      FLA_Eig_gest_iu_opc_var5( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_y = FLA_DOUBLE_COMPLEX_PTR( yL );
      dcomplex* buff_B = FLA_DOUBLE_COMPLEX_PTR( B );

      FLA_Eig_gest_iu_opz_var5( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_ops_var5( integer m_AB,
                                    float* buff_A, integer rs_A, integer cs_A, 
                                    float* buff_y, integer inc_y, 
                                    float* buff_B, integer rs_B, integer cs_B )
{
  float*    buff_m1  = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  float*    buff_m1h = FLA_FLOAT_PTR( FLA_MINUS_ONE_HALF );
  integer       i;

  for ( i = 0; i < m_AB; ++i )
  {
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    float*    b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;
    float*    B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    float     psi11;

    integer       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bl1_sinvscals( beta11, alpha11 );
    bl1_sinvscals( beta11, alpha11 );

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_MINUS_ONE_HALF, psi11 );
    bl1_smult3( buff_m1h, alpha11, &psi11 );

    // FLA_Inv_scal_external( beta11, a12t );
    bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, a12t, b12t, A22 );
    bl1_sher2( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_ahead,
               buff_m1,
               a12t, cs_A,
               b12t, cs_B,
               A22,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a12t );
    bl1_strsv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_ahead,
               B22,  rs_B, cs_B,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opd_var5( integer m_AB,
                                    double* buff_A, integer rs_A, integer cs_A, 
                                    double* buff_y, integer inc_y, 
                                    double* buff_B, integer rs_B, integer cs_B )
{
  double*   buff_m1  = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  double*   buff_m1h = FLA_DOUBLE_PTR( FLA_MINUS_ONE_HALF );
  integer       i;

  for ( i = 0; i < m_AB; ++i )
  {
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    double*   b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;
    double*   B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    double    psi11;

    integer       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bl1_dinvscals( beta11, alpha11 );
    bl1_dinvscals( beta11, alpha11 );

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_MINUS_ONE_HALF, psi11 );
    bl1_dmult3( buff_m1h, alpha11, &psi11 );

    // FLA_Inv_scal_external( beta11, a12t );
    bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, a12t, b12t, A22 );
    bl1_dher2( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_ahead,
               buff_m1,
               a12t, cs_A,
               b12t, cs_B,
               A22,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a12t );
    bl1_dtrsv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_ahead,
               B22,  rs_B, cs_B,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opc_var5( integer m_AB,
                                    scomplex* buff_A, integer rs_A, integer cs_A, 
                                    scomplex* buff_y, integer inc_y, 
                                    scomplex* buff_B, integer rs_B, integer cs_B )
{
  scomplex* buff_m1  = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  scomplex* buff_m1h = FLA_COMPLEX_PTR( FLA_MINUS_ONE_HALF );
  integer       i;

  for ( i = 0; i < m_AB; ++i )
  {
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    scomplex* b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;
    scomplex* B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    scomplex  psi11;

    integer       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bl1_cinvscals( beta11, alpha11 );
    bl1_cinvscals( beta11, alpha11 );

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_MINUS_ONE_HALF, psi11 );
    bl1_cmult3( buff_m1h, alpha11, &psi11 );

    // FLA_Inv_scal_external( beta11, a12t );
    bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, a12t, b12t, A22 );
    bl1_cher2( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_ahead,
               buff_m1,
               a12t, cs_A,
               b12t, cs_B,
               A22,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a12t );
    bl1_ctrsv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_ahead,
               B22,  rs_B, cs_B,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opz_var5( integer m_AB,
                                    dcomplex* buff_A, integer rs_A, integer cs_A, 
                                    dcomplex* buff_y, integer inc_y, 
                                    dcomplex* buff_B, integer rs_B, integer cs_B )
{
  dcomplex* buff_m1  = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  dcomplex* buff_m1h = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE_HALF );
  integer       i;

  for ( i = 0; i < m_AB; ++i )
  {
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    dcomplex* b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;
    dcomplex* B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    dcomplex  psi11;

    integer       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bl1_zinvscals( beta11, alpha11 );
    bl1_zinvscals( beta11, alpha11 );

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_MINUS_ONE_HALF, psi11 );
    bl1_zmult3( buff_m1h, alpha11, &psi11 );

    // FLA_Inv_scal_external( beta11, a12t );
    bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, a12t, b12t, A22 );
    bl1_zher2( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_ahead,
               buff_m1,
               a12t, cs_A,
               b12t, cs_B,
               A22,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a12t );
    bl1_ztrsv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_ahead,
               B22,  rs_B, cs_B,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

