/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Eig_gest_nu_opt_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Datatype datatype;
  int          m_AB;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  int          inc_y;
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

      FLA_Eig_gest_nu_ops_var4( m_AB,
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

      FLA_Eig_gest_nu_opd_var4( m_AB,
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

      FLA_Eig_gest_nu_opc_var4( m_AB,
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

      FLA_Eig_gest_nu_opz_var4( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nu_ops_var4( int m_AB,
                                    float* buff_A, int rs_A, int cs_A, 
                                    float* buff_y, int inc_y, 
                                    float* buff_B, int rs_B, int cs_B )
{
  float*    buff_1   = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_1h  = FLA_FLOAT_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    float*    A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    float*    b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    float*    beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    float     psi11;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bl1_smult3( buff_1h, alpha11, &psi11 );

    // FLA_Axpy_external( psi11, b01, a01 );
    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b01, rs_B,
                a01, rs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_NO_CONJUGATE,
    //                     FLA_ONE, a01, b01, A00 );
    bl1_sher2( BLIS1_UPPER_TRIANGULAR,
               BLIS1_NO_CONJUGATE,
               m_behind,
               buff_1,
               a01, rs_A,
               b01, rs_B,
               A00, rs_A, cs_A );

    // FLA_Axpy_external( psi11, b01, a01 );
    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b01, rs_B,
                a01, rs_A );

    // FLA_Scal_external( beta11, a01 );
    bl1_sscalv( BLIS1_NO_CONJUGATE,
                m_behind,
                beta11,
                a01, rs_A );

    // FLA_Scal_external( beta11, alpha11 );
    bl1_sscals( beta11, alpha11 );
    bl1_sscals( beta11, alpha11 );

    // FLA_Ger_external( FLA_ONE, b01, a12t, A02 );
    bl1_sger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_behind,
              m_ahead,
              buff_1,
              b01,  rs_B,
              a12t, cs_A,
              A02,  rs_A, cs_A );

    // FLA_Scal_external( beta11, a12t );
    bl1_sscalv( BLIS1_NO_CONJUGATE,
                m_ahead,
                beta11,
                a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nu_opd_var4( int m_AB,
                                    double* buff_A, int rs_A, int cs_A, 
                                    double* buff_y, int inc_y, 
                                    double* buff_B, int rs_B, int cs_B )
{
  double*   buff_1   = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_1h  = FLA_DOUBLE_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    double*   A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    double*   b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    double*   beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    double    psi11;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bl1_dmult3( buff_1h, alpha11, &psi11 );

    // FLA_Axpy_external( psi11, b01, a01 );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b01, rs_B,
                a01, rs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_NO_CONJUGATE,
    //                     FLA_ONE, a01, b01, A00 );
    bl1_dher2( BLIS1_UPPER_TRIANGULAR,
               BLIS1_NO_CONJUGATE,
               m_behind,
               buff_1,
               a01, rs_A,
               b01, rs_B,
               A00, rs_A, cs_A );

    // FLA_Axpy_external( psi11, b01, a01 );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b01, rs_B,
                a01, rs_A );

    // FLA_Scal_external( beta11, a01 );
    bl1_dscalv( BLIS1_NO_CONJUGATE,
                m_behind,
                beta11,
                a01, rs_A );

    // FLA_Scal_external( beta11, alpha11 );
    bl1_dscals( beta11, alpha11 );
    bl1_dscals( beta11, alpha11 );

    // FLA_Ger_external( FLA_ONE, b01, a12t, A02 );
    bl1_dger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_behind,
              m_ahead,
              buff_1,
              b01,  rs_B,
              a12t, cs_A,
              A02,  rs_A, cs_A );

    // FLA_Scal_external( beta11, a12t );
    bl1_dscalv( BLIS1_NO_CONJUGATE,
                m_ahead,
                beta11,
                a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nu_opc_var4( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B )
{
  scomplex* buff_1   = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_1h  = FLA_COMPLEX_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    scomplex* A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    scomplex* b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    scomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    scomplex  psi11;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bl1_cmult3( buff_1h, alpha11, &psi11 );

    // FLA_Axpy_external( psi11, b01, a01 );
    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b01, rs_B,
                a01, rs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_NO_CONJUGATE,
    //                     FLA_ONE, a01, b01, A00 );
    bl1_cher2( BLIS1_UPPER_TRIANGULAR,
               BLIS1_NO_CONJUGATE,
               m_behind,
               buff_1,
               a01, rs_A,
               b01, rs_B,
               A00, rs_A, cs_A );

    // FLA_Axpy_external( psi11, b01, a01 );
    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b01, rs_B,
                a01, rs_A );

    // FLA_Scal_external( beta11, a01 );
    bl1_cscalv( BLIS1_NO_CONJUGATE,
                m_behind,
                beta11,
                a01, rs_A );

    // FLA_Scal_external( beta11, alpha11 );
    bl1_cscals( beta11, alpha11 );
    bl1_cscals( beta11, alpha11 );

    // FLA_Ger_external( FLA_ONE, b01, a12t, A02 );
    bl1_cger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_behind,
              m_ahead,
              buff_1,
              b01,  rs_B,
              a12t, cs_A,
              A02,  rs_A, cs_A );

    // FLA_Scal_external( beta11, a12t );
    bl1_cscalv( BLIS1_NO_CONJUGATE,
                m_ahead,
                beta11,
                a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nu_opz_var4( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B )
{
  dcomplex* buff_1   = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_1h  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    dcomplex* A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    dcomplex* b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    dcomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    dcomplex  psi11;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bl1_zmult3( buff_1h, alpha11, &psi11 );

    // FLA_Axpy_external( psi11, b01, a01 );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b01, rs_B,
                a01, rs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_NO_CONJUGATE,
    //                     FLA_ONE, a01, b01, A00 );
    bl1_zher2( BLIS1_UPPER_TRIANGULAR,
               BLIS1_NO_CONJUGATE,
               m_behind,
               buff_1,
               a01, rs_A,
               b01, rs_B,
               A00, rs_A, cs_A );

    // FLA_Axpy_external( psi11, b01, a01 );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b01, rs_B,
                a01, rs_A );

    // FLA_Scal_external( beta11, a01 );
    bl1_zscalv( BLIS1_NO_CONJUGATE,
                m_behind,
                beta11,
                a01, rs_A );

    // FLA_Scal_external( beta11, alpha11 );
    bl1_zscals( beta11, alpha11 );
    bl1_zscals( beta11, alpha11 );

    // FLA_Ger_external( FLA_ONE, b01, a12t, A02 );
    bl1_zger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_behind,
              m_ahead,
              buff_1,
              b01,  rs_B,
              a12t, cs_A,
              A02,  rs_A, cs_A );

    // FLA_Scal_external( beta11, a12t );
    bl1_zscalv( BLIS1_NO_CONJUGATE,
                m_ahead,
                beta11,
                a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
