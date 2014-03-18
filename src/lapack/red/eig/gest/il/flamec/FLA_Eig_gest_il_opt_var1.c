
#include "FLAME.h"

FLA_Error FLA_Eig_gest_il_opt_var1( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Datatype datatype;
  int          m_AB;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  int          inc_y;
  FLA_Obj      yT, yB;

  datatype = FLA_Obj_datatype( A );

  m_AB     = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );
 
  FLA_Part_2x1( Y,    &yT,
                      &yB,     1, FLA_TOP );

  inc_y    = FLA_Obj_vector_inc( yT );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_y = FLA_FLOAT_PTR( yT );
      float* buff_B = FLA_FLOAT_PTR( B );

      FLA_Eig_gest_il_ops_var1( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_y = FLA_DOUBLE_PTR( yT );
      double* buff_B = FLA_DOUBLE_PTR( B );

      FLA_Eig_gest_il_opd_var1( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_y = FLA_COMPLEX_PTR( yT );
      scomplex* buff_B = FLA_COMPLEX_PTR( B );

      FLA_Eig_gest_il_opc_var1( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_y = FLA_DOUBLE_COMPLEX_PTR( yT );
      dcomplex* buff_B = FLA_DOUBLE_COMPLEX_PTR( B );

      FLA_Eig_gest_il_opz_var1( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_il_ops_var1( int m_AB,
                                    float* buff_A, int rs_A, int cs_A, 
                                    float* buff_y, int inc_y, 
                                    float* buff_B, int rs_B, int cs_B )
{
  float*    buff_1   = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0   = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1  = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  float*    buff_m1h = FLA_FLOAT_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    float*    A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;

    float*    y10t     = buff_y + (0  )*inc_y;

    float*    B00      = buff_B + (0  )*cs_B + (0  )*rs_B;
    float*    b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    float*    beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Hemvc_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, A00, b10t, FLA_ZERO, y10t_t );
    bl1_shemv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_behind,
               buff_1,
               A00,  rs_A, cs_A,
               b10t, cs_B,
               buff_0,
               y10t, inc_y );

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B00, a10t );
    bl1_strsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_behind,
               B00,  rs_B, cs_B,
               a10t, cs_A );

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t_t, a10t );
    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, inc_y,
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

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t_t, a10t );
    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, inc_y,
                a10t, cs_A );

    // FLA_Inv_scal_external( beta11, a10t );
    bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a10t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_il_opd_var1( int m_AB,
                                    double* buff_A, int rs_A, int cs_A, 
                                    double* buff_y, int inc_y, 
                                    double* buff_B, int rs_B, int cs_B )
{
  double*   buff_1   = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0   = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1  = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  double*   buff_m1h = FLA_DOUBLE_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    double*   A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;

    double*   y10t     = buff_y + (0  )*inc_y;

    double*   B00      = buff_B + (0  )*cs_B + (0  )*rs_B;
    double*   b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    double*   beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Hemvc_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, A00, b10t, FLA_ZERO, y10t_t );
    bl1_dhemv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_behind,
               buff_1,
               A00,  rs_A, cs_A,
               b10t, cs_B,
               buff_0,
               y10t, inc_y );

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B00, a10t );
    bl1_dtrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_behind,
               B00,  rs_B, cs_B,
               a10t, cs_A );

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t_t, a10t );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, inc_y,
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

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t_t, a10t );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, inc_y,
                a10t, cs_A );

    // FLA_Inv_scal_external( beta11, a10t );
    bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a10t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_il_opc_var1( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B )
{
  scomplex* buff_1   = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0   = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1  = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  scomplex* buff_m1h = FLA_COMPLEX_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    scomplex* A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;

    scomplex* y10t     = buff_y + (0  )*inc_y;

    scomplex* B00      = buff_B + (0  )*cs_B + (0  )*rs_B;
    scomplex* b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    scomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Hemvc_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, A00, b10t, FLA_ZERO, y10t_t );
    bl1_chemv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_behind,
               buff_1,
               A00,  rs_A, cs_A,
               b10t, cs_B,
               buff_0,
               y10t, inc_y );

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B00, a10t );
    bl1_ctrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_behind,
               B00,  rs_B, cs_B,
               a10t, cs_A );

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t_t, a10t );
    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, inc_y,
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

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t_t, a10t );
    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, inc_y,
                a10t, cs_A );

    // FLA_Inv_scal_external( beta11, a10t );
    bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a10t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_il_opz_var1( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B )
{
  dcomplex* buff_1   = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0   = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1  = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  dcomplex* buff_m1h = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    dcomplex* A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;

    dcomplex* y10t     = buff_y + (0  )*inc_y;

    dcomplex* B00      = buff_B + (0  )*cs_B + (0  )*rs_B;
    dcomplex* b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    dcomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Hemvc_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, A00, b10t, FLA_ZERO, y10t_t );
    bl1_zhemv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_behind,
               buff_1,
               A00,  rs_A, cs_A,
               b10t, cs_B,
               buff_0,
               y10t, inc_y );

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B00, a10t );
    bl1_ztrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJ_NO_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_behind,
               B00,  rs_B, cs_B,
               a10t, cs_A );

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t_t, a10t );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, inc_y,
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

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t_t, a10t );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y10t, inc_y,
                a10t, cs_A );

    // FLA_Inv_scal_external( beta11, a10t );
    bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a10t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

