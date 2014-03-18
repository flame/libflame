
#include "FLAME.h"

FLA_Error FLA_Eig_gest_nl_opt_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
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

      FLA_Eig_gest_nl_ops_var5( m_AB,
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

      FLA_Eig_gest_nl_opd_var5( m_AB,
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

      FLA_Eig_gest_nl_opc_var5( m_AB,
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

      FLA_Eig_gest_nl_opz_var5( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_ops_var5( int m_AB,
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
    float*    a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;

    float*    B00      = buff_B + (0  )*cs_B + (0  )*rs_B;
    float*    b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    float*    beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    float     psi11;

    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bl1_smult3( buff_1h, alpha11, &psi11 );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B00, a10t );
    bl1_strmv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_behind,
               B00,  rs_B, cs_B,
               a10t, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Her2c_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, a10t, b10t, A00 );
    bl1_sher2( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_behind,
               buff_1,
               a10t, cs_A,
               b10t, cs_B,
               A00,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Scal_external( beta11, a10t );
    bl1_sscalv( BLIS1_NO_CONJUGATE,
                m_behind,
                beta11,
                a10t, cs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bl1_sscals( beta11, alpha11 );
    bl1_sscals( beta11, alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_opd_var5( int m_AB,
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
    double*   a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;

    double*   B00      = buff_B + (0  )*cs_B + (0  )*rs_B;
    double*   b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    double*   beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    double    psi11;

    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bl1_dmult3( buff_1h, alpha11, &psi11 );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B00, a10t );
    bl1_dtrmv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_behind,
               B00,  rs_B, cs_B,
               a10t, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Her2c_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, a10t, b10t, A00 );
    bl1_dher2( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_behind,
               buff_1,
               a10t, cs_A,
               b10t, cs_B,
               A00,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Scal_external( beta11, a10t );
    bl1_dscalv( BLIS1_NO_CONJUGATE,
                m_behind,
                beta11,
                a10t, cs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bl1_dscals( beta11, alpha11 );
    bl1_dscals( beta11, alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_opc_var5( int m_AB,
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
    scomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;

    scomplex* B00      = buff_B + (0  )*cs_B + (0  )*rs_B;
    scomplex* b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    scomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    scomplex  psi11;

    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bl1_cmult3( buff_1h, alpha11, &psi11 );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B00, a10t );
    bl1_ctrmv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_behind,
               B00,  rs_B, cs_B,
               a10t, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Her2c_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, a10t, b10t, A00 );
    bl1_cher2( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_behind,
               buff_1,
               a10t, cs_A,
               b10t, cs_B,
               A00,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Scal_external( beta11, a10t );
    bl1_cscalv( BLIS1_NO_CONJUGATE,
                m_behind,
                beta11,
                a10t, cs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bl1_cscals( beta11, alpha11 );
    bl1_cscals( beta11, alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_opz_var5( int m_AB,
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
    dcomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;

    dcomplex* B00      = buff_B + (0  )*cs_B + (0  )*rs_B;
    dcomplex* b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    dcomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    dcomplex  psi11;

    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bl1_zmult3( buff_1h, alpha11, &psi11 );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B00, a10t );
    bl1_ztrmv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               m_behind,
               B00,  rs_B, cs_B,
               a10t, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Her2c_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, a10t, b10t, A00 );
    bl1_zher2( BLIS1_LOWER_TRIANGULAR,
               BLIS1_CONJUGATE,
               m_behind,
               buff_1,
               a10t, cs_A,
               b10t, cs_B,
               A00,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Scal_external( beta11, a10t );
    bl1_zscalv( BLIS1_NO_CONJUGATE,
                m_behind,
                beta11,
                a10t, cs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bl1_zscals( beta11, alpha11 );
    bl1_zscals( beta11, alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

