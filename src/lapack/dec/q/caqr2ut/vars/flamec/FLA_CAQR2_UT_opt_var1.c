
#include "FLAME.h"

FLA_Error FLA_CAQR2_UT_opt_var1( FLA_Obj U,
                                 FLA_Obj D, FLA_Obj T )
{
  FLA_Datatype datatype;
  int          mn_UT, m_D;
  int          rs_U, cs_U;
  int          rs_D, cs_D;
  int          rs_T, cs_T;

  datatype = FLA_Obj_datatype( U );

  mn_UT    = FLA_Obj_width( U );
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

      FLA_CAQR2_UT_ops_var1( mn_UT,
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

      FLA_CAQR2_UT_opd_var1( mn_UT,
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

      FLA_CAQR2_UT_opc_var1( mn_UT,
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

      FLA_CAQR2_UT_opz_var1( mn_UT,
                             m_D,
                             buff_U, rs_U, cs_U,
                             buff_D, rs_D, cs_D,
                             buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_CAQR2_UT_ops_var1( int mn_UT,
                                 int m_D,
                                 float* buff_U, int rs_U, int cs_U,
                                 float* buff_D, int rs_D, int cs_D,
                                 float* buff_T, int rs_T, int cs_T )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  int       i, j;
  int       m_DT = m_D - mn_UT;

  for ( i = m_DT, j = 0; j < mn_UT; ++i, ++j )
  {
    float*    upsilon11 = buff_U + (j  )*cs_U + (j  )*rs_U;
    float*    u12t      = buff_U + (j+1)*cs_U + (j  )*rs_U;

    float*    D00       = buff_D + (0  )*cs_D + (0  )*rs_D;
    float*    d1        = buff_D + (j  )*cs_D + (0  )*rs_D;
    float*    D2        = buff_D + (j+1)*cs_D + (0  )*rs_D;

    float*    tau11     = buff_T + (j  )*cs_T + (j  )*rs_T;
    float*    t01       = buff_T + (j  )*cs_T + (0  )*rs_T;

    float*    d1B       = d1  + (m_DT)*rs_D;
    float*    D00B      = D00 + (m_DT)*rs_D;

    int       m_behind = i;
    int       n_behind = j;
    int       mn_ahead = mn_UT - j - 1;

    //------------------------------------------------------------//

    // FLA_Househ2_UT( FLA_LEFT,
    //                 upsilon11,
    //                 d1, tau11 );
    FLA_Househ2_UT_l_ops( m_behind + 1,
                          upsilon11,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, d1, u12t,
    //                                       D2 );
    FLA_Apply_H2_UT_l_ops_var1( m_behind + 1,
                                mn_ahead,
                                tau11,
                                d1, rs_D,
                                u12t, cs_U,
                                D2, rs_D, cs_D );

    // FLA_Copy_external( d01B, t01 );
    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    D00B, t01 );
    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, D00T, d01T, FLA_ONE, t01 );
    bl1_scopyv( BLIS1_NO_CONJUGATE,
                n_behind,
                d1B, rs_D,
                t01, rs_T );
    bl1_strmv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJ_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               n_behind,
               D00B, rs_D, cs_D,
               t01, rs_T );
    bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_DT,
               n_behind,
               buff_1,
               D00, rs_D, cs_D,
               d1,  rs_D,
               buff_1,
               t01, rs_T );

    //------------------------------------------------------------//

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_CAQR2_UT_opd_var1( int mn_UT,
                                 int m_D,
                                 double* buff_U, int rs_U, int cs_U,
                                 double* buff_D, int rs_D, int cs_D,
                                 double* buff_T, int rs_T, int cs_T )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  int       i, j;
  int       m_DT = m_D - mn_UT;

  for ( i = m_DT, j = 0; j < mn_UT; ++i, ++j )
  {
    double*   upsilon11 = buff_U + (j  )*cs_U + (j  )*rs_U;
    double*   u12t      = buff_U + (j+1)*cs_U + (j  )*rs_U;

    double*   D00       = buff_D + (0  )*cs_D + (0  )*rs_D;
    double*   d1        = buff_D + (j  )*cs_D + (0  )*rs_D;
    double*   D2        = buff_D + (j+1)*cs_D + (0  )*rs_D;

    double*   tau11     = buff_T + (j  )*cs_T + (j  )*rs_T;
    double*   t01       = buff_T + (j  )*cs_T + (0  )*rs_T;

    double*   d1B       = d1  + (m_DT)*rs_D;
    double*   D00B      = D00 + (m_DT)*rs_D;

    int       m_behind = i;
    int       n_behind = j;
    int       mn_ahead = mn_UT - j - 1;

    //------------------------------------------------------------//

    // FLA_Househ2_UT( FLA_LEFT,
    //                 upsilon11,
    //                 d1, tau11 );
    FLA_Househ2_UT_l_opd( m_behind + 1,
                          upsilon11,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, d1, u12t,
    //                                       D2 );
    FLA_Apply_H2_UT_l_opd_var1( m_behind + 1,
                                mn_ahead,
                                tau11,
                                d1, rs_D,
                                u12t, cs_U,
                                D2, rs_D, cs_D );

    // FLA_Copy_external( d01B, t01 );
    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    D00B, t01 );
    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, D00T, d01T, FLA_ONE, t01 );
    bl1_dcopyv( BLIS1_NO_CONJUGATE,
                n_behind,
                d1B, rs_D,
                t01, rs_T );
    bl1_dtrmv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJ_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               n_behind,
               D00B, rs_D, cs_D,
               t01, rs_T );
    bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_DT,
               n_behind,
               buff_1,
               D00, rs_D, cs_D,
               d1,  rs_D,
               buff_1,
               t01, rs_T );

    //------------------------------------------------------------//

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_CAQR2_UT_opc_var1( int mn_UT,
                                 int m_D,
                                 scomplex* buff_U, int rs_U, int cs_U,
                                 scomplex* buff_D, int rs_D, int cs_D,
                                 scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  int       i, j;
  int       m_DT = m_D - mn_UT;

  for ( i = m_DT, j = 0; j < mn_UT; ++i, ++j )
  {
    scomplex* upsilon11 = buff_U + (j  )*cs_U + (j  )*rs_U;
    scomplex* u12t      = buff_U + (j+1)*cs_U + (j  )*rs_U;

    scomplex* D00       = buff_D + (0  )*cs_D + (0  )*rs_D;
    scomplex* d1        = buff_D + (j  )*cs_D + (0  )*rs_D;
    scomplex* D2        = buff_D + (j+1)*cs_D + (0  )*rs_D;

    scomplex* tau11     = buff_T + (j  )*cs_T + (j  )*rs_T;
    scomplex* t01       = buff_T + (j  )*cs_T + (0  )*rs_T;

    scomplex* d1B       = d1  + (m_DT)*rs_D;
    scomplex* D00B      = D00 + (m_DT)*rs_D;

    int       m_behind = i;
    int       n_behind = j;
    int       mn_ahead = mn_UT - j - 1;

    //------------------------------------------------------------//

    // FLA_Househ2_UT( FLA_LEFT,
    //                 upsilon11,
    //                 d1, tau11 );
    FLA_Househ2_UT_l_opc( m_behind + 1,
                          upsilon11,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, d1, u12t,
    //                                       D2 );
    FLA_Apply_H2_UT_l_opc_var1( m_behind + 1,
                                mn_ahead,
                                tau11,
                                d1, rs_D,
                                u12t, cs_U,
                                D2, rs_D, cs_D );

    // FLA_Copy_external( d01B, t01 );
    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    D00B, t01 );
    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, D00T, d01T, FLA_ONE, t01 );
    bl1_ccopyv( BLIS1_NO_CONJUGATE,
                n_behind,
                d1B, rs_D,
                t01, rs_T );
    bl1_ctrmv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJ_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               n_behind,
               D00B, rs_D, cs_D,
               t01, rs_T );
    bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_DT,
               n_behind,
               buff_1,
               D00, rs_D, cs_D,
               d1,  rs_D,
               buff_1,
               t01, rs_T );

    //------------------------------------------------------------//

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_CAQR2_UT_opz_var1( int mn_UT,
                                 int m_D,
                                 dcomplex* buff_U, int rs_U, int cs_U,
                                 dcomplex* buff_D, int rs_D, int cs_D,
                                 dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  int       i, j;
  int       m_DT = m_D - mn_UT;

  for ( i = m_DT, j = 0; j < mn_UT; ++i, ++j )
  {
    dcomplex* upsilon11 = buff_U + (j  )*cs_U + (j  )*rs_U;
    dcomplex* u12t      = buff_U + (j+1)*cs_U + (j  )*rs_U;

    dcomplex* D00       = buff_D + (0  )*cs_D + (0  )*rs_D;
    dcomplex* d1        = buff_D + (j  )*cs_D + (0  )*rs_D;
    dcomplex* D2        = buff_D + (j+1)*cs_D + (0  )*rs_D;

    dcomplex* tau11     = buff_T + (j  )*cs_T + (j  )*rs_T;
    dcomplex* t01       = buff_T + (j  )*cs_T + (0  )*rs_T;

    dcomplex* d1B       = d1  + (m_DT)*rs_D;
    dcomplex* D00B      = D00 + (m_DT)*rs_D;

    int       m_behind = i;
    int       n_behind = j;
    int       mn_ahead = mn_UT - j - 1;

    //------------------------------------------------------------//

    // FLA_Househ2_UT( FLA_LEFT,
    //                 upsilon11,
    //                 d1, tau11 );
    FLA_Househ2_UT_l_opz( m_behind + 1,
                          upsilon11,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau11, d1, u12t,
    //                                       D2 );
    FLA_Apply_H2_UT_l_opz_var1( m_behind + 1,
                                mn_ahead,
                                tau11,
                                d1, rs_D,
                                u12t, cs_U,
                                D2, rs_D, cs_D );

    // FLA_Copy_external( d01B, t01 );
    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    D00B, t01 );
    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, D00T, d01T, FLA_ONE, t01 );
    bl1_zcopyv( BLIS1_NO_CONJUGATE,
                n_behind,
                d1B, rs_D,
                t01, rs_T );
    bl1_ztrmv( BLIS1_UPPER_TRIANGULAR,
               BLIS1_CONJ_TRANSPOSE,
               BLIS1_NONUNIT_DIAG,
               n_behind,
               D00B, rs_D, cs_D,
               t01, rs_T );
    bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_DT,
               n_behind,
               buff_1,
               D00, rs_D, cs_D,
               d1,  rs_D,
               buff_1,
               t01, rs_T );

    //------------------------------------------------------------//

  }

  return FLA_SUCCESS;
}

