
#include "FLAME.h"

FLA_Error FLA_Accum_T_UT_fc_opt_var1( FLA_Obj A, FLA_Obj t, FLA_Obj T )
{
  FLA_Datatype datatype;
  int          m_A, n_AT;
  int          rs_A, cs_A;
  int          m_t, inc_t;
  int          rs_T, cs_T;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_AT     = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_t      = FLA_Obj_vector_dim( t );
  inc_t    = FLA_Obj_vector_inc( t );

  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_t = FLA_FLOAT_PTR( t );
      float* buff_T = FLA_FLOAT_PTR( T );

      FLA_Accum_T_UT_fc_ops_var1( m_A,
                                  n_AT,
                                  buff_A, rs_A, cs_A,
                                  m_t, 
                                  buff_t, inc_t,
                                  buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_t = FLA_DOUBLE_PTR( t );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_Accum_T_UT_fc_opd_var1( m_A,
                                  n_AT,
                                  buff_A, rs_A, cs_A,
                                  m_t, 
                                  buff_t, inc_t,
                                  buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_t = FLA_COMPLEX_PTR( t );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_Accum_T_UT_fc_opc_var1( m_A,
                                  n_AT,
                                  buff_A, rs_A, cs_A,
                                  m_t, 
                                  buff_t, inc_t,
                                  buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_t = FLA_DOUBLE_COMPLEX_PTR( t );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_Accum_T_UT_fc_opz_var1( m_A,
                                  n_AT,
                                  buff_A, rs_A, cs_A,
                                  m_t,
                                  buff_t, inc_t,
                                  buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Accum_T_UT_fc_ops_var1( int m_A,
                                      int n_AT,
                                      float* buff_A, int rs_A, int cs_A,
                                      int m_t,
                                      float* buff_t, int inc_t,
                                      float* buff_T, int rs_T, int cs_T )
{
  float* buff_1   = FLA_FLOAT_PTR( FLA_ONE );
  int    i;

  for ( i = 0; i < m_t; ++i )
  {
    float* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    float* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    float* tau1     = buff_t + (i  )*inc_t;

    float* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    float* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int    m_ahead  = m_A  - i - 1;
    int    n_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( tau1, tau11 );
    *tau11 = *tau1;

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    bl1_scopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01, rs_T );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_sgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Accum_T_UT_fc_opd_var1( int m_A,
                                      int n_AT,
                                      double* buff_A, int rs_A, int cs_A,
                                      int m_t, 
                                      double* buff_t, int inc_t,
                                      double* buff_T, int rs_T, int cs_T )
{
  double* buff_1   = FLA_DOUBLE_PTR( FLA_ONE );
  int     i;

  for ( i = 0; i < m_t; ++i )
  {
    double* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    double* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    double* tau1      = buff_t + (i  )*inc_t;

    double* tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    double* t01       = buff_T + (i  )*cs_T + (0  )*rs_T;

    int     m_ahead   = m_A  - i - 1;
    int     n_behind  = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( tau1, tau11 );
    *tau11 = *tau1;

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    bl1_dcopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01, rs_T );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_dgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               n_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Accum_T_UT_fc_opc_var1( int m_A,
                                      int n_AT,
                                      scomplex* buff_A, int rs_A, int cs_A,
                                      int m_t,
                                      scomplex* buff_t, int inc_t,
                                      scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_1   = FLA_COMPLEX_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < m_t; ++i )
  {
    scomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    scomplex* tau1      = buff_t + (i  )*inc_t;

    scomplex* tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    scomplex* t01       = buff_T + (i  )*cs_T + (0  )*rs_T;

    int       m_ahead   = m_A  - i - 1;
    int       n_behind  = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( tau1, tau11 );
    *tau11 = *tau1;

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    bl1_ccopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01, rs_T );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_cgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead, n_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Accum_T_UT_fc_opz_var1( int m_A,
                                      int n_AT,
                                      dcomplex* buff_A, int rs_A, int cs_A,
                                      int m_t,
                                      dcomplex* buff_t, int inc_t,
                                      dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_1   = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < m_t; ++i )
  {
    dcomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    dcomplex* tau1      = buff_t + (i  )*inc_t;

    dcomplex* tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    dcomplex* t01       = buff_T + (i  )*cs_T + (0  )*rs_T;

    int       m_ahead   = m_A  - i - 1;
    int       n_behind  = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( tau1, tau11 );
    *tau11 = *tau1;

    // FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    bl1_zcopyv( BLIS1_CONJUGATE,
                n_behind,
                a10t, cs_A,
                t01, rs_T );

    // FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );
    bl1_zgemv( BLIS1_CONJ_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead, n_behind,
               buff_1,
               A20, rs_A, cs_A,
               a21, rs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}
