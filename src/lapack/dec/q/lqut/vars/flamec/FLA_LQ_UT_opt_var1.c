
#include "FLAME.h"

FLA_Error FLA_LQ_UT_opt_var1( FLA_Obj A, FLA_Obj t )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_t;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_t    = FLA_Obj_vector_inc( t );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_t = FLA_FLOAT_PTR( t );

      FLA_LQ_UT_ops_var1( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_t, inc_t );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_t = FLA_DOUBLE_PTR( t );

      FLA_LQ_UT_opd_var1( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_t, inc_t );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_t = FLA_COMPLEX_PTR( t );

      FLA_LQ_UT_opc_var1( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_t, inc_t );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_t = FLA_DOUBLE_COMPLEX_PTR( t );

      FLA_LQ_UT_opz_var1( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_t, inc_t );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LQ_UT_ops_var1( int m_A,
                              int n_A,
                              float* buff_A, int rs_A, int cs_A,
                              float* buff_t, int inc_t )
{
  int min_m_n = min( m_A, n_A );
  int i;

  for ( i = 0; i < min_m_n; ++i )
  {
    float* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float* tau1     = buff_t + (i  )*inc_t;

    int    m_ahead  = m_A - i - 1;
    int    n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_RIGHT, alpha11, a12t
    //                 tau1 );
    FLA_Househ2_UT_r_ops( n_ahead,
                          alpha11,
                          a12t, cs_A,
                          tau1 );

    // FLA_Apply_H2_UT( FLA_RIGHT, tau1, a12t, a21, A22 );
    FLA_Apply_H2_UT_r_ops_var1( m_ahead,
                                n_ahead,
                                tau1,
                                a12t, cs_A,
                                a21, rs_A,
                                A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LQ_UT_opd_var1( int m_A,
                              int n_A,
                              double* buff_A, int rs_A, int cs_A,
                              double* buff_t, int inc_t )
{
  int min_m_n = min( m_A, n_A );
  int i;

  for ( i = 0; i < min_m_n; ++i )
  {
    double* alpha11 = buff_A + (i  )*cs_A + (i  )*rs_A;
    double* a21     = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double* a12t    = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double* A22     = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double* tau1    = buff_t + (i  )*inc_t;

    int     m_ahead  = m_A - i - 1;
    int     n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_RIGHT, alpha11, a12t
    //                 tau1 );
    FLA_Househ2_UT_r_opd( n_ahead,
                          alpha11,
                          a12t, cs_A,
                          tau1 );

    // FLA_Apply_H2_UT( FLA_RIGHT, tau1, a12t, a21, A22 );
    FLA_Apply_H2_UT_r_opd_var1( m_ahead,
                                n_ahead,
                                tau1,
                                a12t, cs_A,
                                a21, rs_A,
                                A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LQ_UT_opc_var1( int m_A,
                              int n_A,
                              scomplex* buff_A, int rs_A, int cs_A,
                              scomplex* buff_t, int inc_t )
{
  int min_m_n = min( m_A, n_A );
  int i;

  for ( i = 0; i < min_m_n; ++i )
  {
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* tau1     = buff_t + (i  )*inc_t;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_RIGHT, alpha11, a12t
    //                 tau1 );
    FLA_Househ2_UT_r_opc( n_ahead,
                          alpha11,
                          a12t, cs_A,
                          tau1 );

    // FLA_Apply_H2_UT( FLA_RIGHT, tau1, a12t, a21, A22 );
    FLA_Apply_H2_UT_r_opc_var1( m_ahead,
                                n_ahead,
                                tau1,
                                a12t, cs_A,
                                a21, rs_A,
                                A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LQ_UT_opz_var1( int m_A,
                              int n_A,
                              dcomplex* buff_A, int rs_A, int cs_A,
                              dcomplex* buff_t, int inc_t )
{
  int min_m_n = min( m_A, n_A );
  int i;

  for ( i = 0; i < min_m_n; ++i )
  {
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* tau1     = buff_t + (i  )*inc_t;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_RIGHT, alpha11, a12t
    //                 tau1 );
    FLA_Househ2_UT_r_opz( n_ahead,
                          alpha11,
                          a12t, cs_A,
                          tau1 );

    // FLA_Apply_H2_UT( FLA_RIGHT, tau1, a12t, a21, A22 );
    FLA_Apply_H2_UT_r_opz_var1( m_ahead,
                                n_ahead,
                                tau1,
                                a12t, cs_A,
                                a21, rs_A,
                                A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}
