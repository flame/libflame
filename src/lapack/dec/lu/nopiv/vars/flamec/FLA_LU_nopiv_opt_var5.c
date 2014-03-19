/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_LU_nopiv_opt_var5( FLA_Obj A )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );

      FLA_LU_nopiv_ops_var5( m_A,
                             n_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_LU_nopiv_opd_var5( m_A,
                             n_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_LU_nopiv_opc_var5( m_A,
                             n_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_LU_nopiv_opz_var5( m_A,
                             n_A,
                             buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_nopiv_ops_var5( int m_A,
                                 int n_A,
                                 float* buff_A, int rs_A, int cs_A )
{
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   alpha11,
                   a21, rs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a21, a12t, A22 );
    bl1_sger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_ahead,
              n_ahead,
              buff_m1,
              a21, rs_A,
              a12t, cs_A,
              A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_nopiv_opd_var5( int m_A,
                                 int n_A,
                                 double* buff_A, int rs_A, int cs_A )
{
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   alpha11,
                   a21, rs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a21, a12t, A22 );
    bl1_dger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_ahead,
              n_ahead,
              buff_m1,
              a21, rs_A,
              a12t, cs_A,
              A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_nopiv_opc_var5( int m_A,
                                 int n_A,
                                 scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   alpha11,
                   a21, rs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a21, a12t, A22 );
    bl1_cger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_ahead,
              n_ahead,
              buff_m1,
              a21, rs_A,
              a12t, cs_A,
              A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_nopiv_opz_var5( int m_A,
                                 int n_A,
                                 dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( alpha11, a21 );
    bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                   m_ahead,
                   alpha11,
                   a21, rs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a21, a12t, A22 );
    bl1_zger( BLIS1_NO_CONJUGATE,
              BLIS1_NO_CONJUGATE,
              m_ahead,
              n_ahead,
              buff_m1,
              a21, rs_A,
              a12t, cs_A,
              A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

