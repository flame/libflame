/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Sylv_hn_opt_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  FLA_Datatype datatype;
  int          m_C, n_C;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  int          rs_C, cs_C;
  int          info;

  datatype = FLA_Obj_datatype( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );

  m_C      = FLA_Obj_length( C );
  n_C      = FLA_Obj_width( C );
  rs_C     = FLA_Obj_row_stride( C );
  cs_C     = FLA_Obj_col_stride( C );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      int*   buff_isgn  = FLA_INT_PTR( isgn );
      float* buff_A     = FLA_FLOAT_PTR( A );
      float* buff_B     = FLA_FLOAT_PTR( B );
      float* buff_C     = FLA_FLOAT_PTR( C );
      float* buff_scale = FLA_FLOAT_PTR( scale );
      float  sgn        = ( float ) *buff_isgn;

      FLA_Sylv_hn_ops_var1( sgn,
                            m_C,
                            n_C,
                            buff_A, rs_A, cs_A,
                            buff_B, rs_B, cs_B,
                            buff_C, rs_C, cs_C,
                            buff_scale,
                            &info );

      break;
    }

    case FLA_DOUBLE:
    {
      int*    buff_isgn  = FLA_INT_PTR( isgn );
      double* buff_A     = FLA_DOUBLE_PTR( A );
      double* buff_B     = FLA_DOUBLE_PTR( B );
      double* buff_C     = FLA_DOUBLE_PTR( C );
      double* buff_scale = FLA_DOUBLE_PTR( scale );
      double  sgn        = ( double ) *buff_isgn;

      FLA_Sylv_hn_opd_var1( sgn,
                            m_C,
                            n_C,
                            buff_A, rs_A, cs_A,
                            buff_B, rs_B, cs_B,
                            buff_C, rs_C, cs_C,
                            buff_scale,
                            &info );

      break;
    }

    case FLA_COMPLEX:
    {
      int*      buff_isgn  = FLA_INT_PTR( isgn );
      scomplex* buff_A     = FLA_COMPLEX_PTR( A );
      scomplex* buff_B     = FLA_COMPLEX_PTR( B );
      scomplex* buff_C     = FLA_COMPLEX_PTR( C );
      scomplex* buff_scale = FLA_COMPLEX_PTR( scale );
      float     sgn        = ( float ) *buff_isgn;

      FLA_Sylv_hn_opc_var1( sgn,
                            m_C,
                            n_C,
                            buff_A, rs_A, cs_A,
                            buff_B, rs_B, cs_B,
                            buff_C, rs_C, cs_C,
                            buff_scale,
                            &info );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      int*      buff_isgn  = FLA_INT_PTR( isgn );
      dcomplex* buff_A     = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_B     = FLA_DOUBLE_COMPLEX_PTR( B );
      dcomplex* buff_C     = FLA_DOUBLE_COMPLEX_PTR( C );
      dcomplex* buff_scale = FLA_DOUBLE_COMPLEX_PTR( scale );
      double    sgn        = ( double ) *buff_isgn;

      FLA_Sylv_hn_opz_var1( sgn,
                            m_C,
                            n_C,
                            buff_A, rs_A, cs_A,
                            buff_B, rs_B, cs_B,
                            buff_C, rs_C, cs_C,
                            buff_scale,
                            &info );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Sylv_hn_ops_var1( float sgn,
                                int m_C,
                                int n_C,
                                float* buff_A, int rs_A, int cs_A,
                                float* buff_B, int rs_B, int cs_B,
                                float* buff_C, int rs_C, int cs_C,
                                float* buff_scale,
                                int* info )
{
  int l, k;

  for ( l = 0; l < n_C; l++ )
  {
    for ( k = 0; k < m_C; k++ )
    {
      float*    a01      = buff_A + (k  )*cs_A + (0  )*rs_A;
      float*    b01      = buff_B + (l  )*cs_B + (0  )*rs_B;
      float*    c10t     = buff_C + (0  )*cs_C + (k  )*rs_C;
      float*    c01      = buff_C + (l  )*cs_C + (0  )*rs_C;
      float*    alpha11  = buff_A + (k  )*cs_A + (k  )*rs_A;
      float*    beta11   = buff_B + (l  )*cs_B + (l  )*rs_B;
      float*    ckl      = buff_C + (l  )*cs_C + (k  )*rs_C;
      float     suml;
      float     sumr;
      float     vec;
      float     a11;
      float     x11;

      int       m_behind = k;
      int       n_behind = l;

      /*------------------------------------------------------------*/

      bl1_sdot( BLIS1_CONJUGATE,
                m_behind,
                a01, rs_A,
                c01, rs_C,
                &suml );

      bl1_sdot( BLIS1_NO_CONJUGATE,
                n_behind,
                c10t, cs_C,
                b01, rs_B,
                &sumr );

      vec = (*ckl) - ( suml + sgn * sumr );

      a11 = (*alpha11) + sgn * (*beta11);

      bl1_sdiv3( &vec, &a11, &x11 );

      *ckl = x11;

      /*------------------------------------------------------------*/

    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Sylv_hn_opd_var1( double sgn,
                                int m_C,
                                int n_C,
                                double* buff_A, int rs_A, int cs_A,
                                double* buff_B, int rs_B, int cs_B,
                                double* buff_C, int rs_C, int cs_C,
                                double* buff_scale,
                                int* info )
{
  int l, k;

  for ( l = 0; l < n_C; l++ )
  {
    for ( k = 0; k < m_C; k++ )
    {
      double*   a01      = buff_A + (k  )*cs_A + (0  )*rs_A;
      double*   b01      = buff_B + (l  )*cs_B + (0  )*rs_B;
      double*   c10t     = buff_C + (0  )*cs_C + (k  )*rs_C;
      double*   c01      = buff_C + (l  )*cs_C + (0  )*rs_C;
      double*   alpha11  = buff_A + (k  )*cs_A + (k  )*rs_A;
      double*   beta11   = buff_B + (l  )*cs_B + (l  )*rs_B;
      double*   ckl      = buff_C + (l  )*cs_C + (k  )*rs_C;
      double    suml;
      double    sumr;
      double    vec;
      double    a11;
      double    x11;

      int       m_behind = k;
      int       n_behind = l;

      /*------------------------------------------------------------*/

      bl1_ddot( BLIS1_CONJUGATE,
                m_behind,
                a01, rs_A,
                c01, rs_C,
                &suml );

      bl1_ddot( BLIS1_NO_CONJUGATE,
                n_behind,
                c10t, cs_C,
                b01, rs_B,
                &sumr );

      vec = (*ckl) - ( suml + sgn * sumr );

      a11 = (*alpha11) + sgn * (*beta11);

      bl1_ddiv3( &vec, &a11, &x11 );

      *ckl = x11;

      /*------------------------------------------------------------*/

    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Sylv_hn_opc_var1( float sgn,
                                int m_C,
                                int n_C,
                                scomplex* buff_A, int rs_A, int cs_A,
                                scomplex* buff_B, int rs_B, int cs_B,
                                scomplex* buff_C, int rs_C, int cs_C,
                                scomplex* buff_scale,
                                int* info )
{
  int l, k;

  for ( l = 0; l < n_C; l++ )
  {
    for ( k = 0; k < m_C; k++ )
    {
      scomplex* a01      = buff_A + (k  )*cs_A + (0  )*rs_A;
      scomplex* b01      = buff_B + (l  )*cs_B + (0  )*rs_B;
      scomplex* c10t     = buff_C + (0  )*cs_C + (k  )*rs_C;
      scomplex* c01      = buff_C + (l  )*cs_C + (0  )*rs_C;
      scomplex* alpha11  = buff_A + (k  )*cs_A + (k  )*rs_A;
      scomplex* beta11   = buff_B + (l  )*cs_B + (l  )*rs_B;
      scomplex* ckl      = buff_C + (l  )*cs_C + (k  )*rs_C;
      scomplex  suml;
      scomplex  sumr;
      scomplex  vec;
      scomplex  a11;
      scomplex  x11;

      int       m_behind = k;
      int       n_behind = l;

      /*------------------------------------------------------------*/

      bl1_cdot( BLIS1_CONJUGATE,
                m_behind,
                a01, rs_A,
                c01, rs_C,
                &suml );

      bl1_cdot( BLIS1_NO_CONJUGATE,
                n_behind,
                c10t, cs_C,
                b01, rs_B,
                &sumr );

      vec.real = ckl->real - ( suml.real + sgn * sumr.real );
      vec.imag = ckl->imag - ( suml.imag + sgn * sumr.imag );

      a11.real =  alpha11->real + sgn * beta11->real;
      a11.imag = -alpha11->imag + sgn * beta11->imag;

      bl1_cdiv3( &vec, &a11, &x11 );

      *ckl = x11;

      /*------------------------------------------------------------*/

    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Sylv_hn_opz_var1( double sgn,
                                int m_C,
                                int n_C,
                                dcomplex* buff_A, int rs_A, int cs_A,
                                dcomplex* buff_B, int rs_B, int cs_B,
                                dcomplex* buff_C, int rs_C, int cs_C,
                                dcomplex* buff_scale,
                                int* info )
{
  int l, k;

  for ( l = 0; l < n_C; l++ )
  {
    for ( k = 0; k < m_C; k++ )
    {
      dcomplex* a01      = buff_A + (k  )*cs_A + (0  )*rs_A;
      dcomplex* b01      = buff_B + (l  )*cs_B + (0  )*rs_B;
      dcomplex* c10t     = buff_C + (0  )*cs_C + (k  )*rs_C;
      dcomplex* c01      = buff_C + (l  )*cs_C + (0  )*rs_C;
      dcomplex* alpha11  = buff_A + (k  )*cs_A + (k  )*rs_A;
      dcomplex* beta11   = buff_B + (l  )*cs_B + (l  )*rs_B;
      dcomplex* ckl      = buff_C + (l  )*cs_C + (k  )*rs_C;
      dcomplex  suml;
      dcomplex  sumr;
      dcomplex  vec;
      dcomplex  a11;
      dcomplex  x11;

      int       m_behind = k;
      int       n_behind = l;

      /*------------------------------------------------------------*/

      bl1_zdot( BLIS1_CONJUGATE,
                m_behind,
                a01, rs_A,
                c01, rs_C,
                &suml );

      bl1_zdot( BLIS1_NO_CONJUGATE,
                n_behind,
                c10t, cs_C,
                b01, rs_B,
                &sumr );

      vec.real = ckl->real - ( suml.real + sgn * sumr.real );
      vec.imag = ckl->imag - ( suml.imag + sgn * sumr.imag );

      a11.real =  alpha11->real + sgn * beta11->real;
      a11.imag = -alpha11->imag + sgn * beta11->imag;

      bl1_zdiv3( &vec, &a11, &x11 );

      *ckl = x11;

      /*------------------------------------------------------------*/

    }
  }

  return FLA_SUCCESS;
}

