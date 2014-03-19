/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fused_Ahx_Ax_opt_var1( FLA_Obj A, FLA_Obj x, FLA_Obj v, FLA_Obj w )
{
/*
   Effective computation:
   v = A' * x;
   w = A  * x;
*/
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_x, inc_v, inc_w;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_x    = FLA_Obj_vector_inc( x );

  inc_v    = FLA_Obj_vector_inc( v );

  inc_w    = FLA_Obj_vector_inc( w );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_x = FLA_FLOAT_PTR( x );
      float* buff_v = FLA_FLOAT_PTR( v );
      float* buff_w = FLA_FLOAT_PTR( w );

      FLA_Fused_Ahx_Ax_ops_var1( m_A,
                                 n_A,
                                 buff_A, rs_A, cs_A,
                                 buff_x, inc_x,
                                 buff_v, inc_v,
                                 buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_x = FLA_DOUBLE_PTR( x );
      double* buff_v = FLA_DOUBLE_PTR( v );
      double* buff_w = FLA_DOUBLE_PTR( w );

      FLA_Fused_Ahx_Ax_opd_var1( m_A,
                                 n_A,
                                 buff_A, rs_A, cs_A,
                                 buff_x, inc_x,
                                 buff_v, inc_v,
                                 buff_w, inc_w );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_x = FLA_COMPLEX_PTR( x );
      scomplex* buff_v = FLA_COMPLEX_PTR( v );
      scomplex* buff_w = FLA_COMPLEX_PTR( w );

      FLA_Fused_Ahx_Ax_opc_var1( m_A,
                                 n_A,
                                 buff_A, rs_A, cs_A,
                                 buff_x, inc_x,
                                 buff_v, inc_v,
                                 buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_x = FLA_DOUBLE_COMPLEX_PTR( x );
      dcomplex* buff_v = FLA_DOUBLE_COMPLEX_PTR( v );
      dcomplex* buff_w = FLA_DOUBLE_COMPLEX_PTR( w );

      FLA_Fused_Ahx_Ax_opz_var1( m_A,
                                 n_A,
                                 buff_A, rs_A, cs_A,
                                 buff_x, inc_x,
                                 buff_v, inc_v,
                                 buff_w, inc_w );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Ahx_Ax_ops_var1( int m_A,
                                     int n_A,
                                     float* buff_A, int rs_A, int cs_A, 
                                     float* buff_x, int inc_x, 
                                     float* buff_v, int inc_v, 
                                     float* buff_w, int inc_w )
{
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  int       i;

  bl1_ssetv( m_A,
             buff_0,
             buff_w, inc_w );

  for ( i = 0; i < n_A; ++i )
  {
    float*    a1       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    nu1      = buff_v + (i  )*inc_v;
    float*    x        = buff_x;
    float*    chi1     = buff_x + (i  )*inc_x;
    float*    w        = buff_w;

    /*------------------------------------------------------------*/

    bl1_sdot( BLIS1_CONJUGATE,
              m_A,
              a1, rs_A,
              x,  inc_x,
              nu1 );
/*
    *nu1 = F77_sdot( &m_A,
                     a1, &rs_A,
                     x,  &inc_x );
*/

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                chi1,
                a1, rs_A,
                w,  inc_w );
/*
    F77_saxpy( &m_A,
               chi1,
               a1, &rs_A,
               w,  &inc_w );
*/

    /*------------------------------------------------------------*/

  }


  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Ahx_Ax_opd_var1( int m_A,
                                     int n_A,
                                     double* buff_A, int rs_A, int cs_A, 
                                     double* buff_x, int inc_x, 
                                     double* buff_v, int inc_v, 
                                     double* buff_w, int inc_w )
{
  double    zero  = bl1_d0();
  int       i;

  double*   restrict w = buff_w;
  double*   restrict x = buff_x;

  double*   restrict a1;
  double*   restrict a2;
  double*   restrict nu1;
  double*   restrict nu2;
  double*   restrict chi1;
  double*   restrict chi2;

  int       n_run     = n_A / 2;
  int       n_left    = n_A % 2;
  int       step_a1   = 2*cs_A;
  int       step_nu1  = 2*inc_v;
  int       step_chi1 = 2*inc_x;

  bl1_dsetv( m_A,
             &zero,
             buff_w, inc_w );

  a1   = buff_A;
  a2   = buff_A + cs_A;
  nu1  = buff_v;
  nu2  = buff_v + inc_v;
  chi1 = buff_x;
  chi2 = buff_x + inc_x;

  for ( i = 0; i < n_run; ++i )
  {
    /*------------------------------------------------------------*/

    bl1_ddotv2axpyv2b( m_A,
                       a1, rs_A,
                       a2, rs_A,
                       x,  inc_x,
                       chi1,
                       chi2,
                       nu1,
                       nu2,
                       w,  inc_w );

    /*------------------------------------------------------------*/

    a1   += step_a1;
    a2   += step_a1;
    nu1  += step_nu1;
    nu2  += step_nu1;
    chi1 += step_chi1;
    chi2 += step_chi1;
  }

  if ( n_left > 0 )
  {
    for ( i = 0; i < n_left; ++i )
    {
      bl1_ddotaxpy( m_A,
                    a1, rs_A,
                    x,  inc_x,
                    chi1,
                    nu1,
                    w,  inc_w );

      a1   += rs_A;
      nu1  += inc_v;
      chi1 += inc_x;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Ahx_Ax_opc_var1( int m_A,
                                     int n_A,
                                     scomplex* buff_A, int rs_A, int cs_A, 
                                     scomplex* buff_x, int inc_x, 
                                     scomplex* buff_v, int inc_v, 
                                     scomplex* buff_w, int inc_w )
{
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  int       i;

  bl1_csetv( m_A,
             buff_0,
             buff_w, inc_w );

  for ( i = 0; i < n_A; ++i )
  {
    scomplex* a1       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* nu1      = buff_v + (i  )*inc_v;
    scomplex* x        = buff_x;
    scomplex* chi1     = buff_x + (i  )*inc_x;
    scomplex* w        = buff_w;

    /*------------------------------------------------------------*/

    bl1_cdot( BLIS1_CONJUGATE,
              m_A,
              a1, rs_A,
              x,  inc_x,
              nu1 );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                chi1,
                a1, rs_A,
                w,  inc_w );
/*
    F77_caxpy( &m_A,
               chi1,
               a1, &rs_A,
               w,  &inc_w );
*/

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Ahx_Ax_opz_var1( int m_A,
                                     int n_A,
                                     dcomplex* buff_A, int rs_A, int cs_A, 
                                     dcomplex* buff_x, int inc_x, 
                                     dcomplex* buff_v, int inc_v, 
                                     dcomplex* buff_w, int inc_w )
{
  dcomplex  zero  = bl1_z0();
  int       i;

  dcomplex* restrict w = buff_w;
  dcomplex* restrict x = buff_x;

  dcomplex* restrict a1;
  dcomplex* restrict a2;
  dcomplex* restrict nu1;
  dcomplex* restrict nu2;
  dcomplex* restrict chi1;
  dcomplex* restrict chi2;

  int       n_run     = n_A / 2;
  int       n_left    = n_A % 2;
  int       step_a1   = 2*cs_A;
  int       step_nu1  = 2*inc_v;
  int       step_chi1 = 2*inc_x;

  bl1_zsetv( m_A,
             &zero,
             buff_w, inc_w );

  a1   = buff_A;
  a2   = buff_A + cs_A;
  nu1  = buff_v;
  nu2  = buff_v + inc_v;
  chi1 = buff_x;
  chi2 = buff_x + inc_x;

  for ( i = 0; i < n_run; ++i )
  {
    /*------------------------------------------------------------*/

/*
    bl1_zdotaxpy( m_A,
                  a1, rs_A,
                  x,  inc_x,
                  chi1,
                  nu1,
                  w,  inc_w );
*/

    bl1_zdotv2axpyv2b( m_A,
                       a1, rs_A,
                       a2, rs_A,
                       x,  inc_x,
                       chi1,
                       chi2,
                       nu1,
                       nu2,
                       w,  inc_w );

    /*------------------------------------------------------------*/

    a1   += step_a1;
    a2   += step_a1;
    nu1  += step_nu1;
    nu2  += step_nu1;
    chi1 += step_chi1;
    chi2 += step_chi1;
  }

  if ( n_left > 0 )
  {
    for ( i = 0; i < n_left; ++i )
    {
      bl1_zdotaxpy( m_A,
                    a1, rs_A,
                    x,  inc_x,
                    chi1,
                    nu1,
                    w,  inc_w );

      a1   += rs_A;
      nu1  += inc_v;
      chi1 += inc_x;
    }
  }

  return FLA_SUCCESS;
}

