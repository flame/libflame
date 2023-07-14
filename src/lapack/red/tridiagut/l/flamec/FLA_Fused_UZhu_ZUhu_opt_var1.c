/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fused_UZhu_ZUhu_opt_var1( FLA_Obj delta, FLA_Obj U, FLA_Obj Z, FLA_Obj t, FLA_Obj u, FLA_Obj w )
{
/*
   Effective computation:
   w = w + delta * ( U ( Z' u  ) + Z ( U' u  ) );
   t = U' u;
*/
  FLA_Datatype datatype;
  integer          m_U, n_U;
  integer          rs_U, cs_U;
  integer          rs_Z, cs_Z;
  integer          inc_u, inc_w, inc_t;

  datatype = FLA_Obj_datatype( U );

  m_U      = FLA_Obj_length( U );
  n_U      = FLA_Obj_width( U );

  rs_U     = FLA_Obj_row_stride( U );
  cs_U     = FLA_Obj_col_stride( U );

  rs_Z     = FLA_Obj_row_stride( Z );
  cs_Z     = FLA_Obj_col_stride( Z );

  inc_u    = FLA_Obj_vector_inc( u );
  
  inc_w    = FLA_Obj_vector_inc( w );

  inc_t    = FLA_Obj_vector_inc( t );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float*    buff_U     = FLA_FLOAT_PTR( U );
      float*    buff_Z     = FLA_FLOAT_PTR( Z );
      float*    buff_t     = FLA_FLOAT_PTR( t );
      float*    buff_u     = FLA_FLOAT_PTR( u );
      float*    buff_w     = FLA_FLOAT_PTR( w );
      float*    buff_delta = FLA_FLOAT_PTR( delta );

      FLA_Fused_UZhu_ZUhu_ops_var1( m_U,
                                    n_U,
                                    buff_delta,
                                    buff_U, rs_U, cs_U,
                                    buff_Z, rs_Z, cs_Z,
                                    buff_t, inc_t,
                                    buff_u, inc_u,
                                    buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_U     = FLA_DOUBLE_PTR( U );
      double*   buff_Z     = FLA_DOUBLE_PTR( Z );
      double*   buff_t     = FLA_DOUBLE_PTR( t );
      double*   buff_u     = FLA_DOUBLE_PTR( u );
      double*   buff_w     = FLA_DOUBLE_PTR( w );
      double*   buff_delta = FLA_DOUBLE_PTR( delta );

      FLA_Fused_UZhu_ZUhu_opd_var1( m_U,
                                    n_U,
                                    buff_delta,
                                    buff_U, rs_U, cs_U,
                                    buff_Z, rs_Z, cs_Z,
                                    buff_t, inc_t,
                                    buff_u, inc_u,
                                    buff_w, inc_w );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_U     = FLA_COMPLEX_PTR( U );
      scomplex* buff_Z     = FLA_COMPLEX_PTR( Z );
      scomplex* buff_t     = FLA_COMPLEX_PTR( t );
      scomplex* buff_u     = FLA_COMPLEX_PTR( u );
      scomplex* buff_w     = FLA_COMPLEX_PTR( w );
      scomplex* buff_delta = FLA_COMPLEX_PTR( delta );

      FLA_Fused_UZhu_ZUhu_opc_var1( m_U,
                                    n_U,
                                    buff_delta,
                                    buff_U, rs_U, cs_U,
                                    buff_Z, rs_Z, cs_Z,
                                    buff_u, inc_u,
                                    buff_t, inc_t,
                                    buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_U     = FLA_DOUBLE_COMPLEX_PTR( U );
      dcomplex* buff_Z     = FLA_DOUBLE_COMPLEX_PTR( Z );
      dcomplex* buff_t     = FLA_DOUBLE_COMPLEX_PTR( t );
      dcomplex* buff_u     = FLA_DOUBLE_COMPLEX_PTR( u );
      dcomplex* buff_w     = FLA_DOUBLE_COMPLEX_PTR( w );
      dcomplex* buff_delta = FLA_DOUBLE_COMPLEX_PTR( delta );

      FLA_Fused_UZhu_ZUhu_opz_var1( m_U,
                                    n_U,
                                    buff_delta,
                                    buff_U, rs_U, cs_U,
                                    buff_Z, rs_Z, cs_Z,
                                    buff_t, inc_t,
                                    buff_u, inc_u,
                                    buff_w, inc_w );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_UZhu_ZUhu_ops_var1( integer m_U,
                                        integer n_U,
                                        float* buff_delta, 
                                        float* buff_U, integer rs_U, integer cs_U, 
                                        float* buff_Z, integer rs_Z, integer cs_Z, 
                                        float* buff_t, integer inc_t, 
                                        float* buff_u, integer inc_u, 
                                        float* buff_w, integer inc_w ) 
{
  integer i;

  for ( i = 0; i < n_U; ++i )
  {
    float*    u1       = buff_U + (i  )*cs_U + (0  )*rs_U;
    float*    z1       = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    float*    delta    = buff_delta;
    float*    tau1     = buff_t + (i  )*inc_t;
    float*    u        = buff_u;
    float*    w        = buff_w;
    float     alpha;
    float     beta;

    /*------------------------------------------------------------*/

    bl1_sdot( BLIS1_CONJUGATE,
              m_U,
              z1, rs_Z,
              u,  inc_u,
              &alpha );

    bl1_sdot( BLIS1_CONJUGATE,
              m_U,
              u1, rs_U,
              u,  inc_u,
              &beta );

    *tau1 = beta;
    alpha *= *delta;
    beta  *= *delta;

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                u1, rs_U,
                w,  inc_w );

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &beta,
                z1, rs_U,
                w,  inc_w );

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_UZhu_ZUhu_opd_var1( integer m_U,
                                        integer n_U,
                                        double* buff_delta, 
                                        double* buff_U, integer rs_U, integer cs_U, 
                                        double* buff_Z, integer rs_Z, integer cs_Z, 
                                        double* buff_t, integer inc_t, 
                                        double* buff_u, integer inc_u, 
                                        double* buff_w, integer inc_w ) 
{
  double    zero  = bl1_d0();

  integer       n_run    = n_U / 2;
  integer       n_left   = n_U % 2;
  integer       step_u   = 2*cs_U;
  integer       step_z   = 2*cs_Z;
  integer       step_tau = 2*inc_t;
  integer       i;

  double*   u     = buff_u;
  double*   w     = buff_w;
  double*   u1;
  double*   u2;
  double*   z1;
  double*   z2;
  double*   tau1;
  double*   tau2;
  
  u1   = buff_U;
  u2   = buff_U +   cs_U;
  z1   = buff_Z;
  z2   = buff_Z +   cs_Z;
  tau1 = buff_t;
  tau2 = buff_t +   inc_t;

  for ( i = 0; i < n_run; ++i )
  {
    double    rho_z1u;
    double    rho_z2u;
    double    rho_u1u;
    double    rho_u2u;

    bl1_ddotsv2( BLIS1_CONJUGATE,
                 m_U,
                 z1, rs_Z,
                 z2, rs_Z,
                 u,  inc_u,
                 &zero,
                 &rho_z1u,
                 &rho_z2u );
    bl1_dneg1( &rho_z1u );
    bl1_dneg1( &rho_z2u );

    bl1_ddotv2axpyv2b( m_U,
                       u1, rs_U,
                       u2, rs_U,
                       u,  inc_u,
                       &rho_z1u,
                       &rho_z2u,
                       &rho_u1u,
                       &rho_u2u,
                       w,  inc_w );

    *tau1 = rho_u1u;
    *tau2 = rho_u2u;

    bl1_dneg1( &rho_u1u );
    bl1_dneg1( &rho_u2u );

    bl1_daxpyv2b( m_U,
                  &rho_u1u,
                  &rho_u2u,
                  z1, rs_Z,
                  z2, rs_Z,
                  w,  inc_w );

    u1   += step_u;
    u2   += step_u;
    z1   += step_z;
    z2   += step_z;
    tau1 += step_tau;
    tau2 += step_tau;
  }

  if ( n_left > 0 )
  {
    for ( i = 0; i < n_left; ++i )
    {
      double    rho_z1u;
      double    rho_u1u;

      bl1_ddot( BLIS1_CONJUGATE,
                m_U,
                z1, rs_Z,
                u,  inc_u,
                &rho_z1u );
      bl1_dneg1( &rho_z1u );

      bl1_ddotaxpy( m_U,
                    u1, rs_U,
                    u,  inc_u,
                    &rho_z1u,
                    &rho_u1u,
                    w,  inc_w );

      *tau1 = rho_u1u;

      bl1_dneg1( &rho_u1u );
      bl1_daxpyv( BLIS1_NO_CONJUGATE,
                  m_U,
                  &rho_u1u,
                  z1, rs_Z,
                  w,  inc_w );

      u1   += cs_U;
      z1   += cs_Z;
      tau1 += inc_t;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_UZhu_ZUhu_opc_var1( integer m_U,
                                        integer n_U,
                                        scomplex* buff_delta, 
                                        scomplex* buff_U, integer rs_U, integer cs_U, 
                                        scomplex* buff_Z, integer rs_Z, integer cs_Z, 
                                        scomplex* buff_t, integer inc_t, 
                                        scomplex* buff_u, integer inc_u, 
                                        scomplex* buff_w, integer inc_w ) 
{
  integer i;

  for ( i = 0; i < n_U; ++i )
  {
    scomplex* u1       = buff_U + (i  )*cs_U + (0  )*rs_U;
    scomplex* z1       = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    scomplex* delta    = buff_delta;
    scomplex* tau1     = buff_t + (i  )*inc_t;
    scomplex* u        = buff_u;
    scomplex* w        = buff_w;
    scomplex  alpha;
    scomplex  beta;

    bl1_cdot( BLIS1_CONJUGATE,
              m_U,
              z1, rs_Z,
              u,  inc_u,
              &alpha );

    bl1_cdot( BLIS1_CONJUGATE,
              m_U,
              u1, rs_U,
              u,  inc_u,
              &beta );

    *tau1 = beta;

    bl1_cscals( delta, &alpha );
    bl1_cscals( delta, &beta );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                u1, rs_U,
                w,  inc_w );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &beta,
                z1, rs_U,
                w,  inc_w );

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_UZhu_ZUhu_opz_var1( integer m_U,
                                        integer n_U,
                                        dcomplex* buff_delta, 
                                        dcomplex* buff_U, integer rs_U, integer cs_U, 
                                        dcomplex* buff_Z, integer rs_Z, integer cs_Z, 
                                        dcomplex* buff_t, integer inc_t, 
                                        dcomplex* buff_u, integer inc_u, 
                                        dcomplex* buff_w, integer inc_w ) 
{
  integer       n_run    = n_U / 1;
  integer       n_left   = n_U % 1;
  integer       step_u   = 1*cs_U;
  integer       step_z   = 1*cs_Z;
  integer       step_tau = 1*inc_t;
  integer       i;

  dcomplex* u     = buff_u;
  dcomplex* w     = buff_w;
  dcomplex* u1;
  dcomplex* z1;
  dcomplex* tau1;

  u1   = buff_U;
  z1   = buff_Z;
  tau1 = buff_t;

  for ( i = 0; i < n_run; ++i )
  {
    dcomplex  rho_z1u;
    dcomplex  rho_u1u;

/*
   Effective computation:
   w = w + delta * ( U ( Z' u  ) + Z ( U' u  ) );
*/

    bl1_zdot( BLIS1_CONJUGATE,
              m_U,
              z1, rs_Z,
              u,  inc_u,
              &rho_z1u );
    bl1_zneg1( &rho_z1u );

    bl1_zdotaxpy( m_U,
                  u1, rs_U,
                  u,  inc_u,
                  &rho_z1u,
                  &rho_u1u,
                  w,  inc_w );

    *tau1 = rho_u1u;

    bl1_zneg1( &rho_u1u );

    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &rho_u1u,
                z1, rs_Z,
                w,  inc_w );

    u1   += step_u;
    z1   += step_z;
    tau1 += step_tau;
  }

  if ( n_left == 1 )
  {
    dcomplex  rho_z1u;
    dcomplex  rho_u1u;

    bl1_zdot( BLIS1_CONJUGATE,
              m_U,
              z1, rs_Z,
              u,  inc_u,
              &rho_z1u );
    bl1_zneg1( &rho_z1u );

    bl1_zdotaxpy( m_U,
                  u1, rs_U,
                  u,  inc_u,
                  &rho_z1u,
                  &rho_u1u,
                  w,  inc_w );

    *tau1 = rho_u1u;

    bl1_zneg1( &rho_u1u );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &rho_u1u,
                z1, rs_Z,
                w,  inc_w );
  }

  return FLA_SUCCESS;
}

