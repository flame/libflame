
#include "FLAME.h"

FLA_Error FLA_Fused_UZhu_ZUhu_opt_var1( FLA_Obj delta, FLA_Obj U, FLA_Obj Z, FLA_Obj t, FLA_Obj u, FLA_Obj w )
{
/*
   Effective computation:
   w = w + delta * ( U ( Z' u  ) + Z ( U' u  ) );
   t = U' u;
*/
  FLA_Datatype datatype;
  int          m_U, n_U;
  int          rs_U, cs_U;
  int          rs_Z, cs_Z;
  int          inc_u, inc_w, inc_t;

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



FLA_Error FLA_Fused_UZhu_ZUhu_ops_var1( int m_U,
                                        int n_U,
                                        float* buff_delta, 
                                        float* buff_U, int rs_U, int cs_U, 
                                        float* buff_Z, int rs_Z, int cs_Z, 
                                        float* buff_t, int inc_t, 
                                        float* buff_u, int inc_u, 
                                        float* buff_w, int inc_w ) 
{
  int i;

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
/*
    alpha = F77_sdot( &m_U,
                      z1, &rs_Z,
                      u,  &inc_u );
*/

    bl1_sdot( BLIS1_CONJUGATE,
              m_U,
              u1, rs_U,
              u,  inc_u,
              &beta );
/*
    beta = F77_sdot( &m_U,
                     u1, &rs_U,
                     u,  &inc_u );
*/

    *tau1 = beta;

    // bl1_sscals( delta, &alpha );
    // bl1_sscals( delta, &beta );
    alpha *= *delta;
    beta  *= *delta;

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                u1, rs_U,
                w,  inc_w );
/*
    F77_saxpy( &m_U,
               &alpha,
               u1, &rs_U,
               w,  &inc_w );
*/

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &beta,
                z1, rs_U,
                w,  inc_w );
/*
    F77_saxpy( &m_U,
               &beta,
               z1, &rs_Z,
               w,  &inc_w );
*/

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_UZhu_ZUhu_opd_var1( int m_U,
                                        int n_U,
                                        double* buff_delta, 
                                        double* buff_U, int rs_U, int cs_U, 
                                        double* buff_Z, int rs_Z, int cs_Z, 
                                        double* buff_t, int inc_t, 
                                        double* buff_u, int inc_u, 
                                        double* buff_w, int inc_w ) 
{
  double    zero  = bl1_d0();

  int       n_run    = n_U / 2;
  int       n_left   = n_U % 2;
  int       step_u   = 2*cs_U;
  int       step_z   = 2*cs_Z;
  int       step_tau = 2*inc_t;
  int       i;

  double*   u     = buff_u;
  double*   w     = buff_w;
  //double*   delta = buff_delta;

  double*   u1;
  double*   u2;
  double*   u3;
  double*   z1;
  double*   z2;
  double*   z3;
  double*   tau1;
  double*   tau2;
  double*   tau3;

  u1   = buff_U;
  u2   = buff_U +   cs_U;
  u3   = buff_U + 2*cs_U;
  z1   = buff_Z;
  z2   = buff_Z +   cs_Z;
  z3   = buff_Z + 2*cs_Z;
  tau1 = buff_t;
  tau2 = buff_t +   inc_t;
  tau3 = buff_t + 2*inc_t;

  for ( i = 0; i < n_run; ++i )
  {
    double    rho_z1u;
    double    rho_z2u;
    //double    rho_z3u;
    double    rho_u1u;
    double    rho_u2u;
    //double    rho_u3u;

    /*------------------------------------------------------------*/
/*
    bl1_ddotsv3( BLIS1_CONJUGATE,
                 m_U,
                 z1, rs_Z,
                 z2, rs_Z,
                 z3, rs_Z,
                 u,  inc_u,
                 &zero,
                 &rho_z1u,
                 &rho_z2u,
                 &rho_z3u );
    bl1_dneg1( &rho_z1u );
    bl1_dneg1( &rho_z2u );
    bl1_dneg1( &rho_z3u );

    bl1_ddotv2axpyv2b( m_U,
                       u1, rs_U,
                       u2, rs_U,
                       u,  inc_u,
                       &rho_z1u,
                       &rho_z2u,
                       &rho_u1u,
                       &rho_u2u,
                       w,  inc_w );
    bl1_ddotaxpy( m_U,
                  u3, rs_U,
                  u,  inc_u,
                  &rho_z3u,
                  &rho_u3u,
                  w,  inc_w );

    *tau1 = rho_u1u;
    *tau2 = rho_u2u;
    *tau3 = rho_u3u;

    bl1_dneg1( &rho_u1u );
    bl1_dneg1( &rho_u2u );
    bl1_dneg1( &rho_u3u );

    bl1_daxpyv3b( m_U,
                  &rho_u1u,
                  &rho_u2u,
                  &rho_u3u,
                  z1, rs_Z,
                  z2, rs_Z,
                  z3, rs_Z,
                  w,  inc_w );
*/
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


    /*------------------------------------------------------------*/

    u1   += step_u;
    u2   += step_u;
    u3   += step_u;
    z1   += step_z;
    z2   += step_z;
    z3   += step_z;
    tau1 += step_tau;
    tau2 += step_tau;
    tau3 += step_tau;
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



FLA_Error FLA_Fused_UZhu_ZUhu_opc_var1( int m_U,
                                        int n_U,
                                        scomplex* buff_delta, 
                                        scomplex* buff_U, int rs_U, int cs_U, 
                                        scomplex* buff_Z, int rs_Z, int cs_Z, 
                                        scomplex* buff_t, int inc_t, 
                                        scomplex* buff_u, int inc_u, 
                                        scomplex* buff_w, int inc_w ) 
{
  int i;

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

    /*------------------------------------------------------------*/

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
/*
    F77_caxpy( &m_U,
               &alpha,
               u1, &rs_U,
               w,  &inc_w );
*/

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &beta,
                z1, rs_U,
                w,  inc_w );
/*
    F77_caxpy( &m_U,
               &beta,
               z1, &rs_Z,
               w,  &inc_w );
*/

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_UZhu_ZUhu_opz_var1( int m_U,
                                        int n_U,
                                        dcomplex* buff_delta, 
                                        dcomplex* buff_U, int rs_U, int cs_U, 
                                        dcomplex* buff_Z, int rs_Z, int cs_Z, 
                                        dcomplex* buff_t, int inc_t, 
                                        dcomplex* buff_u, int inc_u, 
                                        dcomplex* buff_w, int inc_w ) 
{
  //dcomplex  zero  = bl1_z0();

  int       n_run    = n_U / 1;
  int       n_left   = n_U % 1;
  int       step_u   = 1*cs_U;
  int       step_z   = 1*cs_Z;
  int       step_tau = 1*inc_t;
  int       i;

  dcomplex* u     = buff_u;
  dcomplex* w     = buff_w;
  //dcomplex* delta = buff_delta;

  dcomplex* u1;
  dcomplex* u2;
  dcomplex* z1;
  dcomplex* z2;
  dcomplex* tau1;
  dcomplex* tau2;

  u1   = buff_U;
  u2   = buff_U + cs_U;
  z1   = buff_Z;
  z2   = buff_Z + cs_Z;
  tau1 = buff_t;
  tau2 = buff_t + inc_t;

  for ( i = 0; i < n_run; ++i )
  {
    dcomplex  rho_z1u;
    //dcomplex  rho_z2u;
    dcomplex  rho_u1u;
    //dcomplex  rho_u2u;

    /*------------------------------------------------------------*/

/*
   Effective computation:
   w = w + delta * ( U ( Z' u  ) + Z ( U' u  ) );
*/

/*
    bl1_zdotsv2( BLIS1_CONJUGATE,
                 m_U,
                 z1, rs_Z,
                 u1, rs_U,
                 u,  inc_u,
                 &zero,
                 &rho_z1u,
                 &rho_u1u );

    *tau1 = rho_u1u;

    //bl1_zscals( delta, &rho_z1u );
    //bl1_zscals( delta, &rho_u1u );
    bl1_zneg1( &rho_z1u );
    bl1_zneg1( &rho_u1u );

    bl1_zaxpyv2b( m_U,
                  &rho_z1u,
                  &rho_u1u,
                  u1, rs_U,
                  z1, rs_Z,
                  w,  inc_w );
*/
/*
    bl1_zdotsv2( BLIS1_CONJUGATE,
                 m_U,
                 z1, rs_Z,
                 z2, rs_Z,
                 u,  inc_u,
                 &zero,
                 &rho_z1u,
                 &rho_z2u );
    bl1_zneg1( &rho_z1u );
    bl1_zneg1( &rho_z2u );

    bl1_zdotv2axpyv2b( m_U,
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

    bl1_zneg1( &rho_u1u );
    bl1_zneg1( &rho_u2u );

    bl1_zaxpyv2b( m_U,
                  &rho_u1u,
                  &rho_u2u,
                  z1, rs_Z,
                  z2, rs_Z,
                  w,  inc_w );
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

    /*------------------------------------------------------------*/

    u1   += step_u;
    u2   += step_u;
    z1   += step_z;
    z2   += step_z;
    tau1 += step_tau;
    tau2 += step_tau;
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

