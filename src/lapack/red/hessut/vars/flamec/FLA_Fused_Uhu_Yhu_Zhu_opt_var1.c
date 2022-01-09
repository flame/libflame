/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opt_var1( FLA_Obj delta, FLA_Obj U, FLA_Obj Y, FLA_Obj Z, FLA_Obj t, FLA_Obj u, FLA_Obj y, FLA_Obj z )
{
/*
   Effective computation:
   y = y + delta * ( Y ( U' u ) + U ( Z' u ) );
   z = z + delta * ( U ( Y' u ) + Z ( U' u ) );
   t = U' u;
*/
  FLA_Datatype datatype;
  integer          m_U, n_U;
  integer          rs_U, cs_U;
  integer          rs_Y, cs_Y;
  integer          rs_Z, cs_Z;
  integer          inc_u, inc_y, inc_z, inc_t;

  datatype = FLA_Obj_datatype( U );

  m_U      = FLA_Obj_length( U );
  n_U      = FLA_Obj_width( U );

  rs_U     = FLA_Obj_row_stride( U );
  cs_U     = FLA_Obj_col_stride( U );

  rs_Y     = FLA_Obj_row_stride( Y );
  cs_Y     = FLA_Obj_col_stride( Y );

  rs_Z     = FLA_Obj_row_stride( Z );
  cs_Z     = FLA_Obj_col_stride( Z );

  inc_u    = FLA_Obj_vector_inc( u );
  inc_y    = FLA_Obj_vector_inc( y );
  inc_z    = FLA_Obj_vector_inc( z );
  inc_t    = FLA_Obj_vector_inc( t );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float*    buff_U     = FLA_FLOAT_PTR( U );
      float*    buff_Y     = FLA_FLOAT_PTR( Y );
      float*    buff_Z     = FLA_FLOAT_PTR( Z );
      float*    buff_t     = FLA_FLOAT_PTR( t );
      float*    buff_u     = FLA_FLOAT_PTR( u );
      float*    buff_y     = FLA_FLOAT_PTR( y );
      float*    buff_z     = FLA_FLOAT_PTR( z );
      float*    buff_delta = FLA_FLOAT_PTR( delta );

      FLA_Fused_Uhu_Yhu_Zhu_ops_var1( m_U,
                                      n_U,
                                      buff_delta,
                                      buff_U, rs_U, cs_U,
                                      buff_Y, rs_Y, cs_Y,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_t, inc_t,
                                      buff_u, inc_u,
                                      buff_y, inc_y,
                                      buff_z, inc_z );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_U     = FLA_DOUBLE_PTR( U );
      double*   buff_Y     = FLA_DOUBLE_PTR( Y );
      double*   buff_Z     = FLA_DOUBLE_PTR( Z );
      double*   buff_t     = FLA_DOUBLE_PTR( t );
      double*   buff_u     = FLA_DOUBLE_PTR( u );
      double*   buff_y     = FLA_DOUBLE_PTR( y );
      double*   buff_z     = FLA_DOUBLE_PTR( z );
      double*   buff_delta = FLA_DOUBLE_PTR( delta );

      FLA_Fused_Uhu_Yhu_Zhu_opd_var1( m_U,
                                      n_U,
                                      buff_delta,
                                      buff_U, rs_U, cs_U,
                                      buff_Y, rs_Y, cs_Y,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_t, inc_t,
                                      buff_u, inc_u,
                                      buff_y, inc_y,
                                      buff_z, inc_z );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_U     = FLA_COMPLEX_PTR( U );
      scomplex* buff_Y     = FLA_COMPLEX_PTR( Y );
      scomplex* buff_Z     = FLA_COMPLEX_PTR( Z );
      scomplex* buff_t     = FLA_COMPLEX_PTR( t );
      scomplex* buff_u     = FLA_COMPLEX_PTR( u );
      scomplex* buff_y     = FLA_COMPLEX_PTR( y );
      scomplex* buff_z     = FLA_COMPLEX_PTR( z );
      scomplex* buff_delta = FLA_COMPLEX_PTR( delta );

      FLA_Fused_Uhu_Yhu_Zhu_opc_var1( m_U,
                                      n_U,
                                      buff_delta,
                                      buff_U, rs_U, cs_U,
                                      buff_Y, rs_Y, cs_Y,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_t, inc_t,
                                      buff_u, inc_u,
                                      buff_y, inc_y,
                                      buff_z, inc_z );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_U     = FLA_DOUBLE_COMPLEX_PTR( U );
      dcomplex* buff_Y     = FLA_DOUBLE_COMPLEX_PTR( Y );
      dcomplex* buff_Z     = FLA_DOUBLE_COMPLEX_PTR( Z );
      dcomplex* buff_t     = FLA_DOUBLE_COMPLEX_PTR( t );
      dcomplex* buff_u     = FLA_DOUBLE_COMPLEX_PTR( u );
      dcomplex* buff_y     = FLA_DOUBLE_COMPLEX_PTR( y );
      dcomplex* buff_z     = FLA_DOUBLE_COMPLEX_PTR( z );
      dcomplex* buff_delta = FLA_DOUBLE_COMPLEX_PTR( delta );

      FLA_Fused_Uhu_Yhu_Zhu_opz_var1( m_U,
                                      n_U,
                                      buff_delta,
                                      buff_U, rs_U, cs_U,
                                      buff_Y, rs_Y, cs_Y,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_t, inc_t,
                                      buff_u, inc_u,
                                      buff_y, inc_y,
                                      buff_z, inc_z );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Uhu_Yhu_Zhu_ops_var1( integer m_U,
                                          integer n_U,
                                          float* buff_delta, 
                                          float* buff_U, integer rs_U, integer cs_U, 
                                          float* buff_Y, integer rs_Y, integer cs_Y, 
                                          float* buff_Z, integer rs_Z, integer cs_Z, 
                                          float* buff_t, integer inc_t, 
                                          float* buff_u, integer inc_u, 
                                          float* buff_y, integer inc_y, 
                                          float* buff_z, integer inc_z ) 
{
  integer       i;

  for ( i = 0; i < n_U; ++i )
  {
    float*    u1       = buff_U + (i  )*cs_U + (0  )*rs_U;
    float*    y1       = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    float*    z1       = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    float*    delta    = buff_delta;
    float*    tau1     = buff_t + (i  )*inc_t;
    float*    u        = buff_u;
    float*    y        = buff_y;
    float*    z        = buff_z;
    float     alpha;
    float     beta;
    float     gamma;

    /*------------------------------------------------------------*/

    bl1_sdot( BLIS1_CONJUGATE,
              m_U,
              u1, rs_U,
              u,  inc_u,
              &alpha );
    //alpha = F77_sdot( &m_U,
    //                  u1, &rs_U,
    //                  u,  &inc_u );

    bl1_sdot( BLIS1_CONJUGATE,
              m_U,
              z1, rs_Z,
              u,  inc_u,
              &beta );
    //beta  = F77_sdot( &m_U,
    //                  z1, &rs_Z,
    //                  u,  &inc_u );

    bl1_sdot( BLIS1_CONJUGATE,
              m_U,
              y1, rs_Y,
              u,  inc_u,
              &gamma );
    //gamma = F77_sdot( &m_U,
    //                  y1, &rs_Y,
    //                  u,  &inc_u );

    *tau1 = alpha;

    // bl1_sscals( delta, &alpha );
    // bl1_sscals( delta, &beta );
    // bl1_sscals( delta, &gamma );
    alpha *= *delta;
    beta  *= *delta;
    gamma *= *delta;

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                y1, rs_Y,
                y,  inc_y );
    //F77_saxpy( &m_U,
    //           &alpha,
    //           y1, &rs_Y,
    //           y,  &inc_y );

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &beta,
                u1, rs_U,
                y,  inc_y );
    //F77_saxpy( &m_U,
    //           &beta,
    //           u1, &rs_U,
    //           y,  &inc_y );

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                z1, rs_Z,
                z,  inc_z );
    //F77_saxpy( &m_U,
    //           &alpha,
    //           z1, &rs_Z,
    //           z,  &inc_z );

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &gamma,
                u1, rs_U,
                z,  inc_z );
    //F77_saxpy( &m_U,
    //           &gamma,
    //           u1, &rs_U,
    //           z,  &inc_z );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opd_var1( integer m_U,
                                          integer n_U,
                                          double* buff_delta, 
                                          double* buff_U, integer rs_U, integer cs_U, 
                                          double* buff_Y, integer rs_Y, integer cs_Y, 
                                          double* buff_Z, integer rs_Z, integer cs_Z, 
                                          double* buff_t, integer inc_t, 
                                          double* buff_u, integer inc_u, 
                                          double* buff_y, integer inc_y, 
                                          double* buff_z, integer inc_z ) 
{
  double             zero  = bl1_d0();

  double*   restrict delta = buff_delta;
  double*   restrict u     = buff_u;
  double*   restrict y     = buff_y;
  double*   restrict z     = buff_z;

  double*   restrict u1;
  double*   restrict y1;
  double*   restrict z1;
  double*   restrict upsilon1;
  double*   restrict tau1;

  double    alpha;
  double    beta;
  double    gamma;

  integer       i;

  integer       n_run         = n_U / 1;
  //integer       n_left        = n_U % 1;
  integer       step_u1       = 1*cs_U;
  integer       step_y1       = 1*cs_Y;
  integer       step_z1       = 1*cs_Z;
  integer       step_upsilon1 = 1*inc_u;
  integer       step_tau1     = 1*inc_t;

  u1       = buff_U;
  y1       = buff_Y;
  z1       = buff_Z;
  upsilon1 = buff_u;
  tau1     = buff_t;

  for ( i = 0; i < n_run; ++i )
  {
    /*------------------------------------------------------------*/

/*
    bl1_ddotsv3( BLIS1_CONJUGATE,
                 m_U,
                 u1, rs_U,
                 z1, rs_Z,
                 y1, rs_Y,
                 u,  inc_u,
                 &zero,
                 &alpha,
                 &beta,
                 &gamma );

    *tau1 = alpha;

    bl1_dscals( delta, &alpha );
    bl1_dscals( delta, &beta );
    bl1_dscals( delta, &gamma );

    bl1_daxpyv2b( m_U,
                  &alpha,
                  &beta,
                  y1, rs_Y,
                  u1, rs_U,
                  y, inc_y );
    bl1_daxpyv2b( m_U,
                  &alpha,
                  &gamma,
                  z1, rs_Z,
                  u1, rs_U,
                  z, inc_z );
*/

    bl1_ddotsv2( BLIS1_CONJUGATE,
                 m_U,
                 y1, rs_Y,
                 z1, rs_Z,
                 u,  inc_u,
                 &zero,
                 &beta,
                 &gamma );

    bl1_ddotaxmyv2( m_U,
                    &gamma,
                    &beta,
                    u1, rs_U,
                    u,  inc_u,
                    &alpha,
                    y,  inc_y,
                    z,  inc_z );

    *tau1 = alpha;

    bl1_dscals( delta, &alpha );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                y1, rs_Y,
                y,  inc_y );
    bl1_daxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                z1, rs_Z,
                z,  inc_z );


    /*------------------------------------------------------------*/

    u1       += step_u1;
    y1       += step_y1;
    z1       += step_z1;
    upsilon1 += step_upsilon1;
    tau1     += step_tau1;
  }


  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opc_var1( integer m_U,
                                          integer n_U,
                                          scomplex* buff_delta, 
                                          scomplex* buff_U, integer rs_U, integer cs_U, 
                                          scomplex* buff_Y, integer rs_Y, integer cs_Y, 
                                          scomplex* buff_Z, integer rs_Z, integer cs_Z, 
                                          scomplex* buff_t, integer inc_t, 
                                          scomplex* buff_u, integer inc_u, 
                                          scomplex* buff_y, integer inc_y, 
                                          scomplex* buff_z, integer inc_z ) 
{
  integer       i;

  for ( i = 0; i < n_U; ++i )
  {
    scomplex* u1       = buff_U + (i  )*cs_U + (0  )*rs_U;
    scomplex* y1       = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    scomplex* z1       = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    scomplex* delta    = buff_delta;
    scomplex* tau1     = buff_t + (i  )*inc_t;
    scomplex* u        = buff_u;
    scomplex* y        = buff_y;
    scomplex* z        = buff_z;
    scomplex  alpha;
    scomplex  beta;
    scomplex  gamma;

    /*------------------------------------------------------------*/

    bl1_cdot( BLIS1_CONJUGATE,
              m_U,
              u1, rs_U,
              u,  inc_u,
              &alpha );

    bl1_cdot( BLIS1_CONJUGATE,
              m_U,
              z1, rs_Z,
              u,  inc_u,
              &beta );

    bl1_cdot( BLIS1_CONJUGATE,
              m_U,
              y1, rs_Y,
              u,  inc_u,
              &gamma );

    *tau1 = alpha;

    bl1_cscals( delta, &alpha );
    bl1_cscals( delta, &beta );
    bl1_cscals( delta, &gamma );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                y1, rs_Y,
                y,  inc_y );
    //F77_caxpy( &m_U,
    //           &alpha,
    //           y1, &rs_Y,
    //           y,  &inc_y );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &beta,
                u1, rs_U,
                y,  inc_y );
    //F77_caxpy( &m_U,
    //           &beta,
    //           u1, &rs_U,
    //           y,  &inc_y );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                z1, rs_Z,
                z,  inc_z );
    //F77_caxpy( &m_U,
    //           &alpha,
    //           z1, &rs_Z,
    //           z,  &inc_z );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &gamma,
                u1, rs_U,
                z,  inc_z );
    //F77_caxpy( &m_U,
    //           &gamma,
    //           u1, &rs_U,
    //           z,  &inc_z );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opz_var1( integer m_U,
                                          integer n_U,
                                          dcomplex* buff_delta, 
                                          dcomplex* buff_U, integer rs_U, integer cs_U, 
                                          dcomplex* buff_Y, integer rs_Y, integer cs_Y, 
                                          dcomplex* buff_Z, integer rs_Z, integer cs_Z, 
                                          dcomplex* buff_t, integer inc_t, 
                                          dcomplex* buff_u, integer inc_u, 
                                          dcomplex* buff_y, integer inc_y, 
                                          dcomplex* buff_z, integer inc_z ) 
{
  dcomplex           zero  = bl1_z0();

  dcomplex* restrict delta = buff_delta;
  dcomplex* restrict u     = buff_u;
  dcomplex* restrict y     = buff_y;
  dcomplex* restrict z     = buff_z;

  dcomplex* restrict u1;
  dcomplex* restrict y1;
  dcomplex* restrict z1;
  dcomplex* restrict upsilon1;
  dcomplex* restrict tau1;

  dcomplex  alpha;
  dcomplex  beta;
  dcomplex  gamma;

  integer       i;

  integer       n_run         = n_U / 1;
  //integer       n_left        = n_U % 1;
  integer       step_u1       = 1*cs_U;
  integer       step_y1       = 1*cs_Y;
  integer       step_z1       = 1*cs_Z;
  integer       step_upsilon1 = 1*inc_u;
  integer       step_tau1     = 1*inc_t;

  u1       = buff_U;
  y1       = buff_Y;
  z1       = buff_Z;
  upsilon1 = buff_u;
  tau1     = buff_t;

  for ( i = 0; i < n_run; ++i )
  {
    /*------------------------------------------------------------*/


    bl1_zdotsv3( BLIS1_CONJUGATE,
                 m_U,
                 u1, rs_U,
                 z1, rs_Z,
                 y1, rs_Y,
                 u,  inc_u,
                 &zero,
                 &alpha,
                 &beta,
                 &gamma );

    *tau1 = alpha;

    bl1_zscals( delta, &alpha );
    bl1_zscals( delta, &beta );
    bl1_zscals( delta, &gamma );

    bl1_zaxpyv2b( m_U,
                  &alpha,
                  &beta,
                  y1, rs_Y,
                  u1, rs_U,
                  y, inc_y );
    bl1_zaxpyv2b( m_U,
                  &alpha,
                  &gamma,
                  z1, rs_Z,
                  u1, rs_U,
                  z, inc_z );


/*
    bl1_zdotsv2( BLIS1_CONJUGATE,
                 m_U,
                 y1, rs_Y,
                 z1, rs_Z,
                 u,  inc_u,
                 &zero,
                 &beta,
                 &gamma );

    bl1_zdotaxmyv2( m_U,
                    &gamma,
                    &beta,
                    u1, rs_U,
                    u,  inc_u,
                    &alpha,
                    y,  inc_y,
                    z,  inc_z );

    *tau1 = alpha;

    bl1_zscals( delta, &alpha );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                y1, rs_Y,
                y,  inc_y );
    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
                m_U,
                &alpha,
                z1, rs_Z,
                z,  inc_z );
*/

    /*------------------------------------------------------------*/

    u1       += step_u1;
    y1       += step_y1;
    z1       += step_z1;
    upsilon1 += step_upsilon1;
    tau1     += step_tau1;
  }

  return FLA_SUCCESS;
}

