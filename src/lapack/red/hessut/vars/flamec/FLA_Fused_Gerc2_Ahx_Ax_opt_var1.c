/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fused_Gerc2_Ahx_Ax_opt_var1( FLA_Obj alpha, FLA_Obj u, FLA_Obj y, FLA_Obj z, FLA_Obj A, FLA_Obj x, FLA_Obj v, FLA_Obj w )
{
/*
   Effective computation:
   A = A + alpha * ( u * y' + z * u' );
   v = A' * x;
   w = A  * x;
*/
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_u, inc_y, inc_z, inc_x, inc_v, inc_w;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_u    = FLA_Obj_vector_inc( u );
  inc_y    = FLA_Obj_vector_inc( y );
  inc_z    = FLA_Obj_vector_inc( z );
  inc_x    = FLA_Obj_vector_inc( x );
  inc_v    = FLA_Obj_vector_inc( v );
  inc_w    = FLA_Obj_vector_inc( w );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_u = FLA_FLOAT_PTR( u );
      float* buff_y = FLA_FLOAT_PTR( y );
      float* buff_z = FLA_FLOAT_PTR( z );
      float* buff_x = FLA_FLOAT_PTR( x );
      float* buff_v = FLA_FLOAT_PTR( v );
      float* buff_w = FLA_FLOAT_PTR( w );
      float* buff_alpha = FLA_FLOAT_PTR( alpha );

      FLA_Fused_Gerc2_Ahx_Ax_ops_var1( m_A,
                                       n_A,
                                       buff_alpha,
                                       buff_u, inc_u,
                                       buff_y, inc_y,
                                       buff_z, inc_z,
                                       buff_A, rs_A, cs_A,
                                       buff_x, inc_x,
                                       buff_v, inc_v,
                                       buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_u = FLA_DOUBLE_PTR( u );
      double* buff_y = FLA_DOUBLE_PTR( y );
      double* buff_z = FLA_DOUBLE_PTR( z );
      double* buff_x = FLA_DOUBLE_PTR( x );
      double* buff_v = FLA_DOUBLE_PTR( v );
      double* buff_w = FLA_DOUBLE_PTR( w );
      double* buff_alpha = FLA_DOUBLE_PTR( alpha );

      FLA_Fused_Gerc2_Ahx_Ax_opd_var1( m_A,
                                       n_A,
                                       buff_alpha,
                                       buff_u, inc_u,
                                       buff_y, inc_y,
                                       buff_z, inc_z,
                                       buff_A, rs_A, cs_A,
                                       buff_x, inc_x,
                                       buff_v, inc_v,
                                       buff_w, inc_w );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_u = FLA_COMPLEX_PTR( u );
      scomplex* buff_y = FLA_COMPLEX_PTR( y );
      scomplex* buff_z = FLA_COMPLEX_PTR( z );
      scomplex* buff_x = FLA_COMPLEX_PTR( x );
      scomplex* buff_v = FLA_COMPLEX_PTR( v );
      scomplex* buff_w = FLA_COMPLEX_PTR( w );
      scomplex* buff_alpha = FLA_COMPLEX_PTR( alpha );

      FLA_Fused_Gerc2_Ahx_Ax_opc_var1( m_A,
                                       n_A,
                                       buff_alpha,
                                       buff_u, inc_u,
                                       buff_y, inc_y,
                                       buff_z, inc_z,
                                       buff_A, rs_A, cs_A,
                                       buff_x, inc_x,
                                       buff_v, inc_v,
                                       buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_u = FLA_DOUBLE_COMPLEX_PTR( u );
      dcomplex* buff_y = FLA_DOUBLE_COMPLEX_PTR( y );
      dcomplex* buff_z = FLA_DOUBLE_COMPLEX_PTR( z );
      dcomplex* buff_x = FLA_DOUBLE_COMPLEX_PTR( x );
      dcomplex* buff_v = FLA_DOUBLE_COMPLEX_PTR( v );
      dcomplex* buff_w = FLA_DOUBLE_COMPLEX_PTR( w );
      dcomplex* buff_alpha = FLA_DOUBLE_COMPLEX_PTR( alpha );

      FLA_Fused_Gerc2_Ahx_Ax_opz_var1( m_A,
                                       n_A,
                                       buff_alpha,
                                       buff_u, inc_u,
                                       buff_y, inc_y,
                                       buff_z, inc_z,
                                       buff_A, rs_A, cs_A,
                                       buff_x, inc_x,
                                       buff_v, inc_v,
                                       buff_w, inc_w );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Gerc2_Ahx_Ax_ops_var1( int m_A,
                                           int n_A,
                                           float* buff_alpha, 
                                           float* buff_u, int inc_u, 
                                           float* buff_y, int inc_y, 
                                           float* buff_z, int inc_z, 
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
    float*    psi1     = buff_y + (i  )*inc_y;
    float*    u        = buff_u;
    float*    upsilon1 = buff_u + (i  )*inc_u;
    float*    w        = buff_w;
    float*    z        = buff_z;
    float*    alpha    = buff_alpha;
    float     temp1;
    float     temp2;

    /*------------------------------------------------------------*/

    // bl1_scopyconj( psi1, &conj_psi1 );
    // bl1_smult3( alpha, &conj_psi1, &temp1 );
    temp1 = *alpha * *psi1;

    // bl1_scopyconj( upsilon1, &conj_upsilon1 );
    // bl1_smult3( alpha, &conj_upsilon1, &temp2 );
    temp2 = *alpha * *upsilon1;

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp1,
                u,  inc_u,
                a1, rs_A );
/*
    F77_saxpy( &m_A,
               &temp1,
               u,  &inc_u,
               a1, &rs_A );
*/

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp2,
                z,  inc_z,
                a1, rs_A );
/*
    F77_saxpy( &m_A,
               &temp2,
               z,  &inc_z,
               a1, &rs_A );
*/

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



FLA_Error FLA_Fused_Gerc2_Ahx_Ax_opd_var1( int m_A,
                                           int n_A,
                                           double* buff_alpha, 
                                           double* buff_u, int inc_u, 
                                           double* buff_y, int inc_y, 
                                           double* buff_z, int inc_z, 
                                           double* buff_A, int rs_A, int cs_A, 
                                           double* buff_x, int inc_x, 
                                           double* buff_v, int inc_v, 
                                           double* buff_w, int inc_w )
{
  double             zero  = bl1_d0();

  double*   restrict alpha = buff_alpha;
  double*   restrict u     = buff_u;
  double*   restrict z     = buff_z;
  double*   restrict x     = buff_x;
  double*   restrict w     = buff_w;

  double*   restrict a1;
  double*   restrict nu1;
  double*   restrict chi1;
  double*   restrict psi1;
  double*   restrict upsilon1;

  double             alpha_psi1;
  double             alpha_upsilon1;

  int       n_run         = n_A / 1;
  //int       n_left        = n_A % 1;
  int       step_a1       = 1*cs_A;
  int       step_nu1      = 1*inc_v;
  int       step_chi1     = 1*inc_x;
  int       step_psi1     = 1*inc_y;
  int       step_upsilon1 = 1*inc_u;
  int       i;

  bl1_dsetv( m_A,
             &zero,
             buff_w, inc_w );

  a1       = buff_A;
  nu1      = buff_v;
  chi1     = buff_x;
  psi1     = buff_y;
  upsilon1 = buff_u;

  for ( i = 0; i < n_run; ++i )
  {
    /*------------------------------------------------------------*/

    bl1_dmult3( alpha, psi1, &alpha_psi1 );
    bl1_dmult3( alpha, upsilon1, &alpha_upsilon1 );

    bl1_daxpyv2bdotaxpy( m_A,
                         &alpha_psi1,
                         u,  inc_u,
                         &alpha_upsilon1,
                         z,  inc_z,
                         a1, rs_A,
                         x,  inc_x,
                         chi1,
                         nu1,
                         w,  inc_w );

    /*------------------------------------------------------------*/

    a1       += step_a1;
    nu1      += step_nu1;
    chi1     += step_chi1;
    psi1     += step_psi1;
    upsilon1 += step_upsilon1;
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Gerc2_Ahx_Ax_opc_var1( int m_A,
                                           int n_A,
                                           scomplex* buff_alpha, 
                                           scomplex* buff_u, int inc_u, 
                                           scomplex* buff_y, int inc_y, 
                                           scomplex* buff_z, int inc_z, 
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
    scomplex* psi1     = buff_y + (i  )*inc_y;
    scomplex* u        = buff_u;
    scomplex* upsilon1 = buff_u + (i  )*inc_u;
    scomplex* w        = buff_w;
    scomplex* z        = buff_z;
    scomplex* alpha    = buff_alpha;
    scomplex  temp1;
    scomplex  temp2;
    scomplex  conj_psi1;
    scomplex  conj_upsilon1;

    /*------------------------------------------------------------*/

    bl1_ccopyconj( psi1, &conj_psi1 );
    bl1_cmult3( alpha, &conj_psi1, &temp1 );

    bl1_ccopyconj( upsilon1, &conj_upsilon1 );
    bl1_cmult3( alpha, &conj_upsilon1, &temp2 );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp1,
                u,  inc_u,
                a1, rs_A );
/*
    F77_caxpy( &m_A,
               &temp1,
               u,  &inc_u,
               a1, &rs_A );
*/

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp2,
                z,  inc_z,
                a1, rs_A );
/*
    F77_caxpy( &m_A,
               &temp2,
               z,  &inc_z,
               a1, &rs_A );
*/

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



FLA_Error FLA_Fused_Gerc2_Ahx_Ax_opz_var1( int m_A,
                                           int n_A,
                                           dcomplex* buff_alpha, 
                                           dcomplex* buff_u, int inc_u, 
                                           dcomplex* buff_y, int inc_y, 
                                           dcomplex* buff_z, int inc_z, 
                                           dcomplex* buff_A, int rs_A, int cs_A, 
                                           dcomplex* buff_x, int inc_x, 
                                           dcomplex* buff_v, int inc_v, 
                                           dcomplex* buff_w, int inc_w )
{
  dcomplex           zero  = bl1_z0();

  dcomplex* restrict alpha = buff_alpha;
  dcomplex* restrict u     = buff_u;
  dcomplex* restrict z     = buff_z;
  dcomplex* restrict x     = buff_x;
  dcomplex* restrict w     = buff_w;

  dcomplex* restrict a1;
  dcomplex* restrict nu1;
  dcomplex* restrict chi1;
  dcomplex* restrict psi1;
  dcomplex* restrict upsilon1;

  dcomplex           temp1;
  dcomplex           temp2;
  dcomplex           conj_psi1;
  dcomplex           conj_upsilon1;

  int       n_run         = n_A / 1;
  //int       n_left        = n_A % 1;
  int       step_a1       = 1*cs_A;
  int       step_nu1      = 1*inc_v;
  int       step_chi1     = 1*inc_x;
  int       step_psi1     = 1*inc_y;
  int       step_upsilon1 = 1*inc_u;
  int       i;

  bl1_zsetv( m_A,
             &zero,
             buff_w, inc_w );

  a1       = buff_A;
  nu1      = buff_v;
  chi1     = buff_x;
  psi1     = buff_y;
  upsilon1 = buff_u;

  for ( i = 0; i < n_run; ++i )
  {
    /*------------------------------------------------------------*/

    bl1_zcopyconj( psi1, &conj_psi1 );
    bl1_zmult3( alpha, &conj_psi1, &temp1 );

    bl1_zcopyconj( upsilon1, &conj_upsilon1 );
    bl1_zmult3( alpha, &conj_upsilon1, &temp2 );

/*
    bl1_zaxpyv2bdotaxpy( m_A,
                         &temp1,
                         u,  inc_u,
                         &temp2,
                         z,  inc_z,
                         a1, rs_A,
                         x,  inc_x,
                         chi1,
                         nu1,
                         w,  inc_w );
*/

    bl1_zaxpyv2b( m_A,
                  &temp1,
                  &temp2,
                  u,  inc_u,
                  z,  inc_z,
                  a1, rs_A );
    bl1_zdotaxpy( m_A,
                  a1, rs_A,
                  x,  inc_x,
                  chi1,
                  nu1,
                  w,  inc_w );

    /*------------------------------------------------------------*/

    a1       += step_a1;
    nu1      += step_nu1;
    chi1     += step_chi1;
    psi1     += step_psi1;
    upsilon1 += step_upsilon1;
  }

  return FLA_SUCCESS;
}

