/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opt_var1( FLA_Obj alpha, FLA_Obj tau, FLA_Obj u, FLA_Obj y, FLA_Obj z, FLA_Obj v, FLA_Obj A, FLA_Obj up, FLA_Obj a, FLA_Obj w )
{
/*
   Effective computation:
   A = A + alpha * ( u * y' + z * v' );
   y = A' * up;
   a = a - conj(y) / tau;
   w = A * conj(a);
*/
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_u, inc_y, inc_z, inc_v;
  int          inc_up, inc_a, inc_w;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_u    = FLA_Obj_vector_inc( u );
  inc_y    = FLA_Obj_vector_inc( y );
  inc_z    = FLA_Obj_vector_inc( z );
  inc_v    = FLA_Obj_vector_inc( v );

  inc_up   = FLA_Obj_vector_inc( up );
  inc_a    = FLA_Obj_vector_inc( a );
  inc_w    = FLA_Obj_vector_inc( w );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A   = FLA_FLOAT_PTR( A );
      float* buff_u   = FLA_FLOAT_PTR( u );
      float* buff_y   = FLA_FLOAT_PTR( y );
      float* buff_z   = FLA_FLOAT_PTR( z );
      float* buff_v   = FLA_FLOAT_PTR( v );
      float* buff_up  = FLA_FLOAT_PTR( up );
      float* buff_a   = FLA_FLOAT_PTR( a );
      float* buff_w   = FLA_FLOAT_PTR( w );
      float* buff_tau = FLA_FLOAT_PTR( tau );
      float* buff_alpha = FLA_FLOAT_PTR( alpha );

      FLA_Fused_Gerc2_Ahx_Axpy_Ax_ops_var1( m_A,
                                            n_A,
                                            buff_tau,
                                            buff_alpha,
                                            buff_u, inc_u,
                                            buff_y, inc_y,
                                            buff_z, inc_z,
                                            buff_v, inc_v,
                                            buff_A, rs_A, cs_A,
                                            buff_up, inc_up,
                                            buff_a, inc_a,
                                            buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A   = FLA_DOUBLE_PTR( A );
      double* buff_u   = FLA_DOUBLE_PTR( u );
      double* buff_y   = FLA_DOUBLE_PTR( y );
      double* buff_z   = FLA_DOUBLE_PTR( z );
      double* buff_v   = FLA_DOUBLE_PTR( v );
      double* buff_up  = FLA_DOUBLE_PTR( up );
      double* buff_a   = FLA_DOUBLE_PTR( a );
      double* buff_w   = FLA_DOUBLE_PTR( w );
      double* buff_tau = FLA_DOUBLE_PTR( tau );
      double* buff_alpha = FLA_DOUBLE_PTR( alpha );

      FLA_Fused_Gerc2_Ahx_Axpy_Ax_opd_var1( m_A,
                                            n_A,
                                            buff_tau,
                                            buff_alpha,
                                            buff_u, inc_u,
                                            buff_y, inc_y,
                                            buff_z, inc_z,
                                            buff_v, inc_v,
                                            buff_A, rs_A, cs_A,
                                            buff_up, inc_up,
                                            buff_a, inc_a,
                                            buff_w, inc_w );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A   = FLA_COMPLEX_PTR( A );
      scomplex* buff_u   = FLA_COMPLEX_PTR( u );
      scomplex* buff_y   = FLA_COMPLEX_PTR( y );
      scomplex* buff_z   = FLA_COMPLEX_PTR( z );
      scomplex* buff_v   = FLA_COMPLEX_PTR( v );
      scomplex* buff_up  = FLA_COMPLEX_PTR( up );
      scomplex* buff_a   = FLA_COMPLEX_PTR( a );
      scomplex* buff_w   = FLA_COMPLEX_PTR( w );
      scomplex* buff_tau = FLA_COMPLEX_PTR( tau );
      scomplex* buff_alpha = FLA_COMPLEX_PTR( alpha );

      FLA_Fused_Gerc2_Ahx_Axpy_Ax_opc_var1( m_A,
                                            n_A,
                                            buff_tau,
                                            buff_alpha,
                                            buff_u, inc_u,
                                            buff_y, inc_y,
                                            buff_z, inc_z,
                                            buff_v, inc_v,
                                            buff_A, rs_A, cs_A,
                                            buff_up, inc_up,
                                            buff_a, inc_a,
                                            buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A   = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_u   = FLA_DOUBLE_COMPLEX_PTR( u );
      dcomplex* buff_y   = FLA_DOUBLE_COMPLEX_PTR( y );
      dcomplex* buff_z   = FLA_DOUBLE_COMPLEX_PTR( z );
      dcomplex* buff_v   = FLA_DOUBLE_COMPLEX_PTR( v );
      dcomplex* buff_up  = FLA_DOUBLE_COMPLEX_PTR( up );
      dcomplex* buff_a   = FLA_DOUBLE_COMPLEX_PTR( a );
      dcomplex* buff_w   = FLA_DOUBLE_COMPLEX_PTR( w );
      dcomplex* buff_tau = FLA_DOUBLE_COMPLEX_PTR( tau );
      dcomplex* buff_alpha = FLA_DOUBLE_COMPLEX_PTR( alpha );

      FLA_Fused_Gerc2_Ahx_Axpy_Ax_opz_var1( m_A,
                                            n_A,
                                            buff_tau,
                                            buff_alpha,
                                            buff_u, inc_u,
                                            buff_y, inc_y,
                                            buff_z, inc_z,
                                            buff_v, inc_v,
                                            buff_A, rs_A, cs_A,
                                            buff_up, inc_up,
                                            buff_a, inc_a,
                                            buff_w, inc_w );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_ops_var1( int m_A,
                                                int n_A,
                                                float* buff_tau, 
                                                float* buff_alpha, 
                                                float* buff_u, int inc_u, 
                                                float* buff_y, int inc_y, 
                                                float* buff_z, int inc_z, 
                                                float* buff_v, int inc_v, 
                                                float* buff_A, int rs_A, int cs_A, 
                                                float* buff_up, int inc_up, 
                                                float* buff_a, int inc_a, 
                                                float* buff_w, int inc_w )
{
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  float     minus_inv_tau;
  int       i;

  bl1_ssetv( m_A,
             buff_0,
             buff_w, inc_w );

  minus_inv_tau = *buff_m1 / *buff_tau;

  for ( i = 0; i < n_A; ++i )
  {
    float*    a1       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    u        = buff_u;
    float*    psi1     = buff_y + (i  )*inc_y;
    float*    nu1      = buff_v + (i  )*inc_v;
    float*    z        = buff_z;
    float*    up       = buff_up;
    float*    alpha1   = buff_a + (i  )*inc_a;
    float*    w        = buff_w;
    float*    alpha    = buff_alpha;
    float     temp1;
    float     temp2;

    /*------------------------------------------------------------*/

    // bl1_smult3( alpha, psi1, &temp1 );
    temp1 = *alpha * *psi1;

    // bl1_smult3( alpha, nu1, &temp2 );
    temp2 = *alpha * *nu1;

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp1,
                u,  inc_u,
                a1, rs_A );
    //F77_saxpy( &m_A,
    //           &temp1,
    //           u,  &inc_u,
    //           a1, &rs_A );

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp2,
                z,  inc_z,
                a1, rs_A );
    //F77_saxpy( &m_A,
    //           &temp2,
    //           z,  &inc_z,
    //           a1, &rs_A );

    bl1_sdot( BLIS1_CONJUGATE,
              m_A,
              a1, rs_A,
              up,  inc_up,
              psi1 );
    //*psi1 = F77_sdot( &m_A,
    //                  a1, &rs_A,
    //                  up, &inc_up );

    // bl1_smult4( &minus_inv_tau, psi1, alpha1, alpha1 );
    *alpha1 = *alpha1 + minus_inv_tau * *psi1;

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                alpha1,
                a1, rs_A,
                w,  inc_w );
    //F77_saxpy( &m_A,
    //           alpha1,
    //           a1, &rs_A,
    //           w,  &inc_w );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opd_var1( int m_A,
                                                int n_A,
                                                double* buff_tau, 
                                                double* buff_alpha, 
                                                double* buff_u, int inc_u, 
                                                double* buff_y, int inc_y, 
                                                double* buff_z, int inc_z, 
                                                double* buff_v, int inc_v, 
                                                double* buff_A, int rs_A, int cs_A, 
                                                double* buff_up, int inc_up, 
                                                double* buff_a, int inc_a, 
                                                double* buff_w, int inc_w )
{
  double    zero      = bl1_d0();
  double    minus_one = bl1_dm1();
  double*   restrict u  = buff_u;
  double*   restrict up = buff_up;
  double*   restrict w = buff_w;
  double*   restrict z  = buff_z;
  double*   restrict alpha = buff_alpha;
  double*   restrict a1;
  double*   restrict a2;
  double*   restrict psi1;
  double*   restrict psi2;
  double*   restrict alpha1;
  double*   restrict alpha2;
  double*   restrict nu1;
  double*   restrict nu2;

  double    minus_inv_tau;
  double    alpha_conj_psi1;
  double    alpha_conj_psi2;
  double    alpha_conj_nu1;
  double    alpha_conj_nu2;
  int       i;
  int       n_run    = n_A / 2;
  int       n_left   = n_A % 2;
  int       twocs_A  = 2*cs_A;
  int       twoinc_y = 2*inc_y;
  int       twoinc_a = 2*inc_a;
  int       twoinc_v = 2*inc_v;


  bl1_dsetv( m_A,
             &zero,
             buff_w, inc_w );

  bl1_ddiv3( &minus_one, buff_tau, &minus_inv_tau );

  a1     = buff_A;
  a2     = buff_A + cs_A;
  psi1   = buff_y;
  psi2   = buff_y + inc_y;
  alpha1 = buff_a;
  alpha2 = buff_a + inc_a;
  nu1    = buff_v;
  nu2    = buff_v + inc_v;

  for ( i = 0; i < n_run; ++i )
  {

    /*------------------------------------------------------------*/

    bl1_dmult3( alpha, psi1, &alpha_conj_psi1 );
    bl1_dmult3( alpha, psi2, &alpha_conj_psi2 );

    bl1_dmult3( alpha, nu1, &alpha_conj_nu1 );
    bl1_dmult3( alpha, nu2, &alpha_conj_nu2 );

/*
   Effective computation:
   A = A + alpha * ( u * y' + z * v' );
   y = A' * up;
   a = a - conj(y) / tau;
   w = A * conj(a);
*/
    bl1_daxpyv2b( m_A,
	              &alpha_conj_psi1,
	              &alpha_conj_nu1,
	              u,  inc_u,
	              z,  inc_z,
	              a1, rs_A );
    bl1_daxpyv2b( m_A,
	              &alpha_conj_psi2,
	              &alpha_conj_nu2,
	              u,  inc_u,
	              z,  inc_z,
	              a2, rs_A );


    bl1_ddotsv2( BLIS1_CONJUGATE,
                 m_A,
                 a1, rs_A,
                 a2, rs_A,
                 up, inc_up,
                 &zero,
                 psi1,
                 psi2 );

    bl1_dmult4( &minus_inv_tau, psi1, alpha1, alpha1 );
    bl1_dmult4( &minus_inv_tau, psi2, alpha2, alpha2 );

    bl1_daxpyv2b( m_A,
                  alpha1,
                  alpha2,
                  a1, rs_A,
                  a2, rs_A,
                  w,  inc_w );

    /*------------------------------------------------------------*/

    a1     += twocs_A;
    a2     += twocs_A;
    psi1   += twoinc_y;
    psi2   += twoinc_y;
    alpha1 += twoinc_a;
    alpha2 += twoinc_a;
    nu1    += twoinc_v;
    nu2    += twoinc_v;
  }

  if ( n_left == 1 )
  {
    double   rho1;

    bl1_dmult3( alpha, psi1, &alpha_conj_psi1 );
    bl1_dmult3( alpha, nu1,  &alpha_conj_nu1 );

    bl1_daxpyv2b( m_A,
	              &alpha_conj_psi1,
	              &alpha_conj_nu1,
	              u,  inc_u,
	              z,  inc_z,
	              a1, rs_A );

    bl1_ddot( BLIS1_CONJUGATE,
              m_A,
              a1, rs_A,
              up, inc_up,
              &rho1 );
    bl1_dscals( &zero, psi1 );
    bl1_dadd3( psi1, &rho1, psi1 );

    bl1_dmult4( &minus_inv_tau, psi1, alpha1, alpha1 );

    bl1_daxpyv( BLIS1_NO_CONJUGATE,
	            m_A,
                alpha1,
                a1, rs_A,
                w,  inc_w );
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opc_var1( int m_A,
                                                int n_A,
                                                scomplex* buff_tau, 
                                                scomplex* buff_alpha, 
                                                scomplex* buff_u, int inc_u, 
                                                scomplex* buff_y, int inc_y, 
                                                scomplex* buff_z, int inc_z, 
                                                scomplex* buff_v, int inc_v, 
                                                scomplex* buff_A, int rs_A, int cs_A, 
                                                scomplex* buff_up, int inc_up, 
                                                scomplex* buff_a, int inc_a, 
                                                scomplex* buff_w, int inc_w )
{
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  scomplex  minus_inv_tau;
  scomplex  conj_psi1;
  scomplex  conj_nu1;
  scomplex  conj_alpha1;
  int       i;

  bl1_csetv( m_A,
             buff_0,
             buff_w, inc_w );

  bl1_cdiv3( buff_m1, buff_tau, &minus_inv_tau );

  for ( i = 0; i < n_A; ++i )
  {
    scomplex* a1       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* u        = buff_u;
    scomplex* psi1     = buff_y + (i  )*inc_y;
    scomplex* nu1      = buff_v + (i  )*inc_v;
    scomplex* z        = buff_z;
    scomplex* up       = buff_up;
    scomplex* alpha1   = buff_a + (i  )*inc_a;
    scomplex* w        = buff_w;
    scomplex* alpha    = buff_alpha;
    scomplex  temp1;
    scomplex  temp2;

    /*------------------------------------------------------------*/

    bl1_ccopyconj( psi1, &conj_psi1 );
    bl1_cmult3( alpha, &conj_psi1, &temp1 );

    bl1_ccopyconj( nu1, &conj_nu1 );
    bl1_cmult3( alpha, &conj_nu1, &temp2 );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp1,
                u,  inc_u,
                a1, rs_A );
    //F77_caxpy( &m_A,
    //           &temp1,
    //           u,  &inc_u,
    //           a1, &rs_A );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp2,
                z,  inc_z,
                a1, rs_A );
    //F77_caxpy( &m_A,
    //           &temp2,
    //           z,  &inc_z,
    //           a1, &rs_A );

    bl1_cdot( BLIS1_CONJUGATE,
              m_A,
              a1, rs_A,
              up, inc_up,
              psi1 );

    bl1_ccopyconj( psi1, &conj_psi1 );
    bl1_cmult4( &minus_inv_tau, &conj_psi1, alpha1, alpha1 );

    bl1_ccopyconj( alpha1, &conj_alpha1 );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &conj_alpha1,
                a1, rs_A,
                w,  inc_w );
    //F77_caxpy( &m_A,
    //           &conj_alpha1,
    //           a1, &rs_A,
    //           w,  &inc_w );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opz_var1( int m_A,
                                                int n_A,
                                                dcomplex* buff_tau, 
                                                dcomplex* buff_alpha, 
                                                dcomplex* buff_u, int inc_u, 
                                                dcomplex* buff_y, int inc_y, 
                                                dcomplex* buff_z, int inc_z, 
                                                dcomplex* buff_v, int inc_v, 
                                                dcomplex* buff_A, int rs_A, int cs_A, 
                                                dcomplex* buff_up, int inc_up, 
                                                dcomplex* buff_a, int inc_a, 
                                                dcomplex* buff_w, int inc_w )
{
  dcomplex  zero      = bl1_z0();
  dcomplex  minus_one = bl1_zm1();
  dcomplex* restrict u  = buff_u;
  dcomplex* restrict up = buff_up;
  dcomplex* restrict w = buff_w;
  dcomplex* restrict z  = buff_z;
  dcomplex* restrict alpha = buff_alpha;
  dcomplex* restrict a1;
  dcomplex* restrict a2;
  dcomplex* restrict psi1;
  dcomplex* restrict psi2;
  dcomplex* restrict alpha1;
  dcomplex* restrict alpha2;
  dcomplex* restrict nu1;
  dcomplex* restrict nu2;

  dcomplex  minus_inv_tau;
  dcomplex  conj_psi1;
  dcomplex  conj_psi2;
  dcomplex  conj_nu1;
  dcomplex  conj_nu2;
  dcomplex  conj_alpha1;
  dcomplex  conj_alpha2;
  dcomplex  alpha_conj_psi1;
  dcomplex  alpha_conj_psi2;
  dcomplex  alpha_conj_nu1;
  dcomplex  alpha_conj_nu2;
  int       i;
  int       n_run    = n_A / 2;
  int       n_left   = n_A % 2;
  int       twocs_A  = 2*cs_A;
  int       twoinc_y = 2*inc_y;
  int       twoinc_a = 2*inc_a;
  int       twoinc_v = 2*inc_v;


  bl1_zsetv( m_A,
             &zero,
             buff_w, inc_w );

  bl1_zdiv3( &minus_one, buff_tau, &minus_inv_tau );

  a1     = buff_A;
  a2     = buff_A + cs_A;
  psi1   = buff_y;
  psi2   = buff_y + inc_y;
  alpha1 = buff_a;
  alpha2 = buff_a + inc_a;
  nu1    = buff_v;
  nu2    = buff_v + inc_v;

  for ( i = 0; i < n_run; ++i )
  {

    /*------------------------------------------------------------*/

    bl1_zcopyconj( psi1, &conj_psi1 );
    bl1_zcopyconj( psi2, &conj_psi2 );
    bl1_zmult3( alpha, &conj_psi1, &alpha_conj_psi1 );
    bl1_zmult3( alpha, &conj_psi2, &alpha_conj_psi2 );

    bl1_zcopyconj( nu1, &conj_nu1 );
    bl1_zcopyconj( nu2, &conj_nu2 );
    bl1_zmult3( alpha, &conj_nu1, &alpha_conj_nu1 );
    bl1_zmult3( alpha, &conj_nu2, &alpha_conj_nu2 );

    bl1_zaxpyv2b( m_A,
	              &alpha_conj_psi1,
	              &alpha_conj_nu1,
	              u,  inc_u,
	              z,  inc_z,
	              a1, rs_A );
    bl1_zaxpyv2b( m_A,
	              &alpha_conj_psi2,
	              &alpha_conj_nu2,
	              u,  inc_u,
	              z,  inc_z,
	              a2, rs_A );


    bl1_zdotsv2( BLIS1_CONJUGATE,
                 m_A,
                 a1, rs_A,
                 a2, rs_A,
                 up, inc_up,
                 &zero,
                 psi1,
                 psi2 );

    bl1_zcopyconj( psi1, &conj_psi1 );
    bl1_zcopyconj( psi2, &conj_psi2 );
    bl1_zmult4( &minus_inv_tau, &conj_psi1, alpha1, alpha1 );
    bl1_zmult4( &minus_inv_tau, &conj_psi2, alpha2, alpha2 );
    bl1_zcopyconj( alpha1, &conj_alpha1 );
    bl1_zcopyconj( alpha2, &conj_alpha2 );

    bl1_zaxpyv2b( m_A,
                  &conj_alpha1,
                  &conj_alpha2,
                  a1, rs_A,
                  a2, rs_A,
                  w,  inc_w );

    /*------------------------------------------------------------*/

    a1     += twocs_A;
    a2     += twocs_A;
    psi1   += twoinc_y;
    psi2   += twoinc_y;
    alpha1 += twoinc_a;
    alpha2 += twoinc_a;
    nu1    += twoinc_v;
    nu2    += twoinc_v;
  }

  if ( n_left == 1 )
  {
    dcomplex rho1;

    bl1_zcopyconj( psi1, &conj_psi1 );
    bl1_zmult3( alpha, &conj_psi1, &alpha_conj_psi1 );
    bl1_zcopyconj( nu1, &conj_nu1 );
    bl1_zmult3( alpha, &conj_nu1, &alpha_conj_nu1 );

    bl1_zaxpyv2b( m_A,
	              &alpha_conj_psi1,
	              &alpha_conj_nu1,
	              u,  inc_u,
	              z,  inc_z,
	              a1, rs_A );

    bl1_zdot( BLIS1_CONJUGATE,
              m_A,
              a1, rs_A,
              up, inc_up,
              &rho1 );
    bl1_zscals( &zero, psi1 );
    bl1_zadd3( psi1, &rho1, psi1 );

    bl1_zcopyconj( psi1, &conj_psi1 );
    bl1_zmult4( &minus_inv_tau, &conj_psi1, alpha1, alpha1 );
    bl1_zcopyconj( alpha1, &conj_alpha1 );

    bl1_zaxpyv( BLIS1_NO_CONJUGATE,
	            m_A,
                &conj_alpha1,
                a1, rs_A,
                w,  inc_w );
  }

  return FLA_SUCCESS;
}

