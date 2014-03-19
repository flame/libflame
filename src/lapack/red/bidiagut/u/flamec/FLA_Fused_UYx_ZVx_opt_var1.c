/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fused_UYx_ZVx_opt_var1( FLA_Obj delta, FLA_Obj a, FLA_Obj U, FLA_Obj Y, FLA_Obj Z, FLA_Obj V, FLA_Obj A, FLA_Obj temp, FLA_Obj t, FLA_Obj w, FLA_Obj al )
{
/*
   Effective computation:
   w  = w      + delta * ( U ( Y' conj(a)  ) + Z ( V' conj(a)  ) );
   al = A * e0 + delta * ( U ( Y' e0       ) + Z ( V' e0       ) );
   t  = V' conj(a);
*/
  FLA_Datatype datatype;
  int          m_U, n_U;
  int          m_V, n_V;
  int          rs_A, cs_A;
  int          rs_U, cs_U;
  int          rs_Y, cs_Y;
  int          rs_Z, cs_Z;
  int          rs_V, cs_V;
  int          inc_a, inc_temp, inc_t, inc_w, inc_al;

  datatype = FLA_Obj_datatype( A );

  m_U      = FLA_Obj_length( U );
  n_U      = FLA_Obj_width( U );

  m_V      = FLA_Obj_length( V );
  n_V      = FLA_Obj_width( V );

  rs_U     = FLA_Obj_row_stride( U );
  cs_U     = FLA_Obj_col_stride( U );

  rs_Y     = FLA_Obj_row_stride( Y );
  cs_Y     = FLA_Obj_col_stride( Y );

  rs_Z     = FLA_Obj_row_stride( Z );
  cs_Z     = FLA_Obj_col_stride( Z );

  rs_V     = FLA_Obj_row_stride( V );
  cs_V     = FLA_Obj_col_stride( V );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_temp = FLA_Obj_vector_inc( temp );
  inc_t    = FLA_Obj_vector_inc( t );
  inc_a    = FLA_Obj_vector_inc( a );
  inc_w    = FLA_Obj_vector_inc( w );
  inc_al   = FLA_Obj_vector_inc( al );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float*    buff_A     = FLA_FLOAT_PTR( A );
      float*    buff_U     = FLA_FLOAT_PTR( U );
      float*    buff_Y     = FLA_FLOAT_PTR( Y );
      float*    buff_Z     = FLA_FLOAT_PTR( Z );
      float*    buff_V     = FLA_FLOAT_PTR( V );
      float*    buff_temp  = FLA_FLOAT_PTR( temp );
      float*    buff_t     = FLA_FLOAT_PTR( t );
      float*    buff_a     = FLA_FLOAT_PTR( a );
      float*    buff_w     = FLA_FLOAT_PTR( w );
      float*    buff_al    = FLA_FLOAT_PTR( al );
      float*    buff_delta = FLA_FLOAT_PTR( delta );

      FLA_Fused_UYx_ZVx_ops_var1( m_U,
                                  n_U,
                                  m_V,
                                  n_V,
                                  buff_delta,
                                  buff_U, rs_U, cs_U,
                                  buff_Y, rs_Y, cs_Y,
                                  buff_Z, rs_Z, cs_Z,
                                  buff_V, rs_V, cs_V,
                                  buff_A, rs_A, cs_A,
                                  buff_temp, inc_temp,
                                  buff_t, inc_t,
                                  buff_a, inc_a,
                                  buff_w, inc_w,
                                  buff_al, inc_al );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A     = FLA_DOUBLE_PTR( A );
      double*   buff_U     = FLA_DOUBLE_PTR( U );
      double*   buff_Y     = FLA_DOUBLE_PTR( Y );
      double*   buff_Z     = FLA_DOUBLE_PTR( Z );
      double*   buff_V     = FLA_DOUBLE_PTR( V );
      double*   buff_temp  = FLA_DOUBLE_PTR( temp );
      double*   buff_t     = FLA_DOUBLE_PTR( t );
      double*   buff_a     = FLA_DOUBLE_PTR( a );
      double*   buff_w     = FLA_DOUBLE_PTR( w );
      double*   buff_al    = FLA_DOUBLE_PTR( al );
      double*   buff_delta = FLA_DOUBLE_PTR( delta );

      FLA_Fused_UYx_ZVx_opd_var1( m_U,
                                  n_U,
                                  m_V,
                                  n_V,
                                  buff_delta,
                                  buff_U, rs_U, cs_U,
                                  buff_Y, rs_Y, cs_Y,
                                  buff_Z, rs_Z, cs_Z,
                                  buff_V, rs_V, cs_V,
                                  buff_A, rs_A, cs_A,
                                  buff_temp, inc_temp,
                                  buff_t, inc_t,
                                  buff_a, inc_a,
                                  buff_w, inc_w,
                                  buff_al, inc_al );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A     = FLA_COMPLEX_PTR( A );
      scomplex* buff_U     = FLA_COMPLEX_PTR( U );
      scomplex* buff_Y     = FLA_COMPLEX_PTR( Y );
      scomplex* buff_Z     = FLA_COMPLEX_PTR( Z );
      scomplex* buff_V     = FLA_COMPLEX_PTR( V );
      scomplex* buff_temp  = FLA_COMPLEX_PTR( temp );
      scomplex* buff_t     = FLA_COMPLEX_PTR( t );
      scomplex* buff_a     = FLA_COMPLEX_PTR( a );
      scomplex* buff_w     = FLA_COMPLEX_PTR( w );
      scomplex* buff_al    = FLA_COMPLEX_PTR( al );
      scomplex* buff_delta = FLA_COMPLEX_PTR( delta );

      FLA_Fused_UYx_ZVx_opc_var1( m_U,
                                  n_U,
                                  m_V,
                                  n_V,
                                  buff_delta,
                                  buff_U, rs_U, cs_U,
                                  buff_Y, rs_Y, cs_Y,
                                  buff_Z, rs_Z, cs_Z,
                                  buff_V, rs_V, cs_V,
                                  buff_A, rs_A, cs_A,
                                  buff_temp, inc_temp,
                                  buff_t, inc_t,
                                  buff_a, inc_a,
                                  buff_w, inc_w,
                                  buff_al, inc_al );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A     = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_U     = FLA_DOUBLE_COMPLEX_PTR( U );
      dcomplex* buff_Y     = FLA_DOUBLE_COMPLEX_PTR( Y );
      dcomplex* buff_Z     = FLA_DOUBLE_COMPLEX_PTR( Z );
      dcomplex* buff_V     = FLA_DOUBLE_COMPLEX_PTR( V );
      dcomplex* buff_temp  = FLA_DOUBLE_COMPLEX_PTR( temp );
      dcomplex* buff_t     = FLA_DOUBLE_COMPLEX_PTR( t );
      dcomplex* buff_a     = FLA_DOUBLE_COMPLEX_PTR( a );
      dcomplex* buff_w     = FLA_DOUBLE_COMPLEX_PTR( w );
      dcomplex* buff_al    = FLA_DOUBLE_COMPLEX_PTR( al );
      dcomplex* buff_delta = FLA_DOUBLE_COMPLEX_PTR( delta );

      FLA_Fused_UYx_ZVx_opz_var1( m_U,
                                  n_U,
                                  m_V,
                                  n_V,
                                  buff_delta,
                                  buff_U, rs_U, cs_U,
                                  buff_Y, rs_Y, cs_Y,
                                  buff_Z, rs_Z, cs_Z,
                                  buff_V, rs_V, cs_V,
                                  buff_A, rs_A, cs_A,
                                  buff_temp, inc_temp,
                                  buff_t, inc_t,
                                  buff_a, inc_a,
                                  buff_w, inc_w,
                                  buff_al, inc_al );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_UYx_ZVx_ops_var1( int m_U,
                                      int n_U,
                                      int m_V,
                                      int n_V,
                                      float* buff_delta, 
                                      float* buff_U, int rs_U, int cs_U, 
                                      float* buff_Y, int rs_Y, int cs_Y, 
                                      float* buff_Z, int rs_Z, int cs_Z, 
                                      float* buff_V, int rs_V, int cs_V, 
                                      float* buff_A, int rs_A, int cs_A, 
                                      float* buff_temp, int inc_temp, 
                                      float* buff_t, int inc_t, 
                                      float* buff_a, int inc_a, 
                                      float* buff_w, int inc_w, 
                                      float* buff_al, int inc_al ) 
{
  int       i;
  int       m_A = m_U;
  int       m_Z = m_U;

  bl1_scopyv( BLIS1_NO_CONJUGATE,
              m_A,
              buff_A,  rs_A,
              buff_al, inc_al );

  for ( i = 0; i < n_U; ++i )
  {
    float*    u1       = buff_U + (i  )*cs_U + (0  )*rs_U;
    float*    y1       = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    float*    z1       = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    float*    v1       = buff_V + (0  )*cs_V + (i  )*rs_V;
    float*    tau1     = buff_t + (i  )*inc_t;
    float*    delta    = buff_delta;
    float*    a        = buff_a;
    float*    w        = buff_w;
    float*    al       = buff_al;
    float*    psi20_l  = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    float*    nu20_l   = buff_V + (0  )*cs_V + (i  )*rs_V;
    float     alpha;
    float     beta;
    float     gamma;
    float     kappa;

    /*------------------------------------------------------------*/

    // bl1_sdot( BLIS1_NO_CONJUGATE,
    //           n_V,
    //           y1, rs_Y,
    //           a,  inc_a,
    //           &alpha );
    alpha = F77_sdot( &n_V,
                      y1, &rs_Y,
                      a,  &inc_a );

    // bl1_sdot( BLIS1_NO_CONJUGATE,
    //           n_V,
    //           v1, cs_V,
    //           a,  inc_a,
    //           &beta );
    beta = F77_sdot( &n_V,
                     v1, &cs_V,
                     a,  &inc_a );

    *tau1 = beta;

    // bl1_sconjs( &alpha );
    // bl1_sconjs( &beta );
    // bl1_scopyconj( psi20_l, &gamma );
    // bl1_scopyconj( nu20_l,  &kappa );
    gamma = *psi20_l;
    kappa = *nu20_l;

    // bl1_dscals( delta, &alpha );
    // bl1_dscals( delta, &beta );
    // bl1_dscals( delta, &gamma );
    // bl1_dscals( delta, &kappa );
    alpha *= *delta;
    beta  *= *delta;
    gamma *= *delta;
    kappa *= *delta;

    // bl1_saxpyv( BLIS1_NO_CONJUGATE,
    //             m_U,
    //             &alpha,
    //             u1, rs_U,
    //             w,  inc_w );
    F77_saxpy( &m_U,
               &alpha,
               u1, &rs_U,
               w,  &inc_w );

    // bl1_saxpyv( BLIS1_NO_CONJUGATE,
    //             m_Z,
    //             &beta,
    //             z1, rs_U,
    //             w,  inc_w );
    F77_saxpy( &m_Z,
               &beta,
               z1, &rs_Z,
               w,  &inc_w );

    // bl1_saxpyv( BLIS1_NO_CONJUGATE,
    //             m_U,
    //             &gamma,
    //             u1, rs_U,
    //             al, inc_al );
    F77_saxpy( &m_U,
               &gamma,
               u1, &rs_U,
               al, &inc_al );

    // bl1_saxpyv( BLIS1_NO_CONJUGATE,
    //             m_Z,
    //             &kappa,
    //             u1, rs_U,
    //             z,  inc_z );
    F77_saxpy( &m_Z,
               &kappa,
               z1, &rs_Z,
               al, &inc_al );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_UYx_ZVx_opd_var1( int m_U,
                                      int n_U,
                                      int m_V,
                                      int n_V,
                                      double* buff_delta, 
                                      double* buff_U, int rs_U, int cs_U, 
                                      double* buff_Y, int rs_Y, int cs_Y, 
                                      double* buff_Z, int rs_Z, int cs_Z, 
                                      double* buff_V, int rs_V, int cs_V, 
                                      double* buff_A, int rs_A, int cs_A, 
                                      double* buff_temp, int inc_temp, 
                                      double* buff_t, int inc_t, 
                                      double* buff_a, int inc_a, 
                                      double* buff_w, int inc_w, 
                                      double* buff_al, int inc_al ) 
{
  double    zero = bl1_d0();
  int       i;
  int       m_A = m_U;
  int       m_Z = m_U;

  bl1_dcopyv( BLIS1_NO_CONJUGATE,
              m_A,
              buff_A,  rs_A,
              buff_al, inc_al );

  if ( m_U == 0 || n_U == 0 ) return 0;
  if ( m_V == 0 || n_V == 0 ) return 0;

  for ( i = 0; i < n_U; ++i )
  {
    double*   restrict u1       = buff_U + (i  )*cs_U + (0  )*rs_U;
    double*   restrict y1       = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    double*   restrict z1       = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    double*   restrict v1       = buff_V + (0  )*cs_V + (i  )*rs_V;
    double*   restrict tau1     = buff_t + (i  )*inc_t;
    double*   restrict t1       = buff_temp;
    double*   restrict a        = buff_a;
    double*   restrict w        = buff_w;
    double*   restrict al       = buff_al;
    double*   restrict psi20_l  = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    double*   restrict nu20_l   = buff_V + (0  )*cs_V + (i  )*rs_V;
    double    alpha;
    double    beta;
    double    gamma;
    double    kappa;

    /*------------------------------------------------------------*/

    bl1_dcopyv( BLIS1_NO_CONJUGATE,
                n_V,
                v1, cs_V,
                t1, inc_t );

    bl1_ddotsv2( BLIS1_NO_CONJUGATE,
                 n_V,
                 y1, rs_Y,
                 t1, inc_t,
                 a,  inc_a,
                 &zero,
                 &alpha,
                 &beta );

    *tau1 = beta;

	bl1_dcopyconj( psi20_l, &gamma );
	bl1_dcopyconj( nu20_l,  &kappa );

	bl1_daxmyv2( BLIS1_NO_CONJUGATE,
                 m_U,
	             &alpha,
	             &gamma,
	             u1, rs_U,
	             w,  inc_w,
	             al, inc_al );

	bl1_daxmyv2( BLIS1_NO_CONJUGATE,
                 m_Z,
	             &beta,
	             &kappa,
	             z1, rs_U,
	             w,  inc_w,
	             al, inc_al );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_UYx_ZVx_opc_var1( int m_U,
                                      int n_U,
                                      int m_V,
                                      int n_V,
                                      scomplex* buff_delta, 
                                      scomplex* buff_U, int rs_U, int cs_U, 
                                      scomplex* buff_Y, int rs_Y, int cs_Y, 
                                      scomplex* buff_Z, int rs_Z, int cs_Z, 
                                      scomplex* buff_V, int rs_V, int cs_V, 
                                      scomplex* buff_A, int rs_A, int cs_A, 
                                      scomplex* buff_temp, int inc_temp, 
                                      scomplex* buff_t, int inc_t, 
                                      scomplex* buff_a, int inc_a, 
                                      scomplex* buff_w, int inc_w, 
                                      scomplex* buff_al, int inc_al ) 
{
  int       i;
  int       m_A = m_U;
  int       m_Z = m_U;

  bl1_ccopyv( BLIS1_NO_CONJUGATE,
              m_A,
              buff_A,  rs_A,
              buff_al, inc_al );

  for ( i = 0; i < n_U; ++i )
  {
    scomplex* u1       = buff_U + (i  )*cs_U + (0  )*rs_U;
    scomplex* y1       = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    scomplex* z1       = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    scomplex* v1       = buff_V + (0  )*cs_V + (i  )*rs_V;
    scomplex* tau1     = buff_t + (i  )*inc_t;
    scomplex* delta    = buff_delta;
    scomplex* a        = buff_a;
    scomplex* w        = buff_w;
    scomplex* al       = buff_al;
    scomplex* psi20_l  = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    scomplex* nu20_l   = buff_V + (0  )*cs_V + (i  )*rs_V;
    scomplex  alpha;
    scomplex  beta;
    scomplex  gamma;
    scomplex  kappa;

    /*------------------------------------------------------------*/

    bl1_cdot( BLIS1_NO_CONJUGATE,
              n_V,
              y1, rs_Y,
              a,  inc_a,
              &alpha );

    bl1_cdot( BLIS1_NO_CONJUGATE,
              n_V,
              v1, cs_V,
              a,  inc_a,
              &beta );

    bl1_cconjs( &alpha );
    bl1_cconjs( &beta );
	bl1_ccopyconj( psi20_l, &gamma );
	bl1_ccopyconj( nu20_l,  &kappa );

    *tau1 = beta;

    bl1_cscals( delta, &alpha );
    bl1_cscals( delta, &beta );
    bl1_cscals( delta, &gamma );
    bl1_cscals( delta, &kappa );

    // bl1_caxpyv( BLIS1_NO_CONJUGATE,
    //             m_U,
    //             &alpha,
    //             u1, rs_U,
    //             w,  inc_w );
    F77_caxpy( &m_U,
               &alpha,
               u1, &rs_U,
               w,  &inc_w );

    // bl1_caxpyv( BLIS1_NO_CONJUGATE,
    //             m_Z,
    //             &beta,
    //             z1, rs_U,
    //             w,  inc_w );
    F77_caxpy( &m_Z,
               &beta,
               z1, &rs_Z,
               w,  &inc_w );

    // bl1_caxpyv( BLIS1_NO_CONJUGATE,
    //             m_U,
    //             &gamma,
    //             u1, rs_U,
    //             al, inc_al );
    F77_caxpy( &m_U,
               &gamma,
               u1, &rs_U,
               al, &inc_al );

    // bl1_caxpyv( BLIS1_NO_CONJUGATE,
    //             m_Z,
    //             &kappa,
    //             u1, rs_U,
    //             z,  inc_z );
    F77_caxpy( &m_Z,
               &kappa,
               z1, &rs_Z,
               al, &inc_al );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_UYx_ZVx_opz_var1( int m_U,
                                      int n_U,
                                      int m_V,
                                      int n_V,
                                      dcomplex* buff_delta, 
                                      dcomplex* buff_U, int rs_U, int cs_U, 
                                      dcomplex* buff_Y, int rs_Y, int cs_Y, 
                                      dcomplex* buff_Z, int rs_Z, int cs_Z, 
                                      dcomplex* buff_V, int rs_V, int cs_V, 
                                      dcomplex* buff_A, int rs_A, int cs_A, 
                                      dcomplex* buff_temp, int inc_temp, 
                                      dcomplex* buff_t, int inc_t, 
                                      dcomplex* buff_a, int inc_a, 
                                      dcomplex* buff_w, int inc_w, 
                                      dcomplex* buff_al, int inc_al ) 
{
  dcomplex  zero = bl1_z0();
  int       i;
  int       m_A = m_U;
  int       m_Z = m_U;

  bl1_zcopyv( BLIS1_NO_CONJUGATE,
              m_A,
              buff_A,  rs_A,
              buff_al, inc_al );

  if ( m_U == 0 || n_U == 0 ) return 0;
  if ( m_V == 0 || n_V == 0 ) return 0;

  for ( i = 0; i < n_U; ++i )
  {
    dcomplex* restrict u1       = buff_U + (i  )*cs_U + (0  )*rs_U;
    dcomplex* restrict y1       = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    dcomplex* restrict z1       = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    dcomplex* restrict v1       = buff_V + (0  )*cs_V + (i  )*rs_V;
    dcomplex* restrict tau1     = buff_t + (i  )*inc_t;
    dcomplex* restrict a        = buff_a;
    dcomplex* restrict w        = buff_w;
    dcomplex* restrict al       = buff_al;
    dcomplex* restrict psi20_l  = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    dcomplex* restrict nu20_l   = buff_V + (0  )*cs_V + (i  )*rs_V;
    dcomplex  alpha;
    dcomplex  beta;
    dcomplex  gamma;
    dcomplex  kappa;

    /*------------------------------------------------------------*/

    bl1_zdotsv2( BLIS1_NO_CONJUGATE,
                 n_V,
                 y1, rs_Y,
                 v1, cs_V,
                 a,  inc_a,
                 &zero,
                 &alpha,
                 &beta );

    bl1_zconjs( &alpha );
    bl1_zconjs( &beta );

    *tau1 = beta;

	bl1_zcopyconj( psi20_l, &gamma );
	bl1_zcopyconj( nu20_l,  &kappa );

	bl1_zaxmyv2( BLIS1_NO_CONJUGATE,
                 m_U,
	             &alpha,
	             &gamma,
	             u1, rs_U,
	             w,  inc_w,
	             al, inc_al );

	bl1_zaxmyv2( BLIS1_NO_CONJUGATE,
                 m_Z,
	             &beta,
	             &kappa,
	             z1, rs_U,
	             w,  inc_w,
	             al, inc_al );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

