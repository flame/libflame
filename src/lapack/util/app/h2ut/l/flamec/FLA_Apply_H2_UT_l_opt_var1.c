/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_H2_UT_l_opt_var1( FLA_Obj tau, FLA_Obj u2, FLA_Obj a1t,
                                                               FLA_Obj A2 )
/*
  Compute:

    / a1t \  :=  / I - 1/tau / 1  \ ( 1  u2' ) \ / a1t \ 
    \ A2  /      \           \ u2 /            / \ A2  / 
 
  w = ( a1t + u2' * A2 ) / tau;

  a1t = a1t - w;
  A2  = A2  - u2 * w;
*/
{
  FLA_Datatype datatype;
  integer          m_u2_A2;
  integer          n_a1t;
  integer          inc_u2;
  integer          inc_a1t;
  integer          rs_A2;
  integer          cs_A2;

  // The house-holder transformation in libFLAME never creates a zero tau value.
  // However, when libFLAME is mixed with LAPACK, zero tau means to apply an 
  // identity matrix that does nothing here.
  if ( FLA_Obj_has_zero_dim( a1t ) ||
       FLA_Obj_equals( tau, FLA_ZERO ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A2 );

  m_u2_A2  = FLA_Obj_length( A2 );
  n_a1t    = FLA_Obj_width( a1t );
  inc_u2   = FLA_Obj_vector_inc( u2 );
  inc_a1t  = FLA_Obj_vector_inc( a1t );
  rs_A2    = FLA_Obj_row_stride( A2 );
  cs_A2    = FLA_Obj_col_stride( A2 );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* tau_p = ( float* ) FLA_FLOAT_PTR( tau );
      float* u2_p  = ( float* ) FLA_FLOAT_PTR( u2 );
      float* a1t_p = ( float* ) FLA_FLOAT_PTR( a1t );
      float* A2_p  = ( float* ) FLA_FLOAT_PTR( A2 );

      FLA_Apply_H2_UT_l_ops_var1( m_u2_A2, n_a1t,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1t_p, inc_a1t,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_DOUBLE:
    {
      double* tau_p = ( double* ) FLA_DOUBLE_PTR( tau );
      double* u2_p  = ( double* ) FLA_DOUBLE_PTR( u2 );
      double* a1t_p = ( double* ) FLA_DOUBLE_PTR( a1t );
      double* A2_p  = ( double* ) FLA_DOUBLE_PTR( A2 );

      FLA_Apply_H2_UT_l_opd_var1( m_u2_A2, n_a1t,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1t_p, inc_a1t,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* tau_p = ( scomplex* ) FLA_COMPLEX_PTR( tau );
      scomplex* u2_p  = ( scomplex* ) FLA_COMPLEX_PTR( u2 );
      scomplex* a1t_p = ( scomplex* ) FLA_COMPLEX_PTR( a1t );
      scomplex* A2_p  = ( scomplex* ) FLA_COMPLEX_PTR( A2 );

      FLA_Apply_H2_UT_l_opc_var1( m_u2_A2, n_a1t,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1t_p, inc_a1t,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* tau_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( tau );
      dcomplex* u2_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( u2 );
      dcomplex* a1t_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( a1t );
      dcomplex* A2_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A2 );

      FLA_Apply_H2_UT_l_opz_var1( m_u2_A2, n_a1t,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1t_p, inc_a1t,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_l_ops_var1( integer m_u2_A2,
                                      integer n_a1t,
                                      float* tau,
                                      float* u2, integer inc_u2,
                                      float* a1t, integer inc_a1t,
                                      float* A2, integer rs_A2, integer cs_A2 )
{
  float*    one_p       = FLA_FLOAT_PTR( FLA_ONE );
  float*    minus_one_p = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  integer       inc_w1t;

  // FLA_Obj w1t;
  float*    w1t;

  if ( n_a1t == 0 || *tau == 0.0F ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  w1t = ( float* ) FLA_malloc( n_a1t * sizeof( *a1t ) );
  inc_w1t = 1;

  // // w1t = a1t;
  // FLA_Copy_external( a1t, w1t );
  bl1_scopyv( BLIS1_NO_CONJUGATE,
              n_a1t,
              a1t, inc_a1t, 
              w1t, inc_w1t ); 

  // // w1t = w1t + u2' * A2;
  // // w1t = w1t + A2^T * conj(u2);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, FLA_ONE, w1t );
  bl1_sgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_u2_A2,
             n_a1t,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1t, inc_w1t );

  // // w1t = w1t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );
  bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                 n_a1t,
                 tau,
                 w1t, inc_w1t );

  // // a1t = - w1t + a1t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );
  bl1_saxpyv( BLIS1_NO_CONJUGATE,
              n_a1t,
              minus_one_p,
              w1t, inc_w1t,
              a1t, inc_a1t );

  // // A2 = - u2 * w1t + A2;
  // FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );
  bl1_sger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_u2_A2,
            n_a1t,
            minus_one_p,
            u2, inc_u2,
            w1t, inc_w1t,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1t );
  FLA_free( w1t );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_l_opd_var1( integer m_u2_A2,
                                      integer n_a1t,
                                      double* tau,
                                      double* u2, integer inc_u2,
                                      double* a1t, integer inc_a1t,
                                      double* A2, integer rs_A2, integer cs_A2 )
{
  double*   one_p       = FLA_DOUBLE_PTR( FLA_ONE );
  double*   minus_one_p = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  integer       inc_w1t;

  // FLA_Obj w1t;
  double*   w1t;

  if ( n_a1t == 0 && *tau == 0.0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  w1t = ( double* ) FLA_malloc( n_a1t * sizeof( *a1t ) );
  inc_w1t = 1;

  // // w1t = a1t;
  // FLA_Copy_external( a1t, w1t );
  bl1_dcopyv( BLIS1_NO_CONJUGATE,
              n_a1t,
              a1t, inc_a1t, 
              w1t, inc_w1t ); 

  // // w1t = w1t + u2' * A2;
  // // w1t = w1t + A2^T * conj(u2);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, FLA_ONE, w1t );
  bl1_dgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_u2_A2,
             n_a1t,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1t, inc_w1t );

  // // w1t = w1t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );
  bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                 n_a1t,
                 tau,
                 w1t, inc_w1t );

  // // a1t = - w1t + a1t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );
  bl1_daxpyv( BLIS1_NO_CONJUGATE,
              n_a1t,
              minus_one_p,
              w1t, inc_w1t,
              a1t, inc_a1t );

  // // A2 = - u2 * w1t + A2;
  // FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );
  bl1_dger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_u2_A2,
            n_a1t,
            minus_one_p,
            u2, inc_u2,
            w1t, inc_w1t,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1t );
  FLA_free( w1t );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_l_opc_var1( integer m_u2_A2,
                                      integer n_a1t,
                                      scomplex* tau,
                                      scomplex* u2, integer inc_u2,
                                      scomplex* a1t, integer inc_a1t,
                                      scomplex* A2, integer rs_A2, integer cs_A2 )
{
  scomplex* one_p       = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* minus_one_p = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  integer       inc_w1t;

  // FLA_Obj w1t;
  scomplex* w1t;

  if ( n_a1t == 0 || ( tau->real == 0.0F &&
                       tau->imag == 0.0F ) ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  w1t = ( scomplex* ) FLA_malloc( n_a1t * sizeof( *a1t ) );
  inc_w1t = 1;

  // // w1t = a1t;
  // FLA_Copy_external( a1t, w1t );
  bl1_ccopyv( BLIS1_NO_CONJUGATE,
              n_a1t,
              a1t, inc_a1t, 
              w1t, inc_w1t ); 

  // // w1t = w1t + u2' * A2;
  // // w1t = w1t + A2^T * conj(u2);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, FLA_ONE, w1t );
  bl1_cgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_u2_A2,
             n_a1t,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1t, inc_w1t );

  // // w1t = w1t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );
  bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                 n_a1t,
                 tau,
                 w1t, inc_w1t );

  // // a1t = - w1t + a1t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );
  bl1_caxpyv( BLIS1_NO_CONJUGATE,
              n_a1t,
              minus_one_p,
              w1t, inc_w1t,
              a1t, inc_a1t );

  // // A2 = - u2 * w1t + A2;
  // FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );
  bl1_cger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_u2_A2,
            n_a1t,
            minus_one_p,
            u2, inc_u2,
            w1t, inc_w1t,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1t );
  FLA_free( w1t );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_l_opz_var1( integer m_u2_A2,
                                      integer n_a1t,
                                      dcomplex* tau,
                                      dcomplex* u2, integer inc_u2,
                                      dcomplex* a1t, integer inc_a1t,
                                      dcomplex* A2, integer rs_A2, integer cs_A2 )
{
  dcomplex* one_p       = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* minus_one_p = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  integer       inc_w1t;

  // FLA_Obj w1t;
  dcomplex* w1t;

  if ( n_a1t == 0 || ( tau->real == 0.0 &&
                       tau->imag == 0.0 ) ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  w1t = ( dcomplex* ) FLA_malloc( n_a1t * sizeof( *a1t ) );
  inc_w1t = 1;

  // // w1t = a1t;
  // FLA_Copy_external( a1t, w1t );
  bl1_zcopyv( BLIS1_NO_CONJUGATE,
              n_a1t,
              a1t, inc_a1t, 
              w1t, inc_w1t ); 

  // // w1t = w1t + u2' * A2;
  // // w1t = w1t + A2^T * conj(u2);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, FLA_ONE, w1t );
  bl1_zgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_u2_A2,
             n_a1t,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1t, inc_w1t );

  // // w1t = w1t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );
  bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                 n_a1t,
                 tau,
                 w1t, inc_w1t );

  // // a1t = - w1t + a1t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );
  bl1_zaxpyv( BLIS1_NO_CONJUGATE,
              n_a1t,
              minus_one_p,
              w1t, inc_w1t,
              a1t, inc_a1t );

  // // A2 = - u2 * w1t + A2;
  // FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );
  bl1_zger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_u2_A2,
            n_a1t,
            minus_one_p,
            u2, inc_u2,
            w1t, inc_w1t,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1t );
  FLA_free( w1t );

  return FLA_SUCCESS;
}
