/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_HUD_UT_l_opt_var1( FLA_Obj tau, FLA_Obj w12t,
                                                    FLA_Obj r12t,
                                       FLA_Obj u1,  FLA_Obj C2,
                                       FLA_Obj v1,  FLA_Obj D2 )
{
  FLA_Datatype datatype;
  int          m_u1_C2;
  int          m_v1_D2;
  int          n_r12t;
  int          inc_u1;
  int          inc_v1;
  int          inc_w12t;
  int          inc_r12t;
  int          rs_C2, cs_C2;
  int          rs_D2, cs_D2;

  if ( FLA_Obj_has_zero_dim( r12t ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( C2 );

  m_u1_C2  = FLA_Obj_length( u1 );
  m_v1_D2  = FLA_Obj_length( v1 );
  n_r12t   = FLA_Obj_width( r12t );
  inc_w12t = FLA_Obj_vector_inc( w12t );
  inc_r12t = FLA_Obj_vector_inc( r12t );
  inc_u1   = FLA_Obj_vector_inc( u1 );
  rs_C2    = FLA_Obj_row_stride( C2 );
  cs_C2    = FLA_Obj_col_stride( C2 );
  inc_v1   = FLA_Obj_vector_inc( v1 );
  rs_D2    = FLA_Obj_row_stride( D2 );
  cs_D2    = FLA_Obj_col_stride( D2 );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* tau_p  = ( float* ) FLA_FLOAT_PTR( tau );
      float* w12t_p = ( float* ) FLA_FLOAT_PTR( w12t );
      float* r12t_p = ( float* ) FLA_FLOAT_PTR( r12t );
      float* u1_p   = ( float* ) FLA_FLOAT_PTR( u1 );
      float* C2_p   = ( float* ) FLA_FLOAT_PTR( C2 );
      float* v1_p   = ( float* ) FLA_FLOAT_PTR( v1 );
      float* D2_p   = ( float* ) FLA_FLOAT_PTR( D2 );

      FLA_Apply_HUD_UT_l_ops_var1( m_u1_C2,
                                   m_v1_D2,
                                   n_r12t,
                                   tau_p,
                                   w12t_p, inc_w12t,
                                   r12t_p, inc_r12t,
                                   u1_p, inc_u1,
                                   C2_p, rs_C2, cs_C2,
                                   v1_p, inc_v1,
                                   D2_p, rs_D2, cs_D2 );
      break;
    }

    case FLA_DOUBLE:
    {
      double* tau_p  = ( double* ) FLA_DOUBLE_PTR( tau );
      double* w12t_p = ( double* ) FLA_DOUBLE_PTR( w12t );
      double* r12t_p = ( double* ) FLA_DOUBLE_PTR( r12t );
      double* u1_p   = ( double* ) FLA_DOUBLE_PTR( u1 );
      double* C2_p   = ( double* ) FLA_DOUBLE_PTR( C2 );
      double* v1_p   = ( double* ) FLA_DOUBLE_PTR( v1 );
      double* D2_p   = ( double* ) FLA_DOUBLE_PTR( D2 );

      FLA_Apply_HUD_UT_l_opd_var1( m_u1_C2,
                                   m_v1_D2,
                                   n_r12t,
                                   tau_p,
                                   w12t_p, inc_w12t,
                                   r12t_p, inc_r12t,
                                   u1_p, inc_u1,
                                   C2_p, rs_C2, cs_C2,
                                   v1_p, inc_v1,
                                   D2_p, rs_D2, cs_D2 );
      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* tau_p  = ( scomplex* ) FLA_COMPLEX_PTR( tau );
      scomplex* w12t_p = ( scomplex* ) FLA_COMPLEX_PTR( w12t );
      scomplex* r12t_p = ( scomplex* ) FLA_COMPLEX_PTR( r12t );
      scomplex* u1_p   = ( scomplex* ) FLA_COMPLEX_PTR( u1 );
      scomplex* C2_p   = ( scomplex* ) FLA_COMPLEX_PTR( C2 );
      scomplex* v1_p   = ( scomplex* ) FLA_COMPLEX_PTR( v1 );
      scomplex* D2_p   = ( scomplex* ) FLA_COMPLEX_PTR( D2 );

      FLA_Apply_HUD_UT_l_opc_var1( m_u1_C2,
                                   m_v1_D2,
                                   n_r12t,
                                   tau_p,
                                   w12t_p, inc_w12t,
                                   r12t_p, inc_r12t,
                                   u1_p, inc_u1,
                                   C2_p, rs_C2, cs_C2,
                                   v1_p, inc_v1,
                                   D2_p, rs_D2, cs_D2 );
      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* tau_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( tau );
      dcomplex* w12t_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( w12t );
      dcomplex* r12t_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( r12t );
      dcomplex* u1_p   = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( u1 );
      dcomplex* C2_p   = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( C2 );
      dcomplex* v1_p   = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( v1 );
      dcomplex* D2_p   = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( D2 );

      FLA_Apply_HUD_UT_l_opz_var1( m_u1_C2,
                                   m_v1_D2,
                                   n_r12t,
                                   tau_p,
                                   w12t_p, inc_w12t,
                                   r12t_p, inc_r12t,
                                   u1_p, inc_u1,
                                   C2_p, rs_C2, cs_C2,
                                   v1_p, inc_v1,
                                   D2_p, rs_D2, cs_D2 );
      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_HUD_UT_l_ops_var1( int m_u1_C2,
                                       int m_v1_D2,
                                       int n_r12t,
                                       float* tau,
                                       float* w12t, int inc_w12t,
                                       float* r12t, int inc_r12t,
                                       float* u1, int inc_u1,
                                       float* C2, int rs_C2, int cs_C2,
                                       float* v1, int inc_v1,
                                       float* D2, int rs_D2, int cs_D2 )
{
  float*    one_p       = FLA_FLOAT_PTR( FLA_ONE );
  float*    minus_one_p = FLA_FLOAT_PTR( FLA_MINUS_ONE );

  // if ( FLA_Obj_has_zero_dim( r12t ) ) return FLA_SUCCESS;
  if ( n_r12t == 0 ) return FLA_SUCCESS;

  // // w12t = r12t;
  // FLA_Copy_external( r12t, w12t );
  bl1_scopyv( BLIS1_NO_CONJUGATE,
              n_r12t,
              r12t, inc_r12t, 
              w12t, inc_w12t ); 

  // // w12t = w12t + u1' * C2;
  // //      = w12t + C2^T * conj(u1);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, C2, u1, FLA_ONE, w12t );
  bl1_sgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_u1_C2,
             n_r12t,
             one_p,
             C2, rs_C2, cs_C2,
             u1, inc_u1,
             one_p,
             w12t, inc_w12t );

  // // w12t = w12t + v1' * D2;
  // //      = w12t + D2^T * conj(v1);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, D2, v1, FLA_ONE, w12t );
  bl1_sgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_v1_D2,
             n_r12t,
             one_p,
             D2, rs_D2, cs_D2,
             v1, inc_v1,
             one_p,
             w12t, inc_w12t );

  // // w12t = w12t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w12t );
  bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                 n_r12t,
                 tau,
                 w12t, inc_w12t );

  // // r12t = - w12t + r12t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w12t, r12t );
  bl1_saxpyv( BLIS1_NO_CONJUGATE,
              n_r12t,
              minus_one_p,
              w12t, inc_w12t,
              r12t, inc_r12t );

  // // C2 = - u1 * w12t + C2;
  // FLA_Ger_external( FLA_MINUS_ONE, u1, w12t, C2 );
  bl1_sger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_u1_C2,
            n_r12t,
            minus_one_p,
            u1, inc_u1,
            w12t, inc_w12t,
            C2, rs_C2, cs_C2 );

  // // D2 = v1 * w12t + D2;
  // FLA_Ger_external( FLA_ONE, v1, w12t, D2 );
  bl1_sger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_v1_D2,
            n_r12t,
            one_p,
            v1, inc_v1,
            w12t, inc_w12t,
            D2, rs_D2, cs_D2 );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_HUD_UT_l_opd_var1( int m_u1_C2,
                                       int m_v1_D2,
                                       int n_r12t,
                                       double* tau,
                                       double* w12t, int inc_w12t,
                                       double* r12t, int inc_r12t,
                                       double* u1, int inc_u1,
                                       double* C2, int rs_C2, int cs_C2,
                                       double* v1, int inc_v1,
                                       double* D2, int rs_D2, int cs_D2 )
{
  double*   one_p       = FLA_DOUBLE_PTR( FLA_ONE );
  double*   minus_one_p = FLA_DOUBLE_PTR( FLA_MINUS_ONE );

  // if ( FLA_Obj_has_zero_dim( r12t ) ) return FLA_SUCCESS;
  if ( n_r12t == 0 ) return FLA_SUCCESS;

  // // w12t = r12t;
  // FLA_Copy_external( r12t, w12t );
  bl1_dcopyv( BLIS1_NO_CONJUGATE,
              n_r12t,
              r12t, inc_r12t, 
              w12t, inc_w12t ); 

  // // w12t = w12t + u1' * C2;
  // //      = w12t + C2^T * conj(u1);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, C2, u1, FLA_ONE, w12t );
  bl1_dgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_u1_C2,
             n_r12t,
             one_p,
             C2, rs_C2, cs_C2,
             u1, inc_u1,
             one_p,
             w12t, inc_w12t );

  // // w12t = w12t + v1' * D2;
  // //      = w12t + D2^T * conj(v1);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, D2, v1, FLA_ONE, w12t );
  bl1_dgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_v1_D2,
             n_r12t,
             one_p,
             D2, rs_D2, cs_D2,
             v1, inc_v1,
             one_p,
             w12t, inc_w12t );

  // // w12t = w12t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w12t );
  bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                 n_r12t,
                 tau,
                 w12t, inc_w12t );

  // // r12t = - w12t + r12t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w12t, r12t );
  bl1_daxpyv( BLIS1_NO_CONJUGATE,
              n_r12t,
              minus_one_p,
              w12t, inc_w12t,
              r12t, inc_r12t );

  // // C2 = - u1 * w12t + C2;
  // FLA_Ger_external( FLA_MINUS_ONE, u1, w12t, C2 );
  bl1_dger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_u1_C2,
            n_r12t,
            minus_one_p,
            u1, inc_u1,
            w12t, inc_w12t,
            C2, rs_C2, cs_C2 );

  // // D2 = v1 * w12t + D2;
  // FLA_Ger_external( FLA_ONE, v1, w12t, D2 );
  bl1_dger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_v1_D2,
            n_r12t,
            one_p,
            v1, inc_v1,
            w12t, inc_w12t,
            D2, rs_D2, cs_D2 );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_HUD_UT_l_opc_var1( int m_u1_C2,
                                       int m_v1_D2,
                                       int n_r12t,
                                       scomplex* tau,
                                       scomplex* w12t, int inc_w12t,
                                       scomplex* r12t, int inc_r12t,
                                       scomplex* u1, int inc_u1,
                                       scomplex* C2, int rs_C2, int cs_C2,
                                       scomplex* v1, int inc_v1,
                                       scomplex* D2, int rs_D2, int cs_D2 )
{
  scomplex* one_p       = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* minus_one_p = FLA_COMPLEX_PTR( FLA_MINUS_ONE );

  // if ( FLA_Obj_has_zero_dim( r12t ) ) return FLA_SUCCESS;
  if ( n_r12t == 0 ) return FLA_SUCCESS;

  // // w12t = r12t;
  // FLA_Copy_external( r12t, w12t );
  bl1_ccopyv( BLIS1_NO_CONJUGATE,
              n_r12t,
              r12t, inc_r12t, 
              w12t, inc_w12t ); 

  // // w12t = w12t + u1' * C2;
  // //      = w12t + C2^T * conj(u1);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, C2, u1, FLA_ONE, w12t );
  bl1_cgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_u1_C2,
             n_r12t,
             one_p,
             C2, rs_C2, cs_C2,
             u1, inc_u1,
             one_p,
             w12t, inc_w12t );

  // // w12t = w12t + v1' * D2;
  // //      = w12t + D2^T * conj(v1);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, D2, v1, FLA_ONE, w12t );
  bl1_cgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_v1_D2,
             n_r12t,
             one_p,
             D2, rs_D2, cs_D2,
             v1, inc_v1,
             one_p,
             w12t, inc_w12t );

  // // w12t = w12t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w12t );
  bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                 n_r12t,
                 tau,
                 w12t, inc_w12t );

  // // r12t = - w12t + r12t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w12t, r12t );
  bl1_caxpyv( BLIS1_NO_CONJUGATE,
              n_r12t,
              minus_one_p,
              w12t, inc_w12t,
              r12t, inc_r12t );

  // // C2 = - u1 * w12t + C2;
  // FLA_Ger_external( FLA_MINUS_ONE, u1, w12t, C2 );
  bl1_cger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_u1_C2,
            n_r12t,
            minus_one_p,
            u1, inc_u1,
            w12t, inc_w12t,
            C2, rs_C2, cs_C2 );

  // // D2 = v1 * w12t + D2;
  // FLA_Ger_external( FLA_ONE, v1, w12t, D2 );
  bl1_cger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_v1_D2,
            n_r12t,
            one_p,
            v1, inc_v1,
            w12t, inc_w12t,
            D2, rs_D2, cs_D2 );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_HUD_UT_l_opz_var1( int m_u1_C2,
                                       int m_v1_D2,
                                       int n_r12t,
                                       dcomplex* tau,
                                       dcomplex* w12t, int inc_w12t,
                                       dcomplex* r12t, int inc_r12t,
                                       dcomplex* u1, int inc_u1,
                                       dcomplex* C2, int rs_C2, int cs_C2,
                                       dcomplex* v1, int inc_v1,
                                       dcomplex* D2, int rs_D2, int cs_D2 )
{
  dcomplex* one_p       = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* minus_one_p = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );

  // if ( FLA_Obj_has_zero_dim( r12t ) ) return FLA_SUCCESS;
  if ( n_r12t == 0 ) return FLA_SUCCESS;

  // // w12t = r12t;
  // FLA_Copy_external( r12t, w12t );
  bl1_zcopyv( BLIS1_NO_CONJUGATE,
              n_r12t,
              r12t, inc_r12t, 
              w12t, inc_w12t ); 

  // // w12t = w12t + u1' * C2;
  // //      = w12t + C2^T * conj(u1);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, C2, u1, FLA_ONE, w12t );
  bl1_zgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_u1_C2,
             n_r12t,
             one_p,
             C2, rs_C2, cs_C2,
             u1, inc_u1,
             one_p,
             w12t, inc_w12t );

  // // w12t = w12t + v1' * D2;
  // //      = w12t + D2^T * conj(v1);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, D2, v1, FLA_ONE, w12t );
  bl1_zgemv( BLIS1_TRANSPOSE,
             BLIS1_CONJUGATE,
             m_v1_D2,
             n_r12t,
             one_p,
             D2, rs_D2, cs_D2,
             v1, inc_v1,
             one_p,
             w12t, inc_w12t );

  // // w12t = w12t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w12t );
  bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                 n_r12t,
                 tau,
                 w12t, inc_w12t );

  // // r12t = - w12t + r12t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w12t, r12t );
  bl1_zaxpyv( BLIS1_NO_CONJUGATE,
              n_r12t,
              minus_one_p,
              w12t, inc_w12t,
              r12t, inc_r12t );

  // // C2 = - u1 * w12t + C2;
  // FLA_Ger_external( FLA_MINUS_ONE, u1, w12t, C2 );
  bl1_zger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_u1_C2,
            n_r12t,
            minus_one_p,
            u1, inc_u1,
            w12t, inc_w12t,
            C2, rs_C2, cs_C2 );

  // // D2 = v1 * w12t + D2;
  // FLA_Ger_external( FLA_ONE, v1, w12t, D2 );
  bl1_zger( BLIS1_NO_CONJUGATE,
            BLIS1_NO_CONJUGATE,
            m_v1_D2,
            n_r12t,
            one_p,
            v1, inc_v1,
            w12t, inc_w12t,
            D2, rs_D2, cs_D2 );

  return FLA_SUCCESS;
}
