/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_H2_UT_r_opt_var1( FLA_Obj tau, FLA_Obj u2, FLA_Obj a1, FLA_Obj A2 )
{
  FLA_Datatype datatype;
  int          n_u2_A2;
  int          m_a1;
  int          inc_u2;
  int          inc_a1;
  int          rs_A2, cs_A2;

  // The house-holder transformation in libFLAME never creates a zero tau value.
  // However, when libFLAME is mixed with LAPACK, zero tau means to apply an 
  // identity matrix that does nothing here.
  if ( FLA_Obj_has_zero_dim( a1 ) || FLA_Obj_equals( tau, FLA_ZERO) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A2 );

  n_u2_A2  = FLA_Obj_width( A2 );
  m_a1     = FLA_Obj_length( a1 );
  inc_u2   = FLA_Obj_vector_inc( u2 );
  inc_a1   = FLA_Obj_vector_inc( a1 );
  rs_A2    = FLA_Obj_row_stride( A2 );
  cs_A2    = FLA_Obj_col_stride( A2 );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* tau_p = ( float* ) FLA_FLOAT_PTR( tau );
      float* u2_p  = ( float* ) FLA_FLOAT_PTR( u2 );
      float* a1_p  = ( float* ) FLA_FLOAT_PTR( a1 );
      float* A2_p  = ( float* ) FLA_FLOAT_PTR( A2 );

      if ( *tau_p != 0.0F )
        FLA_Apply_H2_UT_r_ops_var1( m_a1, n_u2_A2,
                                    tau_p,
                                    u2_p, inc_u2,
                                    a1_p, inc_a1,
                                    A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_DOUBLE:
    {
      double* tau_p = ( double* ) FLA_DOUBLE_PTR( tau );
      double* u2_p  = ( double* ) FLA_DOUBLE_PTR( u2 );
      double* a1_p  = ( double* ) FLA_DOUBLE_PTR( a1 );
      double* A2_p  = ( double* ) FLA_DOUBLE_PTR( A2 );

      if ( *tau_p != 0.0 ) 
        FLA_Apply_H2_UT_r_opd_var1( m_a1, n_u2_A2,
                                    tau_p,
                                    u2_p, inc_u2,
                                    a1_p, inc_a1,
                                    A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* tau_p = ( scomplex* ) FLA_COMPLEX_PTR( tau );
      scomplex* u2_p  = ( scomplex* ) FLA_COMPLEX_PTR( u2 );
      scomplex* a1_p  = ( scomplex* ) FLA_COMPLEX_PTR( a1 );
      scomplex* A2_p  = ( scomplex* ) FLA_COMPLEX_PTR( A2 );

      if ( tau_p->real != 0.0F && tau_p->imag != 0.0F )
        FLA_Apply_H2_UT_r_opc_var1( m_a1, n_u2_A2,
                                    tau_p,
                                    u2_p, inc_u2,
                                    a1_p, inc_a1,
                                    A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* tau_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( tau );
      dcomplex* u2_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( u2 );
      dcomplex* a1_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( a1 );
      dcomplex* A2_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A2 );

      if ( tau_p->real != 0.0 && tau_p->imag != 0.0 )
        FLA_Apply_H2_UT_r_opz_var1( m_a1, n_u2_A2,
                                    tau_p,
                                    u2_p, inc_u2,
                                    a1_p, inc_a1,
                                    A2_p, rs_A2, cs_A2 );
      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_r_ops_var1( int m_a1,
                                      int n_u2_A2,
                                      float* tau,
                                      float* u2, int inc_u2,
                                      float* a1, int inc_a1,
                                      float* A2, int rs_A2, int cs_A2 )
{
  float*    one_p       = FLA_FLOAT_PTR( FLA_ONE );
  float*    minus_one_p = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       inc_w1;

  // FLA_Obj w1;
  float*    w1;

  if ( m_a1 == 0 || *tau == 0.0F ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1, &w1 );
  w1 = ( float* ) FLA_malloc( m_a1 * sizeof( *a1 ) );
  inc_w1 = 1;

  // // w1 = a1;
  // FLA_Copy_external( a1, w1 );
  bl1_scopyv( BLIS1_NO_CONJUGATE,
              m_a1,
              a1, inc_a1, 
              w1, inc_w1 ); 

  // // w1 = w1 + A2 * u2;
  // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A2, u2, FLA_ONE, w1 );
  bl1_sgemv( BLIS1_NO_TRANSPOSE,
             BLIS1_NO_CONJUGATE,
             m_a1,
             n_u2_A2,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1, inc_w1 );

  // // w1 = w1 / tau;
  // FLA_Inv_scal_external( tau, w1 );
  bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                 m_a1,
                 tau,
                 w1, inc_w1 );

  // // a1 = a1 - w1;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1, a1 );
  bl1_saxpyv( BLIS1_NO_CONJUGATE,
              m_a1,
              minus_one_p,
              w1, inc_w1,
              a1, inc_a1 );

  // // A2 = A2 - w1 * u2';
  // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, w1, u2, A2 );
  bl1_sger( BLIS1_NO_CONJUGATE,
            BLIS1_CONJUGATE,
            m_a1,
            n_u2_A2,
            minus_one_p,
            w1, inc_w1,
            u2, inc_u2,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1 );
  FLA_free( w1 );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_r_opd_var1( int m_a1,
                                      int n_u2_A2,
                                      double* tau,
                                      double* u2, int inc_u2,
                                      double* a1, int inc_a1,
                                      double* A2, int rs_A2, int cs_A2 )
{
  double*   one_p       = FLA_DOUBLE_PTR( FLA_ONE );
  double*   minus_one_p = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       inc_w1;

  // FLA_Obj w1;
  double*   w1;

  if ( m_a1 == 0 || *tau == 0.0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1, &w1 );
  w1 = ( double* ) FLA_malloc( m_a1 * sizeof( *a1 ) );
  inc_w1 = 1;

  // // w1 = a1;
  // FLA_Copy_external( a1, w1 );
  bl1_dcopyv( BLIS1_NO_CONJUGATE,
              m_a1,
              a1, inc_a1, 
              w1, inc_w1 ); 

  // // w1 = w1 + A2 * u2;
  // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A2, u2, FLA_ONE, w1 );
  bl1_dgemv( BLIS1_NO_TRANSPOSE,
             BLIS1_NO_CONJUGATE,
             m_a1,
             n_u2_A2,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1, inc_w1 );

  // // w1 = w1 / tau;
  // FLA_Inv_scal_external( tau, w1 );
  bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                 m_a1,
                 tau,
                 w1, inc_w1 );

  // // a1 = a1 - w1;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1, a1 );
  bl1_daxpyv( BLIS1_NO_CONJUGATE,
              m_a1,
              minus_one_p,
              w1, inc_w1,
              a1, inc_a1 );

  // // A2 = A2 - w1 * u2';
  // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, w1, u2, A2 );
  bl1_dger( BLIS1_NO_CONJUGATE,
            BLIS1_CONJUGATE,
            m_a1,
            n_u2_A2,
            minus_one_p,
            w1, inc_w1,
            u2, inc_u2,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1 );
  FLA_free( w1 );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_r_opc_var1( int m_a1,
                                      int n_u2_A2,
                                      scomplex* tau,
                                      scomplex* u2, int inc_u2,
                                      scomplex* a1, int inc_a1,
                                      scomplex* A2, int rs_A2, int cs_A2 )
{
  scomplex* one_p       = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* minus_one_p = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       inc_w1;

  // FLA_Obj w1;
  scomplex* w1;

  if ( m_a1 == 0 || ( tau->real == 0.0F && tau->imag == 0.0F ) ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1, &w1 );
  w1 = ( scomplex* ) FLA_malloc( m_a1 * sizeof( *a1 ) );
  inc_w1 = 1;

  // // w1 = a1;
  // FLA_Copy_external( a1, w1 );
  bl1_ccopyv( BLIS1_NO_CONJUGATE,
              m_a1,
              a1, inc_a1, 
              w1, inc_w1 ); 

  // // w1 = w1 + A2 * u2;
  // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A2, u2, FLA_ONE, w1 );
  bl1_cgemv( BLIS1_NO_TRANSPOSE,
             BLIS1_NO_CONJUGATE,
             m_a1,
             n_u2_A2,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1, inc_w1 );

  // // w1 = w1 / tau;
  // FLA_Inv_scal_external( tau, w1 );
  bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                 m_a1,
                 tau,
                 w1, inc_w1 );

  // // a1 = a1 - w1;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1, a1 );
  bl1_caxpyv( BLIS1_NO_CONJUGATE,
              m_a1,
              minus_one_p,
              w1, inc_w1,
              a1, inc_a1 );

  // // A2 = A2 - w1 * u2';
  // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, w1, u2, A2 );
  bl1_cger( BLIS1_NO_CONJUGATE,
            BLIS1_CONJUGATE,
            m_a1,
            n_u2_A2,
            minus_one_p,
            w1, inc_w1,
            u2, inc_u2,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1 );
  FLA_free( w1 );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_r_opz_var1( int m_a1,
                                      int n_u2_A2,
                                      dcomplex* tau,
                                      dcomplex* u2, int inc_u2,
                                      dcomplex* a1, int inc_a1,
                                      dcomplex* A2, int rs_A2, int cs_A2 )
{
  dcomplex* one_p       = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* minus_one_p = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       inc_w1;

  // FLA_Obj w1;
  dcomplex* w1;

  if ( m_a1 == 0 || ( tau->real == 0.0 && tau->imag == 0.0 ) ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1, &w1 );
  w1 = ( dcomplex* ) FLA_malloc( m_a1 * sizeof( *a1 ) );
  inc_w1 = 1;

  // // w1 = a1;
  // FLA_Copy_external( a1, w1 );
  bl1_zcopyv( BLIS1_NO_CONJUGATE,
              m_a1,
              a1, inc_a1, 
              w1, inc_w1 ); 

  // // w1 = w1 + A2 * u2;
  // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A2, u2, FLA_ONE, w1 );
  bl1_zgemv( BLIS1_NO_TRANSPOSE,
             BLIS1_NO_CONJUGATE,
             m_a1,
             n_u2_A2,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1, inc_w1 );

  // // w1 = w1 / tau;
  // FLA_Inv_scal_external( tau, w1 );
  bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                 m_a1,
                 tau,
                 w1, inc_w1 );

  // // a1 = a1 - w1;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1, a1 );
  bl1_zaxpyv( BLIS1_NO_CONJUGATE,
              m_a1,
              minus_one_p,
              w1, inc_w1,
              a1, inc_a1 );

  // // A2 = A2 - w1 * u2';
  // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, w1, u2, A2 );
  bl1_zger( BLIS1_NO_CONJUGATE,
            BLIS1_CONJUGATE,
            m_a1,
            n_u2_A2,
            minus_one_p,
            w1, inc_w1,
            u2, inc_u2,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1 );
  FLA_free( w1 );

  return FLA_SUCCESS;
}
