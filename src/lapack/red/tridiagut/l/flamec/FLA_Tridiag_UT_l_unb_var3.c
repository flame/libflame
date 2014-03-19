/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_l_unb_var3( FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;
  FLA_Obj   Z;

  //FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Y );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Z );

  r_val = FLA_Tridiag_UT_l_step_unb_var3( A, Z, T );

  //FLA_Obj_free( &Y );
  FLA_Obj_free( &Z );

  return r_val;
}

FLA_Error FLA_Tridiag_UT_l_step_unb_var3( FLA_Obj A, FLA_Obj Z, FLA_Obj T )
{
  FLA_Obj  ATL,   ATR,      A00,  a01,     A02, 
           ABL,   ABR,      a10t, alpha11, a12t,
                            A20,  a21,     A22;
  FLA_Obj  ZTL,   ZTR,      Z00,  z011,    Z02, 
           ZBL,   ZBR,      z10t, zeta11, z12t,
                            Z20,  z21,    Z22;
  FLA_Obj  TTL,   TTR,      T00,  t01,   T02, 
           TBL,   TBR,      t10t, tau11, t12t,
                            T20,  t21,   T22;
  FLA_Obj  dT,              d01,
           dB,              delta11,
                            d21;
  FLA_Obj  fT,              f01,
           fB,              phi11,
                            f21;
  FLA_Obj  d, f;

  FLA_Obj  inv_tau11;
  FLA_Obj  minus_inv_tau11;
  FLA_Obj  beta;
  FLA_Obj  first_elem;
  FLA_Obj  last_elem;

  FLA_Obj  a10t_l, a10t_r;
  FLA_Obj  a21_t,
           a21_b;
  FLA_Obj  a2;

  FLA_Datatype datatype_A;
  dim_t        m_A;
  dim_t        b_alg;


  b_alg      = FLA_Obj_length( T );

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );

  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &inv_tau11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_inv_tau11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &beta );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &first_elem );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &last_elem );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &d );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );

  FLA_Set( FLA_ZERO, Z );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x2( Z,    &ZTL, &ZTR,
                      &ZBL, &ZBR,     0, 0, FLA_TL );
  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );
  FLA_Part_2x1( d,    &dT, 
                      &dB,            0, FLA_TOP );
  FLA_Part_2x1( f,    &fT, 
                      &fB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < b_alg )
  {
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( ZTL, /**/ ZTR,       &Z00,  /**/ &z011,    &Z02,
                        /* ************* */   /* ************************* */
                                                &z10t, /**/ &zeta11, &z12t,
                           ZBL, /**/ ZBR,       &Z20,  /**/ &z21,    &Z22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************ */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );
    FLA_Repart_2x1_to_3x1( dT,                &d01, 
                        /* ** */            /* ******* */
                                              &delta11, 
                           dB,                &d21,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( fT,                &f01, 
                        /* ** */            /* ***** */
                                              &phi11, 
                           fB,                &f21,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // Save first element of a10_r and set it to one so we can use a10t as
    // u10t in subsequent computations. We will restore a10_r later on.
    if ( FLA_Obj_length( ATL ) > 0 )
    {
      FLA_Part_1x2( a10t,   &a10t_l, &a10t_r,     1, FLA_RIGHT );
      FLA_Copy( a10t_r, last_elem );
      FLA_Set( FLA_ONE, a10t_r );
    }
  
    FLA_Merge_2x1( alpha11,
                   a21,      &a2 );

    // alpha11 = alpha11 - u10t * z10t' - z10t * u10t';
    // a21     = a21     - U20  * z10t' - Z20  * u10t';
    FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, z10t, FLA_ONE, a2 );
    FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a10t, FLA_ONE, a2 );

    // Restore last element of a10t.
    if ( FLA_Obj_length( ATL ) > 0 )
    {
      FLA_Copy( last_elem, a10t_r );
    }

    if ( FLA_Obj_length( A22 ) > 0 )
    {
      FLA_Part_2x1( a21,    &a21_t,
                            &a21_b,            1, FLA_TOP );

      // [ u21, tau11, a21 ] = House( a21 );
      FLA_Househ2_UT( FLA_LEFT,
                      a21_t,
                      a21_b, tau11 );
  
      // inv_tau11       =  1 / tau11;
      // minus_inv_tau11 = -1 / tau11;
      FLA_Set( FLA_ONE, inv_tau11 );
      FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      FLA_Copy( inv_tau11, minus_inv_tau11 );
      FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );

      // Save first element of a21_t and set it to one.
      FLA_Copy( a21_t, first_elem );
      FLA_Set( FLA_ONE, a21_t );
  
      // z21 = A22 * u21;
      FLA_Hemv( FLA_LOWER_TRIANGULAR, FLA_ONE, A22, a21, FLA_ZERO, z21 );
 
      // z21 = z21 - U20 * ( Z20' * u21 ) - Z20 * ( U20' * u21 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, d01 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, a21, FLA_ZERO, f01 );
  
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f01, FLA_ONE, z21 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, d01, FLA_ONE, z21 );

      // t01 = U20' * u21;
      FLA_Copy( d01, t01 );
  
      // beta = u21' * z21 / 2;
      FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      FLA_Inv_scal( FLA_TWO, beta );
  
      // z21 = ( z21 - beta / tau * u21 ) / tau;
      FLA_Scal( minus_inv_tau11, beta );
      FLA_Axpy( beta, a21, z21 );
      FLA_Scal( inv_tau11, z21 );
  
      // Restore first element of a21.
      FLA_Copy( first_elem, a21_t );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &ZTL, /**/ &ZTR,       Z00,  z011,    /**/ Z02,
                                                     z10t, zeta11, /**/ z12t,
                            /* ************** */  /* *********************** */
                              &ZBL, /**/ &ZBR,       Z20,  z21,    /**/ Z22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ********************** */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &dT,                d01, 
                                                  delta11, 
                            /* ** */           /* ******* */
                              &dB,                d21,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &fT,                f01, 
                                                  phi11, 
                            /* ** */           /* ***** */
                              &fB,                f21,     FLA_TOP );
  }

  FLA_Obj_free( &inv_tau11 );
  FLA_Obj_free( &minus_inv_tau11 );
  FLA_Obj_free( &beta );
  FLA_Obj_free( &first_elem );
  FLA_Obj_free( &last_elem );
  FLA_Obj_free( &d );
  FLA_Obj_free( &f );

  return FLA_SUCCESS;
}

