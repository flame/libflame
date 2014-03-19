/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_l_unb_var2( FLA_Obj A, FLA_Obj T )
{
  return FLA_Tridiag_UT_l_step_unb_var2( A, T );
}

FLA_Error FLA_Tridiag_UT_l_step_unb_var2( FLA_Obj A, FLA_Obj T )
{
  FLA_Obj  ATL,   ATR,      A00,  a01,     A02, 
           ABL,   ABR,      a10t, alpha11, a12t,
                            A20,  a21,     A22;
  FLA_Obj  TTL,   TTR,      T00,  t01,   T02, 
           TBL,   TBR,      t10t, tau11, t12t,
                            T20,  t21,   T22;
  FLA_Obj  uT,              u01,
           uB,              upsilon11,
                            u21;
  FLA_Obj  zT,              z01,
           zB,              zeta11,
                            z21;
  FLA_Obj  wT,              w01,
           wB,              omega11,
                            w21;
  FLA_Obj  u, z, w;

  FLA_Obj  inv_tau11;
  FLA_Obj  minus_inv_tau11;
  FLA_Obj  first_elem;
  FLA_Obj  beta;
  FLA_Obj  minus_upsilon11;
  FLA_Obj  minus_conj_upsilon11;
  FLA_Obj  minus_zeta11;
  FLA_Obj  minus_conj_zeta11;

  FLA_Obj  a21_t,
           a21_b;

  FLA_Datatype datatype_A;
  dim_t        m_A;
  dim_t        b_alg;


  b_alg      = FLA_Obj_length( T );

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );

  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &inv_tau11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_inv_tau11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &first_elem );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &beta );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_upsilon11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_conj_upsilon11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_zeta11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_conj_zeta11 );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );
  FLA_Part_2x1( u,    &uT, 
                      &uB,            0, FLA_TOP );
  FLA_Part_2x1( z,    &zT, 
                      &zB,            0, FLA_TOP );
  FLA_Part_2x1( w,    &wT, 
                      &wB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < b_alg )
  {
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************ */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );
    FLA_Repart_2x1_to_3x1( uT,                &u01, 
                        /* ** */            /* ********* */
                                              &upsilon11, 
                           uB,                &u21,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( zT,                &z01, 
                        /* ** */            /* ****** */
                                              &zeta11, 
                           zB,                &z21,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( wT,                &w01, 
                        /* ** */            /* ******* */
                                              &omega11, 
                           wB,                &w21,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    
    if ( FLA_Obj_length( ATL ) > 0 )
    {
      FLA_Copy( upsilon11, minus_upsilon11 );
      FLA_Scal( FLA_MINUS_ONE, minus_upsilon11 );
      FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon11, minus_conj_upsilon11 );

      FLA_Copy( zeta11, minus_zeta11 );
      FLA_Scal( FLA_MINUS_ONE, minus_zeta11 );
      FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, minus_zeta11, minus_conj_zeta11 );

      // alpha11 = alpha11 - upsilon11 * conj(zeta11) - zeta11 * conj(upsilon11);
      FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon11, zeta11,    alpha11 );
      FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_zeta11,    upsilon11, alpha11 );

      // a21 = a21 - conj(zeta11) * u21 - conj(upsilon11) * z21;
      FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_zeta11,    u21, a21 );
      FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_upsilon11, z21, a21 );
    }

    if ( FLA_Obj_length( A22 ) > 0 )
    {
      FLA_Part_2x1( a21,    &a21_t,
                            &a21_b,            1, FLA_TOP );

      // [ x21, tau11, a21 ] = House( a21 );
      FLA_Househ2_UT( FLA_LEFT,
                      a21_t,
                      a21_b, tau11 );

      // inv_tau11       =  1 / tau11;
      // minus_inv_tau11 = -1 / tau11;
      FLA_Set( FLA_ONE, inv_tau11 );
      FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      FLA_Copy( inv_tau11, minus_inv_tau11 );
      FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );

      // Save first element of a21_t and set it to one so we can use a21 as
      // u21 in subsequent computations. We will restore a21_t later on.
      FLA_Copy( a21_t, first_elem );
      FLA_Set( FLA_ONE, a21_t );
    }

    if ( FLA_Obj_length( ATL ) > 0 )
    {
      // A22 = A22 - u21 * z21' - z21 * u21';
      FLA_Her2( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, u21, z21, A22 );
    }

    if ( FLA_Obj_length( A22 ) > 0 )
    {
      // w21 = A22 * x21;
      FLA_Hemv( FLA_LOWER_TRIANGULAR, FLA_ONE, A22, a21, FLA_ZERO, w21 );

      // u21 = x21; 
      // z21 = w21;
      FLA_Copy( a21, u21 );
      FLA_Copy( w21, z21 );

      // beta = u21' * z21 / 2;
      FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      FLA_Inv_scal( FLA_TWO, beta );

      // z21 = ( z21 - beta / tau * u21 ) / tau;
      FLA_Scal( minus_inv_tau11, beta );
      FLA_Axpy( beta, a21, z21 );
      FLA_Scal( inv_tau11, z21 );

      // t01 = U20' * u21;
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );

      // Restore first element of a21.
      FLA_Copy( first_elem, a21_t );
    }

    // Update A22 if this is the last iteration; this is needed when we're
    // being called from the blocked routine so A22 is left in a valid state.
    if ( FLA_Obj_length( ATL ) + 1 == b_alg &&
         FLA_Obj_length( A22 ) > 0 )
    {
      // A22 = A22 - u21 * y21' - z21 * u21';
      FLA_Her2( FLA_LOWER_TRIANGULAR, FLA_MINUS_ONE, u21, z21, A22 );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ********************** */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &uT,                u01, 
                                                  upsilon11, 
                            /* ** */           /* ********* */
                              &uB,                u21,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &zT,                z01, 
                                                  zeta11, 
                            /* ** */           /* ****** */
                              &zB,                z21,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &wT,                w01, 
                                                  omega11, 
                            /* ** */           /* ******* */
                              &wB,                w21,     FLA_TOP );
  }

  FLA_Obj_free( &inv_tau11 );
  FLA_Obj_free( &minus_inv_tau11 );
  FLA_Obj_free( &first_elem );
  FLA_Obj_free( &beta );
  FLA_Obj_free( &minus_upsilon11 );
  FLA_Obj_free( &minus_conj_upsilon11 );
  FLA_Obj_free( &minus_zeta11 );
  FLA_Obj_free( &minus_conj_zeta11 );
  FLA_Obj_free( &u );
  FLA_Obj_free( &z );
  FLA_Obj_free( &w );

  return FLA_SUCCESS;
}


