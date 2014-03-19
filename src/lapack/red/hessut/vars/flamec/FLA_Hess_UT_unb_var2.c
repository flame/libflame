/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Hess_UT_unb_var2( FLA_Obj A, FLA_Obj T )
{
  return FLA_Hess_UT_step_unb_var2( A, T );
}

FLA_Error FLA_Hess_UT_step_unb_var2( FLA_Obj A, FLA_Obj T )
{
  FLA_Obj  ATL,   ATR,      A00,  a01,     A02, 
           ABL,   ABR,      a10t, alpha11, a12t,
                            A20,  a21,     A22;
  FLA_Obj  TTL,   TTR,      T00,  t01,   T02, 
           TBL,   TBR,      t10t, tau11, t12t,
                            T20,  t21,   T22;
  FLA_Obj  yT,              y0,
           yB,              psi1,
                            y2;
  FLA_Obj  zT,              z0,
           zB,              zeta1,
                            z2;
  FLA_Obj  y, z;
           
  FLA_Obj  inv_tau11;
  FLA_Obj  minus_inv_tau11;
  FLA_Obj  first_elem;
  FLA_Obj  beta;
  FLA_Obj  conj_beta;
  FLA_Obj  dot_product;

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
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &conj_beta );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &dot_product );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &y );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );
  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_TOP );
  FLA_Part_2x1( z,    &zT, 
                      &zB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < b_alg )
  {
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************** */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );
    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                        /* ** */            /* **** */
                                              &psi1, 
                           yB,                &y2,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( zT,                &z0, 
                        /* ** */            /* ***** */
                                              &zeta1, 
                           zB,                &z2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    if ( FLA_Obj_length( A22 ) > 0 )
    {
      FLA_Part_2x1( a21,    &a21_t,
                            &a21_b,        1, FLA_TOP );

      // [ u21, tau11, a21 ] = House( a21 );
      FLA_Househ2_UT( FLA_LEFT,
                      a21_t,
                      a21_b, tau11 );

      // inv_tau11            =  1 / tau11;
      // minus_inv_tau11      = -1 / tau11;
      FLA_Set( FLA_ONE, inv_tau11 );
      FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      FLA_Copy( inv_tau11, minus_inv_tau11 );
      FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );

      // Save first element of a21_t and set it to one so we can use a21 as
      // u21 in subsequent computations. We will restore a21_t later on.
      FLA_Copy( a21_t, first_elem );
      FLA_Set( FLA_ONE, a21_t );

      // y21 = A22' * u21;
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, y2 );

      // z21 = A22 * u21;
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, z2 );

      // beta      = u21' * z21 / 2;
      // conj_beta = conj(beta);
      FLA_Dotc( FLA_CONJUGATE, a21, z2, beta );
      FLA_Inv_scal( FLA_TWO, beta );
      FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );

      // y21' = ( y21' - beta / tau * u21' ) / tau;
      // y21  = ( y21 - conj(beta) / tau * u21 ) / tau;
      FLA_Scal( minus_inv_tau11, conj_beta );
      FLA_Axpy( conj_beta, a21, y2 );
      FLA_Scal( inv_tau11, y2 );

      // z21 = ( z21 - beta / tau * u21 ) / tau;
      FLA_Scal( minus_inv_tau11, beta );
      FLA_Axpy( beta, a21, z2 );
      FLA_Scal( inv_tau11, z2 );

      // a12t = a12t * ( I - u21 * u21' / tau );
      //      = a12t - ( a12t * u21 ) * u21' / tau;
      FLA_Dot( a12t, a21, dot_product );
      FLA_Scal( minus_inv_tau11, dot_product );
      FLA_Axpyt( FLA_CONJ_TRANSPOSE, dot_product, a21, a12t );

      // A02 = A02 * ( I - u21 * u21' / tau );
      //     = A02 - ( A02 * u21 ) * u21' / tau;
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, y0 );
      FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, y0, a21, A02 );

      // A22 = A22 - u21 * y21' - z21 * u21';
      FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, a21, y2, A22 );
      FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, a21, A22 );

      // t01 = U20' * u21;
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );

      // Restore first element of a21.
      FLA_Copy( first_elem, a21_t );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ************************ */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                                                  psi1, 
                            /* ** */           /* **** */
                              &yB,                y2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &zT,                z0, 
                                                  zeta1, 
                            /* ** */           /* ***** */
                              &zB,                z2,     FLA_TOP );
  }

  FLA_Obj_free( &inv_tau11 );
  FLA_Obj_free( &minus_inv_tau11 );
  FLA_Obj_free( &first_elem );
  FLA_Obj_free( &beta );
  FLA_Obj_free( &conj_beta );
  FLA_Obj_free( &dot_product );
  FLA_Obj_free( &y );
  FLA_Obj_free( &z );

  return FLA_SUCCESS;
}

