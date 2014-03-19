/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_unb_var5( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  FLA_Error    r_val;
  FLA_Obj      Y, Z;
  FLA_Datatype datatype_A;
  dim_t        m_A, n_A;

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );
  n_A        = FLA_Obj_width( A );
  
  FLA_Obj_create( datatype_A, n_A, n_A, 0, 0, &Y );
  FLA_Obj_create( datatype_A, m_A, n_A, 0, 0, &Z );

  r_val = FLA_Bidiag_UT_u_step_unb_var5( A, Y, Z, TU, TV );

  FLA_Obj_free( &Y );
  FLA_Obj_free( &Z );

  return r_val;
}

FLA_Error FLA_Bidiag_UT_u_step_unb_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj S )
{
  FLA_Obj  ATL,   ATR,      A00,  a01,     A02, 
           ABL,   ABR,      a10t, alpha11, a12t,
                            A20,  a21,     A22;
  FLA_Obj  YTL,   YTR,      Y00,  y01,   Y02, 
           YBL,   YBR,      y10t, psi11, y12t,
                            Y20,  y21,   Y22;
  FLA_Obj  ZTL,   ZTR,      Z00,  z01,    Z02, 
           ZBL,   ZBR,      z10t, zeta11, z12t,
                            Z20,  z21,    Z22;
  FLA_Obj  TTL,   TTR,      T00,  t01,   T02, 
           TBL,   TBR,      t10t, tau11, t12t,
                            T20,  t21,   T22;
  FLA_Obj  STL,   STR,      S00,  s01,     S02, 
           SBL,   SBR,      s10t, sigma11, s12t,
                            S20,  s21,     S22;
  FLA_Obj  uT,              u01,
           uB,              upsilon11,
                            u21;
  FLA_Obj  vT,              v01,
           vB,              nu11,
                            v21;
  FLA_Obj  dT,              d0,
           dB,              delta1,
                            d2;
  FLA_Obj  eT,              e0,
           eB,              epsilon1,
                            e2;
  FLA_Obj  fT,              f0,
           fB,              phi1,
                            f2;
  FLA_Obj  gT,              g0,
           gB,              ghi1,
                            g2;
  FLA_Obj  u, v;
  FLA_Obj  d, e, f, g;

  FLA_Obj  last_elem;
  FLA_Obj  beta;
  FLA_Obj  minus_upsilon11;
  FLA_Obj  minus_zeta11;

  FLA_Obj  a01_t,
           a01_b;
  FLA_Obj  a12t_l, a12t_r;
  FLA_Obj  v21_t,
           v21_b;
  FLA_Obj  a2;

  FLA_Datatype datatype_A;
  dim_t        m_A, n_A;
  dim_t        b_alg;


  b_alg      = FLA_Obj_length( T );

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );
  n_A        = FLA_Obj_width( A );

  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &last_elem );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &beta );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_upsilon11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_zeta11 );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &d );
  FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &e );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &g );

  FLA_Set( FLA_ZERO, Y );
  FLA_Set( FLA_ZERO, Z );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x2( Y,    &YTL, &YTR,
                      &YBL, &YBR,     0, 0, FLA_TL );
  FLA_Part_2x2( Z,    &ZTL, &ZTR,
                      &ZBL, &ZBR,     0, 0, FLA_TL );
  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );
  FLA_Part_2x2( S,    &STL, &STR,
                      &SBL, &SBR,     0, 0, FLA_TL );
  FLA_Part_2x1( u,    &uT, 
                      &uB,            0, FLA_TOP );
  FLA_Part_2x1( v,    &vT, 
                      &vB,            0, FLA_TOP );
  FLA_Part_2x1( d,    &dT,
                      &dB,            0, FLA_TOP );
  FLA_Part_2x1( e,    &eT, 
                      &eB,            0, FLA_TOP );
  FLA_Part_2x1( f,    &fT,
                      &fB,            0, FLA_TOP );
  FLA_Part_2x1( g,    &gT,
                      &gB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < b_alg )
  {
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( YTL, /**/ YTR,       &Y00,  /**/ &y01,   &Y02,
                        /* ************* */   /* ************************ */
                                                &y10t, /**/ &psi11, &y12t,
                           YBL, /**/ YBR,       &Y20,  /**/ &y21,   &Y22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( ZTL, /**/ ZTR,       &Z00,  /**/ &z01,    &Z02,
                        /* ************* */   /* ************************* */
                                                &z10t, /**/ &zeta11, &z12t,
                           ZBL, /**/ ZBR,       &Z20,  /**/ &z21,    &Z22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************** */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( STL, /**/ STR,       &S00,  /**/ &s01,     &S02,
                        /* ************* */   /* ************************** */
                                                &s10t, /**/ &sigma11, &s12t,
                           SBL, /**/ SBR,       &S20,  /**/ &s21,     &S22,
                           1, 1, FLA_BR );
    FLA_Repart_2x1_to_3x1( uT,                &u01, 
                        /* ** */            /* ***** */
                                              &upsilon11, 
                           uB,                &u21,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( vT,                &v01, 
                        /* ** */            /* ***** */
                                              &nu11, 
                           vB,                &v21,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( dT,                &d0,
                        /* ** */            /* ****** */
                                              &delta1,
                           dB,                &d2,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( eT,                &e0,
                        /* ** */            /* ******** */
                                              &epsilon1,
                           eB,                &e2,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( fT,                &f0,
                        /* ** */            /* **** */
                                              &phi1,
                           fB,                &f2,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( gT,                &g0,
                        /* ** */            /* **** */
                                              &ghi1,
                           gB,                &g2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // Save last element of a01 and set it to one so we can use a01 as
    // v10t^T in subsequent computations. We will restore a01_b later on.
    // Also note: V20^T is stored in A02.
    if ( FLA_Obj_length( ATL ) > 0 )
    {
      FLA_Part_2x1( a01,    &a01_t,
                            &a01_b,            1, FLA_BOTTOM );
      FLA_Copy( a01_b, last_elem );
      FLA_Set( FLA_ONE, a01_b );
    }
    
    FLA_Merge_2x1( alpha11,
                   a21,      &a2 );

    // alpha11 = alpha11 - u10t * y10t' - z10t * v10t';
    // a21     = a21     - U20  * y10t' - Z20  * v10t';
    FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a01,  FLA_ONE, a2 );

    // a12t = a12t - u10t * Y20' - z10t * V20';
    FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    FLA_Gemv( FLA_CONJ_TRANSPOSE,    FLA_MINUS_ONE, A02, z10t, FLA_ONE, a12t );
    
    // Restore last element of a01.
    if ( FLA_Obj_length( ATL ) > 0 )
    {
      FLA_Copy( last_elem, a01_b );
    }

    // [ alpha11, u21, tau11 ] = House2( alpha11, a21 );
    FLA_Househ2_UT( FLA_LEFT,
                    alpha11,
                    a21, tau11 );
    FLA_Copy( a21, u21 );

    if ( FLA_Obj_width( A22 ) > 0 )
    {
      // y21' = a12t + u21' * A22;
      // y21  = conj(a12t) + A22' * u21;
      FLA_Copyt( FLA_CONJ_TRANSPOSE, a12t, y21 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, u21, FLA_ONE, y21 );

      // y21 = y21 - Y20 * ( U20' * u21 ) - V20 * ( Z20' * u21 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ZERO, d0 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, u21, FLA_ZERO, e0 );

      // t01 = a10t' + U20' * u21;
      FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      FLA_Axpy( FLA_ONE, d0, t01 );

      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      FLA_Gemv( FLA_TRANSPOSE,    FLA_MINUS_ONE, A02, e0, FLA_ONE, y21 );

      // y21 = y21 / tau11;
      FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, y21 );

      // a12t = a12t - conj(y21)^T;
      FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_MINUS_ONE, y21, a12t );

      FLA_Part_1x2( a12t,    &a12t_l, &a12t_r,      1, FLA_LEFT );
      FLA_Part_2x1( v21,     &v21_t, 
                             &v21_b,            1, FLA_TOP );

      // [ a12t_l, v21_b, sigma11 ] = House2( a12t_l, a12t_r );
      FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );

      // v21_t = 1;
      // v21_b = a12t_r^T;
      FLA_Set( FLA_ONE, v21_t );
      FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );

      // beta = - y21' * v21;
      FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      FLA_Scal( FLA_MINUS_ONE, beta );

      // z21 = A22 * v21 + beta * u21;
      FLA_Copy( u21, z21 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, v21, beta, z21 );

      // z21 = z21 - U20 * ( Y20' * v21 ) - Z20 * ( V20' * v21 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE,    FLA_ONE, Y20, v21, FLA_ZERO, f0 );
      FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, g0 );

      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, z21 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, g0, FLA_ONE, z21 );

      // z21 = z21 / sigma11;
      FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );

      // s01 = conj(V02) * v21;
      FLA_Copy( g0, s01 );
    }
    else // if ( FLA_Obj_width( A22 ) == 0 )
    {
      // t01 = a10t' + U20' * u21;
      FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ONE, t01 );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &YTL, /**/ &YTR,       Y00,  y01,   /**/ Y02,
                                                     y10t, psi11, /**/ y12t,
                            /* ************** */  /* ********************** */
                              &YBL, /**/ &YBR,       Y20,  y21,   /**/ Y22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &ZTL, /**/ &ZTR,       Z00,  z01,    /**/ Z02,
                                                     z10t, zeta11, /**/ z12t,
                            /* ************** */  /* *********************** */
                              &ZBL, /**/ &ZBR,       Z20,  z21,    /**/ Z22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ************************ */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &STL, /**/ &STR,       S00,  s01,     /**/ S02,
                                                     s10t, sigma11, /**/ s12t,
                            /* ************** */  /* ************************ */
                              &SBL, /**/ &SBR,       S20,  s21,     /**/ S22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &uT,                u01, 
                                                  upsilon11, 
                            /* ** */           /* ***** */
                              &uB,                u21,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &vT,                v01, 
                                                  nu11, 
                            /* ** */           /* ***** */
                              &vB,                v21,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &dT,                d0,
                                                  delta1,
                            /* ** */           /* ****** */
                              &dB,                d2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &eT,                e0,
                                                  epsilon1,
                            /* ** */           /* ******** */
                              &eB,                e2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &fT,                f0,
                                                  phi1,
                            /* ** */           /* **** */
                              &fB,                f2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &gT,                g0,
                                                  ghi1,
                            /* ** */           /* **** */
                              &gB,                g2,     FLA_TOP );
  }

  FLA_Obj_free( &last_elem );
  FLA_Obj_free( &beta );
  FLA_Obj_free( &minus_upsilon11 );
  FLA_Obj_free( &minus_zeta11 );
  FLA_Obj_free( &u );
  FLA_Obj_free( &v );
  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
  FLA_Obj_free( &f );
  FLA_Obj_free( &g );

  return FLA_SUCCESS;
}

