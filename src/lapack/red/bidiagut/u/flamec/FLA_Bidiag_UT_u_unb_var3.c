
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_unb_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  return FLA_Bidiag_UT_u_step_unb_var3( A, TU, TV );
}

FLA_Error FLA_Bidiag_UT_u_step_unb_var3( FLA_Obj A, FLA_Obj T, FLA_Obj S )
{
  FLA_Obj  ATL,   ATR,      A00,  a01,     A02, 
           ABL,   ABR,      a10t, alpha11, a12t,
                            A20,  a21,     A22;
  FLA_Obj  TTL,   TTR,      T00,  t01,   T02, 
           TBL,   TBR,      t10t, tau11, t12t,
                            T20,  t21,   T22;
  FLA_Obj  STL,   STR,      S00,  s01,     S02, 
           SBL,   SBR,      s10t, sigma11, s12t,
                            S20,  s21,     S22;
  FLA_Obj  wT,              w01,
           wB,              omega11,
                            w21;
  FLA_Obj  apT,             a01p,
           apB,             alpha11p,
                            a12p;
  FLA_Obj  uT,              u01,
           uB,              upsilon11,
                            u21;
  FLA_Obj  uTp,             u01p,
           uBp,             upsilon11p,
                            u21p;
  FLA_Obj  vT,              v01,
           vB,              nu11,
                            v21;
  FLA_Obj  yT,              y01,
           yB,              psi11,
                            y21;
  FLA_Obj  zT,              z01,
           zB,              zeta11,
                            z21;
  FLA_Obj  w, ap, u, up, v, y, z;

  FLA_Obj  minus_inv_tau11;
  FLA_Obj  beta;
  FLA_Obj  alpha12;
  FLA_Obj  minus_conj_alpha12;
  FLA_Obj  psi11_minus_alpha12;
  FLA_Obj  minus_upsilon11;
  FLA_Obj  minus_conj_nu11;
  FLA_Obj  minus_conj_psi11;
  FLA_Obj  minus_zeta11;

  FLA_Obj  a12t_l, a12t_r;
  FLA_Obj  a12p_t,
           a12p_b;
  FLA_Obj  A22_l, A22_r;
  FLA_Obj  v21_t,
           v21_b;

  FLA_Datatype datatype_A;
  dim_t        m_A, n_A;
  dim_t        b_alg;


  b_alg      = FLA_Obj_length( T );

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );
  n_A        = FLA_Obj_width( A );

  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_inv_tau11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &beta );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &alpha12 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_conj_alpha12 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &psi11_minus_alpha12 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_upsilon11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_conj_nu11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_conj_psi11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_zeta11 );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &ap );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &up );
  FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );
  FLA_Part_2x2( S,    &STL, &STR,
                      &SBL, &SBR,     0, 0, FLA_TL );
  FLA_Part_2x1( w,    &wT, 
                      &wB,            0, FLA_TOP );
  FLA_Part_2x1( ap,   &apT, 
                      &apB,           0, FLA_TOP );
  FLA_Part_2x1( u,    &uT, 
                      &uB,            0, FLA_TOP );
  FLA_Part_2x1( up,   &uTp, 
                      &uBp,           0, FLA_TOP );
  FLA_Part_2x1( v,    &vT, 
                      &vB,            0, FLA_TOP );
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
    FLA_Repart_2x2_to_3x3( STL, /**/ STR,       &S00,  /**/ &s01,     &S02,
                        /* ************* */   /* ************************** */
                                                &s10t, /**/ &sigma11, &s12t,
                           SBL, /**/ SBR,       &S20,  /**/ &s21,     &S22,
                           1, 1, FLA_BR );
    FLA_Repart_2x1_to_3x1( wT,                &w01, 
                        /* ** */            /* ***** */
                                              &omega11, 
                           wB,                &w21,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( apT,               &a01p, 
                        /* ** */            /* ***** */
                                              &alpha11p, 
                           apB,               &a12p,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( uT,                &u01, 
                        /* ** */            /* ***** */
                                              &upsilon11, 
                           uB,                &u21,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( uTp,               &u01p, 
                        /* ** */            /* ***** */
                                              &upsilon11p, 
                           uBp,               &u21p,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( vT,                &v01, 
                        /* ** */            /* ***** */
                                              &nu11, 
                           vB,                &v21,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( yT,                &y01, 
                        /* ** */            /* ***** */
                                              &psi11, 
                           yB,                &y21,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( zT,                &z01, 
                        /* ** */            /* ***** */
                                              &zeta11, 
                           zB,                &z21,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    if ( FLA_Obj_length( ATL ) > 0 )
    {
      FLA_Copy( upsilon11, minus_upsilon11 );
      FLA_Scal( FLA_MINUS_ONE, minus_upsilon11 );

      FLA_Copy( zeta11, minus_zeta11 );
      FLA_Scal( FLA_MINUS_ONE, minus_zeta11 );

      FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, psi11, minus_conj_psi11 );
      FLA_Scal( FLA_MINUS_ONE, minus_conj_psi11 );

      FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, nu11, minus_conj_nu11 );
      FLA_Scal( FLA_MINUS_ONE, minus_conj_nu11 );

      // alpha11 = alpha11 - upsilon11 * conj(psi11) - zeta11 * conj(nu1);
      FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi11, upsilon11, alpha11 );
      FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_nu11,  zeta11,    alpha11 );

      // a21 = a21 - u21 * conj(psi11) - z21 * conj(nu11);
      FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi11, u21, a21 );
      FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_nu11,  z21, a21 );

      // a12t = a12t - upsilon11 * y21' - zeta11 * v21';
      FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_upsilon11, y21, a12t );
      FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_zeta11,    v21, a12t );
    }

    // [ alpha11, u21p, tau11 ] = House2( alpha11, a21 );
    FLA_Househ2_UT( FLA_LEFT,
                    alpha11,
                    a21, tau11 );
    FLA_Copy( a21, u21p );

    if ( FLA_Obj_width( A22 ) > 0 )
    {
      // minus_inv_tau11 = - 1 / tau11;
      FLA_Copy( FLA_MINUS_ONE, minus_inv_tau11 );
      FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, minus_inv_tau11 );

      // a12p = ( tau11 - 1 ) * a12t^T / tau11;
      //      = a12t^T - ( 1 / tau11 ) * a12t^T;
      FLA_Copyt( FLA_TRANSPOSE, a12t, a12p );
      FLA_Axpyt( FLA_TRANSPOSE, minus_inv_tau11, a12t, a12p );
    }

    if ( FLA_Obj_length( ATL ) > 0 )
    {
      // A22 = A22 - u21 * y21' - z21 * v21';
      FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u21, y21, A22 );
      FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );
    }

    if ( FLA_Obj_width( A22 ) > 0 )
    {
      // y21 = A22' * u21p;
      FLA_Gemvc( FLA_CONJ_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, u21p, FLA_ZERO, y21 );

      // a12p = a12p - conj(y21) / tau11;
      FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );

      // w21 = A22 * conj(a12p);
      FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );

      // y21 = y21 + conj(a12t)^T;
      FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_ONE, a12t, y21 );

      FLA_Part_1x2( a12t,    &a12t_l, &a12t_r,      1, FLA_LEFT );
      FLA_Part_2x1( v21,     &v21_t, 
                             &v21_b,            1, FLA_TOP );
      FLA_Part_2x1( a12p,    &a12p_t, 
                             &a12p_b,           1, FLA_TOP );

      // [ alpha12, psi11_minus_alpha12, sigma11 ] = House2s( a12p_t, a12p_b );
      FLA_Househ2s_UT( FLA_RIGHT,
                       a12p_t,
                       a12p_b,
                       alpha12, psi11_minus_alpha12, sigma11 );

      // v21 = conj( ( a12p - alpha12 * e0 ) / ( psi11 - alpha12 ) );
      FLA_Copy( a12p, v21 );
      FLA_Mult_add( FLA_MINUS_ONE, alpha12, v21_t );
      FLA_Inv_scalc( FLA_NO_CONJUGATE, psi11_minus_alpha12, v21 );
      FLA_Conjugate( v21_b );

      // a12t_l = alpha12;
      // a12t_r = v21_b^T;
      FLA_Copyt( FLA_NO_TRANSPOSE, alpha12, a12t_l );
      FLA_Copyt( FLA_TRANSPOSE, v21_b, a12t_r );
    }

    // u21 = u21p;
    FLA_Copy( u21p, u21 );

    if ( FLA_Obj_width( A22 ) > 0 )
    {
      // beta = - y21' * v21 / tau11;
      FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      FLA_Scal( FLA_MINUS_ONE, beta );
      FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, beta );

      FLA_Part_1x2( A22,    &A22_l, &A22_r,      1, FLA_LEFT );

      // minus_conj_alpha12 = - conj(alpha12);
      FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha12, minus_conj_alpha12 );
      FLA_Scal( FLA_MINUS_ONE, minus_conj_alpha12 );

      // z21 = ( w21 - conj(alpha12) * A22 * e0 ) / conj(psi11 - alpha12) + beta * u21;
      FLA_Copy( w21, z21 );
      FLA_Axpy( minus_conj_alpha12, A22_l, z21 );
      FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, z21 );
      FLA_Axpy( beta, u21, z21 );

      // y21 = y21 / tau11;
      // z21 = z21 / sigma11;
      FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11,   y21 );
      FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );

      // s01 = conj(V02) * v21;
      FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, s01 );
    }

    // t01 = a10t' + U20' * u21;
    FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
    FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ONE, t01 );

    // Update A22 if this is the last iteration; this is needed when we're
    // being called from the blocked routine so A22 is left in a valid state.
    if ( FLA_Obj_length( ATL ) + 1 == b_alg &&
         FLA_Obj_width( A22 ) > 0 )
    {
      // A22 = A22 - u21 * y21' - z21 * v21';
      FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u21, y21, A22 );
      FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );
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
    FLA_Cont_with_3x3_to_2x2( &STL, /**/ &STR,       S00,  s01,     /**/ S02,
                                                     s10t, sigma11, /**/ s12t,
                            /* ************** */  /* ************************ */
                              &SBL, /**/ &SBR,       S20,  s21,     /**/ S22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &wT,                w01, 
                                                  omega11, 
                            /* ** */           /* ***** */
                              &wB,                w21,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &apT,               a01p, 
                                                  alpha11p, 
                            /* ** */           /* ***** */
                              &apB,               a12p,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &uT,                u01, 
                                                  upsilon11, 
                            /* ** */           /* ***** */
                              &uB,                u21,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &uTp,               u01p, 
                                                  upsilon11p, 
                            /* ** */           /* ***** */
                              &uBp,               u21p,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &vT,                v01, 
                                                  nu11, 
                            /* ** */           /* ***** */
                              &vB,                v21,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &yT,                y01, 
                                                  psi11, 
                            /* ** */           /* ***** */
                              &yB,                y21,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &zT,                z01, 
                                                  zeta11, 
                            /* ** */           /* ***** */
                              &zB,                z21,     FLA_TOP );
  }

  FLA_Obj_free( &minus_inv_tau11 );
  FLA_Obj_free( &beta );
  FLA_Obj_free( &alpha12 );
  FLA_Obj_free( &minus_conj_alpha12 );
  FLA_Obj_free( &psi11_minus_alpha12 );
  FLA_Obj_free( &minus_upsilon11 );
  FLA_Obj_free( &minus_conj_nu11 );
  FLA_Obj_free( &minus_conj_psi11 );
  FLA_Obj_free( &minus_zeta11 );
  FLA_Obj_free( &w );
  FLA_Obj_free( &ap );
  FLA_Obj_free( &u );
  FLA_Obj_free( &up );
  FLA_Obj_free( &v );
  FLA_Obj_free( &y );
  FLA_Obj_free( &z );

  return FLA_SUCCESS;
}

