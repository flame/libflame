
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_unb_var4( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
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

  r_val = FLA_Bidiag_UT_u_step_unb_var4( A, Y, Z, TU, TV );

  FLA_Obj_free( &Y );
  FLA_Obj_free( &Z );

  return r_val;
}

FLA_Error FLA_Bidiag_UT_u_step_unb_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj S )
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
  FLA_Obj  wT,              w01,
           wB,              omega11,
                            w21;
  FLA_Obj  alT,             a01l,
           alB,             alpha11l,
                            a22l;
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
  FLA_Obj  w, al, ap, u, up, v;
  FLA_Obj  d, e, f, g;

  FLA_Obj  minus_inv_tau11;
  FLA_Obj  last_elem;
  FLA_Obj  beta;
  FLA_Obj  alpha12;
  FLA_Obj  minus_alpha12;
  FLA_Obj  minus_conj_alpha12;
  FLA_Obj  psi11_minus_alpha12;
  FLA_Obj  minus_upsilon11;
  FLA_Obj  minus_conj_nu11;
  FLA_Obj  minus_conj_psi11;
  FLA_Obj  minus_zeta11;

  FLA_Obj  a01_t,
           a01_b;
  FLA_Obj  A02_l, A02_r;
  FLA_Obj  a12t_l, a12t_r;
  FLA_Obj  a12p_t,
           a12p_b;
  FLA_Obj  A22_l, A22_r;
  FLA_Obj  v21_t,
           v21_b;
  FLA_Obj  Y20_t,
           Y20_b;
  FLA_Obj  a2;

  FLA_Datatype datatype_A;
  dim_t        m_A, n_A;
  dim_t        b_alg;


  b_alg      = FLA_Obj_length( T );

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );
  n_A        = FLA_Obj_width( A );

  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_inv_tau11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &last_elem );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &beta );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &alpha12 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_alpha12 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_conj_alpha12 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &psi11_minus_alpha12 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_upsilon11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_conj_nu11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_conj_psi11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_zeta11 );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &al );
  FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &ap );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &up );
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
  FLA_Part_2x1( w,    &wT, 
                      &wB,            0, FLA_TOP );
  FLA_Part_2x1( al,   &alT, 
                      &alB,           0, FLA_TOP );
  FLA_Part_2x1( ap,   &apT, 
                      &apB,           0, FLA_TOP );
  FLA_Part_2x1( u,    &uT, 
                      &uB,            0, FLA_TOP );
  FLA_Part_2x1( up,   &uTp, 
                      &uBp,           0, FLA_TOP );
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
    FLA_Repart_2x1_to_3x1( wT,                &w01, 
                        /* ** */            /* ***** */
                                              &omega11, 
                           wB,                &w21,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( alT,               &a01l, 
                        /* ** */            /* ***** */
                                              &alpha11l, 
                           alB,               &a22l,        1, FLA_BOTTOM );
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

      // y21 = - Y20 * ( U20' * u21p ) - V20 * ( Z20' * u21p );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21p, FLA_ZERO, d0 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, u21p, FLA_ZERO, e0 );

      FLA_Set( FLA_ZERO, y21 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      FLA_Gemv( FLA_TRANSPOSE,    FLA_MINUS_ONE, A02, e0, FLA_ONE, y21 );

      // t01 = a10t' + U20' * u21;
      FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      FLA_Axpy( FLA_ONE, d0, t01 );

      // y21 = y21 + A22' * u21p;
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, u21p, FLA_ONE, y21 );

      // a12p = a12p - conj(y21) / tau11;
      FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_inv_tau11, y21, a12p );

      // w21 = A22 * conj(a12p);
      FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A22, a12p, FLA_ZERO, w21 );

      // w21 = w21 - U20 * ( Y20' * conj(a12p) ) - Z20 * ( V20' * conj(a12p) );
      FLA_Gemvc( FLA_CONJ_TRANSPOSE,    FLA_CONJUGATE, FLA_ONE, Y20, a12p, FLA_ZERO, f0 );
      FLA_Gemvc( FLA_CONJ_NO_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A02, a12p, FLA_ZERO, g0 );

      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, w21 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, g0, FLA_ONE, w21 );

      FLA_Part_1x2( A22,    &A22_l, &A22_r,     1, FLA_LEFT );
      FLA_Part_2x1( Y20,    &Y20_t,
                            &Y20_b,             1, FLA_TOP );
      FLA_Part_1x2( A02,    &A02_l, &A02_r,     1, FLA_LEFT );

      // a22l = A22 * e0 - U20 * ( Y20' * e0 ) - Z20 * ( V20' * e0 );
      FLA_Copy( A22_l, a22l );
      FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A20, Y20_t, FLA_ONE, a22l );
      FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, Z20, A02_l, FLA_ONE, a22l );

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

      // minus_conj_alpha12 = - conj(alpha12);
      FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha12, minus_conj_alpha12 );
      FLA_Scal( FLA_MINUS_ONE, minus_conj_alpha12 );

      // s01 = V20' * v21;
      //     = conj(V02) * v21;
      //     = conj(V02) * conj( ( a12p - alpha12 * e0 ) / ( psi11 - alpha12 ) );
      //     = conj(V02) * ( conj(a12p) - conj(alpha12) * e0 ) / conj( psi11 - alpha12 ) );
      //     = ( conj(V02) * conj(a12p) - conj(V02) * conj(alpha12) * e0 ) / conj( psi11 - alpha12 );
      //     = ( g0 - conj(V02) * conj(alpha12) * e0 ) / conj( psi11 - alpha12 );
      FLA_Copy( g0, s01 );
      FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_conj_alpha12, A02_l, s01 );
      FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, s01 );

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

      // z21 = ( w21 - conj(alpha12) * a22l ) / conj(psi11 - alpha12) + beta * u21;
      FLA_Copy( w21, z21 );
      FLA_Axpy( minus_conj_alpha12, a22l, z21 );
      FLA_Inv_scalc( FLA_CONJUGATE, psi11_minus_alpha12, z21 );
      FLA_Axpy( beta, u21, z21 );

      // y21 = y21 / tau11;
      // z21 = z21 / sigma11;
      FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11,   y21 );
      FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
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
    FLA_Cont_with_3x1_to_2x1( &wT,                w01, 
                                                  omega11, 
                            /* ** */           /* ***** */
                              &wB,                w21,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &alT,               a01l, 
                                                  alpha11l, 
                            /* ** */           /* ***** */
                              &alB,               a22l,     FLA_TOP );
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

  FLA_Obj_free( &minus_inv_tau11 );
  FLA_Obj_free( &last_elem );
  FLA_Obj_free( &beta );
  FLA_Obj_free( &alpha12 );
  FLA_Obj_free( &minus_alpha12 );
  FLA_Obj_free( &minus_conj_alpha12 );
  FLA_Obj_free( &psi11_minus_alpha12 );
  FLA_Obj_free( &minus_upsilon11 );
  FLA_Obj_free( &minus_conj_nu11 );
  FLA_Obj_free( &minus_conj_psi11 );
  FLA_Obj_free( &minus_zeta11 );
  FLA_Obj_free( &w );
  FLA_Obj_free( &al );
  FLA_Obj_free( &ap );
  FLA_Obj_free( &u );
  FLA_Obj_free( &up );
  FLA_Obj_free( &v );
  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
  FLA_Obj_free( &f );
  FLA_Obj_free( &g );

  return FLA_SUCCESS;
}

