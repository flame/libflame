
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_unb_var2( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  return FLA_Bidiag_UT_u_step_unb_var2( A, TU, TV );
}

FLA_Error FLA_Bidiag_UT_u_step_unb_var2( FLA_Obj A, FLA_Obj T, FLA_Obj S )
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
  FLA_Obj  yT,              y01,
           yB,              psi11,
                            y21;
  FLA_Obj  zT,              z01,
           zB,              zeta11,
                            z21;
  FLA_Obj  vT,              v01,
           vB,              nu11,
                            v21;
  FLA_Obj  v, y, z;

  FLA_Obj  beta;

  FLA_Obj  a12t_l, a12t_r;
  FLA_Obj  v21_t,
           v21_b;

  FLA_Datatype datatype_A;
  dim_t        m_A, n_A;
  dim_t        b_alg;


  b_alg      = FLA_Obj_length( T );

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );
  n_A        = FLA_Obj_width( A );

  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &beta );
  FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &y );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );
  FLA_Part_2x2( S,    &STL, &STR,
                      &SBL, &SBR,     0, 0, FLA_TL );
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

    // [ alpha11_new, u21, tau11 ] = House2( alpha11, a21 );
    FLA_Househ2_UT( FLA_LEFT,
                    alpha11,
                    a21, tau11 );

    if ( FLA_Obj_width( A22 ) > 0 )
    {
      // y21' = a12t + u21' * A22;
      // y21  = conj(a12t) + A22' * u21;
      FLA_Copyt( FLA_CONJ_TRANSPOSE, a12t, y21 );
      FLA_Gemvc( FLA_CONJ_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, a21, FLA_ONE, y21 );

      // y21 = y21 / tau11;
      FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, y21 );

      // a12t = a12t - conj(y21)^T;
      FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_MINUS_ONE, y21, a12t );

      FLA_Part_1x2( a12t,    &a12t_l, &a12t_r,      1, FLA_LEFT );
      FLA_Part_2x1( v21,     &v21_t, 
                             &v21_b,            1, FLA_TOP );

      // [ a12t_l, v12t_b, sigma11 ] = House2( a12t_l, a12t_r );
      FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );

      // v21_t = 1;
      // v21_b = a12t_r^T;
      FLA_Set( FLA_ONE, v21_t );
      FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );

      // beta = - y21' * v21;
      FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      FLA_Scal( FLA_MINUS_ONE, beta );

      // z21 = ( A22 - u21 * y21' ) * v21 / sigma11;
      //     = ( A22 * v21 - u21 * y21' * v21 ) / sigma11;
      //     = ( A22 * v21 + beta * u21 ) / sigma11;
      FLA_Copy( a21, z21 );
      FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_NO_CONJUGATE, FLA_ONE, A22, v21, beta, z21 );
      FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );

      // A22 = A22 - u21 * y21' - z21 * v21';
      FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, a21, y21, A22 );
      FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z21, v21, A22 );

      // s01 = conj(V02) * v21;
      FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, s01 );
    }

    // t01 = a10t' + U20' * u21;
    FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
    FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );

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

  FLA_Obj_free( &beta );
  FLA_Obj_free( &v );
  FLA_Obj_free( &y );
  FLA_Obj_free( &z );

  return FLA_SUCCESS;
}

