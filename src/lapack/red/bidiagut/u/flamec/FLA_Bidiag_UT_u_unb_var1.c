
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_unb_var1( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  return FLA_Bidiag_UT_u_step_unb_var1( A, TU, TV );
}

FLA_Error FLA_Bidiag_UT_u_step_unb_var1( FLA_Obj A, FLA_Obj T, FLA_Obj S )
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
  FLA_Obj  vT,              v01,
           vB,              nu11,
                            v21;
  FLA_Obj  v;

  FLA_Obj  a12t_l, a12t_r;
  FLA_Obj  A22_l, A22_r;
  FLA_Obj  v21_t,
           v21_b;

  FLA_Datatype datatype_A;
  dim_t        n_A;
  dim_t        b_alg;


  b_alg      = FLA_Obj_length( T );

  datatype_A = FLA_Obj_datatype( A );
  n_A        = FLA_Obj_width( A );

  FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );
  FLA_Part_2x2( S,    &STL, &STR,
                      &SBL, &SBR,     0, 0, FLA_TL );
  FLA_Part_2x1( v,    &vT, 
                      &vB,            0, FLA_TOP );

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

    /*------------------------------------------------------------*/

    // [ alpha11_new, u21, tau11 ] = House2( alpha11, a21 );
    FLA_Househ2_UT( FLA_LEFT,
                    alpha11,
                    a21, tau11 );

    if ( FLA_Obj_width( A22 ) > 0 )
    {
      FLA_Part_1x2( a12t,    &a12t_l, &a12t_r,      1, FLA_LEFT );
      FLA_Part_1x2( A22,     &A22_l,  &A22_r,       1, FLA_LEFT );
      FLA_Part_2x1( v21,     &v21_t, 
                             &v21_b,            1, FLA_TOP );

      // Apply H from the left to a12t and A22.
      FLA_Apply_H2_UT( FLA_LEFT, tau11, a21, a12t, A22 );

      // [ alpha12t, u12t_r, tau11 ] = House2( a12t_l, a12t_r );
      FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );

      // v21_t = 1;
      // v21_b = a12t_r;
      FLA_Set( FLA_ONE, v21_t );
      FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );

      // Apply H from the right to A22.
      FLA_Apply_H2_UT( FLA_RIGHT, sigma11, v21_b, A22_l, A22_r );

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
  }

  FLA_Obj_free( &v );

  return FLA_SUCCESS;
}

