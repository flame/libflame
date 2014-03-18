
#include "FLAME.h"

FLA_Error FLA_QR_UT_unb_var1( FLA_Obj A, FLA_Obj t )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj tLt,   tRt,      t0t,  tau1,  t2t;


  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_1x2( t,    &tLt,  &tRt,      0, FLA_LEFT );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_1x2_to_1x3( tLt,  /**/ tRt,      &t0t, /**/ &tau1, &t2t,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    // Compute tau11 and u21 from alpha11 and a21 such that tau11 and u21
    // determine a Householder transform H such that applying H from the
    // left to the column vector consisting of alpha11 and a21 annihilates
    // the entries in a21 (and updates alpha11).
    FLA_Househ2_UT( FLA_LEFT,
                    alpha11,
                    a21, tau1 );

    // / a12t \ =  H / a12t \
    // \ A22  /      \ A22  /
    //
    // where H is formed from tau11 and u21.
    FLA_Apply_H2_UT( FLA_LEFT, tau1, a21, a12t,
                                          A22 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &tLt,  /**/ &tRt,      t0t, tau1, /**/ t2t,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

