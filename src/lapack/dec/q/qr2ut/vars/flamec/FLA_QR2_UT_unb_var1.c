
#include "FLAME.h"

FLA_Error FLA_QR2_UT_unb_var1( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T )
{
  FLA_Obj UTL,   UTR,      U00,  u01,       U02, 
          UBL,   UBR,      u10t, upsilon11, u12t,
                           U20,  u21,       U22;

  FLA_Obj DL,    DR,       D0,  d1,  D2;

  FLA_Obj TTL,   TTR,      T00,  t01,   T02, 
          TBL,   TBR,      t10t, tau11, t12t,
                           T20,  t21,   T22;


  FLA_Part_2x2( U,    &UTL, &UTR,
                      &UBL, &UBR,     0, 0, FLA_TL );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );

  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  while ( FLA_Obj_min_dim( UBR ) > 0 ){

    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00,  /**/ &u01,       &U02,
                        /* ************* */   /* ************************** */
                                                &u10t, /**/ &upsilon11, &u12t,
                           UBL, /**/ UBR,       &U20,  /**/ &u21,       &U22,
                           1, 1, FLA_BR );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &d1, &D2,
                           1, FLA_RIGHT );

    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************ */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/

    // Compute tau11 and u2 from upsilon11 and d1 such that tau11 and u2
    // determine a Householder transform H such that applying H from the
    // left to the column vector consisting of upsilon11 and d1 annihilates
    // the entries in d1 (and updates upsilon11).
    FLA_Househ2_UT( FLA_LEFT,
                    upsilon11,
                    d1, tau11 );

    // / u12t \ =  H / u12t \
    // \  D2  /      \  D2  / 
    //
    // where H is formed from tau11 and d1.
    FLA_Apply_H2_UT( FLA_LEFT, tau11, d1, u12t,
                                          D2 );

    // t01 = D0' * d1;
    FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, D0, d1, FLA_ZERO, t01 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00,  u01,       /**/ U02,
                                                     u10t, upsilon11, /**/ u12t,
                            /* ************** */  /* ************************ */
                              &UBL, /**/ &UBR,       U20,  u21,       /**/ U22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, d1, /**/ D2,
                              FLA_LEFT );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ********************** */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );
  }

  return FLA_SUCCESS;
}

