/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_CAQR2_UT_unb_var1( FLA_Obj U,
                                 FLA_Obj D, FLA_Obj T )
{
  FLA_Obj UTL,   UTR,      U00,  u01,       U02, 
          UBL,   UBR,      u10t, upsilon11, u12t,
                           U20,  u21,       U22;

  FLA_Obj DTL,   DTR,      D00,  d01,     D02, 
          DBL,   DBR,      d10t, delta11, d12t,
                           D20,  d21,     D22;

  FLA_Obj TTL,   TTR,      T00,  t01,   T02, 
          TBL,   TBR,      t10t, tau11, t12t,
                           T20,  t21,   T22;

  FLA_Obj d1, D2;

  FLA_Obj d01T,
          d01B;

  FLA_Obj D00T,
          D00B;

  dim_t   m_DT;

  // Begin partitioning diagonally through D with m - n rows above
  // the diagonal.
  m_DT = FLA_Obj_length( D ) - FLA_Obj_width( D );

  FLA_Part_2x2( U,    &UTL, &UTR,
                      &UBL, &UBR,     0, 0, FLA_TL );

  FLA_Part_2x2( D,    &DTL, &DTR,
                      &DBL, &DBR,     m_DT, 0, FLA_TL );

  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  while ( FLA_Obj_min_dim( UBR ) > 0 ){

    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00,  /**/ &u01,       &U02,
                        /* ************* */   /* ************************** */
                                                &u10t, /**/ &upsilon11, &u12t,
                           UBL, /**/ UBR,       &U20,  /**/ &u21,       &U22,
                           1, 1, FLA_BR );

    FLA_Repart_2x2_to_3x3( DTL, /**/ DTR,       &D00,  /**/ &d01,     &D02,
                        /* ************* */   /* ************************** */
                                                &d10t, /**/ &delta11, &d12t,
                           DBL, /**/ DBR,       &D20,  /**/ &d21,     &D22,
                           1, 1, FLA_BR );

    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************ */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/

    FLA_Merge_2x1( d01,
                   delta11,  &d1 );

    FLA_Merge_2x1( D02,
                   d12t,     &D2 );

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

    FLA_Part_2x1( d01,   &d01T,
                         &d01B,    m_DT, FLA_TOP );

    FLA_Part_2x1( D00,   &D00T,
                         &D00B,    m_DT, FLA_TOP );

    // t01 = D00' * d01;
    //     = D00T' * d01T + triu( D00B )' * d01B;
    FLA_Copy_external( d01B, t01 );
	FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       D00B, t01 );
    FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, D00T, d01T, FLA_ONE, t01 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00,  u01,       /**/ U02,
                                                     u10t, upsilon11, /**/ u12t,
                            /* ************** */  /* ************************ */
                              &UBL, /**/ &UBR,       U20,  u21,       /**/ U22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &DTL, /**/ &DTR,       D00,  d01,     /**/ D02,
                                                     d10t, delta11, /**/ d12t,
                            /* ************** */  /* ************************ */
                              &DBL, /**/ &DBR,       D20,  d21,     /**/ D22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ********************** */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );
  }

  return FLA_SUCCESS;
}

