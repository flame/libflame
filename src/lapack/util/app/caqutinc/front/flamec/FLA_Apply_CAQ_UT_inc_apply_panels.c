/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_CAQ_UT_inc_apply_panels( dim_t nb_part, FLA_Obj A, FLA_Obj TW, FLA_Obj W, FLA_Obj B )
{
  FLA_Obj AT,              A0, 
          AB,              A1,
                           A2;

  FLA_Obj TWT,             TW0, 
          TWB,             TW1,
                           TW2;

  FLA_Obj WT,              W0,
          WB,              W1,
                           W2;

  FLA_Obj BT,              B0, 
          BB,              B1,
                           B2;

  dim_t b;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( TW,   &TWT, 
                      &TWB,           0, FLA_TOP );

  FLA_Part_2x1( W,    &WT, 
                      &WB,            0, FLA_TOP );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_TOP );

  while ( FLA_Obj_length( AB ) > 0 ){

    b = min( nb_part, FLA_Obj_length( AB ) );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* ** */
                                              &A1, 
                           AB,                &A2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( TWT,               &TW0, 
                        /* ** */            /* ** */
                                              &TW1, 
                           TWB,               &TW2,       b, FLA_BOTTOM );

    // NOTE: we use a blocksize of 1 for W since it has exactly nb_part
    // rows (where each row is a row panels of b_alg x b_flash blocks).
    FLA_Repart_2x1_to_3x1( WT,                &W0, 
                        /* ** */            /* ** */
                                              &W1, 
                           WB,                &W2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                        /* ** */            /* ** */
                                              &B1, 
                           BB,                &B2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // Apply incremental Q's associated with each block A1 to the
    // corresponding block of right-hand side B1.
    FLASH_Apply_Q_UT_inc( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                          A1, TW1, W1, B1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,               A0, 
                                                 A1, 
                            /* ** */          /* ** */
                              &AB,               A2,      FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &TWT,              TW0, 
                                                 TW1, 
                            /* ** */          /* ** */
                              &TWB,              TW2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &WT,                W0, 
                                                  W1, 
                            /* ** */           /* ** */
                              &WB,                W2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                                                  B1, 
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

