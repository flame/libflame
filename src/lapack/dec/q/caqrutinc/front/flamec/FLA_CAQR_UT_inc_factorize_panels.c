/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_CAQR_UT_inc_factorize_panels( dim_t nb_part, FLA_Obj A, FLA_Obj TW )
{
  FLA_Obj AT,              A0, 
          AB,              A1,
                           A2;

  FLA_Obj TWT,             TW0, 
          TWB,             TW1,
                           TW2;

  dim_t b;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( TW,   &TWT, 
                      &TWB,           0, FLA_TOP );

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

    /*------------------------------------------------------------*/

    // Perform an incremental QR factorization on A1, writing triangular
    // block Householder factors to T in TW1.
    FLASH_QR_UT_inc( A1, TW1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,               A0, 
                                                 A1, 
                            /* ** */          /* ** */
                              &AB,               A2,      FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &TWT,              TW0, 
                                                 TW1, 
                            /* ** */          /* ** */
                              &TWB,              TW2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

