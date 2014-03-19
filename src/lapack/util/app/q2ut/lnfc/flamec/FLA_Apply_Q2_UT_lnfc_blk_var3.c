/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_Q2_UT_lnfc_blk_var3( FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, 
                                                                          FLA_Obj E, fla_apq2ut_t* cntl )
{
  FLA_Obj WL,    WR,       W0,  W1,  W2;

  FLA_Obj CL,    CR,       C0,  C1,  C2;

  FLA_Obj EL,    ER,       E0,  E1,  E2;

  dim_t b;

  FLA_Part_1x2( W,    &WL,  &WR,      0, FLA_RIGHT );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_RIGHT );

  FLA_Part_1x2( E,    &EL,  &ER,      0, FLA_RIGHT );

  while ( FLA_Obj_width( CR ) < FLA_Obj_width( C ) ){

    b = FLA_Determine_blocksize( CL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( WL,  /**/ WR,        &W0, &W1, /**/ &W2,
                           b, FLA_LEFT );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, &C1, /**/ &C2,
                           b, FLA_LEFT );

    FLA_Repart_1x2_to_1x3( EL,  /**/ ER,        &E0, &E1, /**/ &E2,
                           b, FLA_LEFT );

    /*------------------------------------------------------------*/

    //  / C1 \ =  Q / C1 \
    //  \ E1 /      \ E1 /
    //
    //  where Q is formed from D and T.

    FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                              D, T, W1, C1,
                                        E1, FLA_Cntl_sub_apq2ut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &WL,  /**/ &WR,        W0, /**/ W1, W2,
                              FLA_RIGHT );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, /**/ C1, C2,
                              FLA_RIGHT );

    FLA_Cont_with_1x3_to_1x2( &EL,  /**/ &ER,        E0, /**/ E1, E2,
                              FLA_RIGHT );
  }

  return FLA_SUCCESS;
}

