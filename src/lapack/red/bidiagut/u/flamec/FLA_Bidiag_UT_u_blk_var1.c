/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_blk_var1( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  FLA_Obj  ATL,   ATR,      A00, A01, A02, 
           ABL,   ABR,      A10, A11, A12,
                            A20, A21, A22;
  FLA_Obj  TUL,   TUR,      TU0, TU1, TU2; 
  FLA_Obj  TVL,   TVR,      TV0, TV1, TV2; 

  FLA_Obj  TU1_tl;
  FLA_Obj  TV1_tl;
  FLA_Obj  none, none2, none3;
  dim_t    b_alg, b;

  b_alg = FLA_Obj_length( TU );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,   0, 0, FLA_TL );
  FLA_Part_1x2( TU,   &TUL, &TUR,      0, FLA_LEFT ); 
  FLA_Part_1x2( TV,   &TVL, &TVR,      0, FLA_LEFT ); 

  while ( FLA_Obj_min_dim( ABR ) > 0 )
  {
    b = min( FLA_Obj_min_dim( ABR ), b_alg );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );
    FLA_Repart_1x2_to_1x3( TUL, /**/ TUR,       &TU0, /**/ &TU1, &TU2,
                           b, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( TVL, /**/ TVR,       &TV0, /**/ &TV1, &TV2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x2( TU1,     &TU1_tl, &none,   
                           &none2,  &none3,   b, b, FLA_TL ); 

    FLA_Part_2x2( TV1,     &TV1_tl, &none,   
                           &none2,  &none3,   b, b, FLA_TL ); 

    // [ ABR, T1 ] = FLA_Bidiag_UT_u_step_unb_var1( ABR, TU1, TV1, b );
    //FLA_Bidiag_UT_u_step_unb_var1( ABR, TU1_tl, TV1_tl );
    FLA_Bidiag_UT_u_step_opt_var1( ABR, TU1_tl, TV1_tl );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
    FLA_Cont_with_1x3_to_1x2( &TUL, /**/ &TUR,       TU0, TU1, /**/ TU2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &TVL, /**/ &TVR,       TV0, TV1, /**/ TV2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

