/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_l_blk_var2( FLA_Obj A, FLA_Obj T )
{
  FLA_Obj  ATL,   ATR,      A00, A01, A02, 
           ABL,   ABR,      A10, A11, A12,
                            A20, A21, A22;
  FLA_Obj  TL,    TR,       T0, T1, T2; 

  FLA_Obj  T1_tl;
  FLA_Obj  none, none2, none3;
  dim_t    b_alg, b;

  b_alg = FLA_Obj_length( T );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT ); 

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) )
  {
    b = min( FLA_Obj_length( ABR ), b_alg );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );
    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &T2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x2( T1,     &T1_tl, &none,   
                          &none2, &none3,   b, b, FLA_TL ); 

    // [ ABR, T1 ] = FLA_Tridiag_UT_l_step_unb_var2( ABR, T1, b );
    //FLA_Tridiag_UT_l_step_unb_var2( ABR, T1_tl );
    //FLA_Tridiag_UT_l_step_ofu_var2( ABR, T1_tl );
    FLA_Tridiag_UT_l_step_opt_var2( ABR, T1_tl );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ T2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

