/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Accum_T_UT_fc_blk_var2( FLA_Obj A, FLA_Obj t, FLA_Obj T )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj tT,              t0,
          tB,              t1,
                           t2;

  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj AB1;

  dim_t b_alg, b;

  b_alg = FLA_Obj_length( T );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( t,    &tT, 
                      &tB,            0, FLA_TOP );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  while ( FLA_Obj_length( tB ) > 0 ) {

    b = fla_min( FLA_Obj_length( tB ), b_alg );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x1_to_3x1( tT,                &t0, 
                        /* ** */            /* ** */
                                              &t1, 
                           tB,                &t2,        b, FLA_BOTTOM );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &T2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Merge_2x1( A11,
                   A21,   &AB1 );

    //FLA_Accum_T_UT_fc_unb_var1( AB1, t1, T1 );
    FLA_Accum_T_UT_fc_opt_var1( AB1, t1, T1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &tT,                t0, 
                                                  t1, 
                            /* ** */           /* ** */
                              &tB,                t2,     FLA_TOP );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ T2,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}

