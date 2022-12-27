/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_Q2_UT_lhfc_blk_var1( FLA_Obj D, FLA_Obj T, FLA_Obj W1, FLA_Obj C, 
                                                                           FLA_Obj E, fla_apq2ut_t* cntl )
{
  FLA_Obj DL,    DR,       D0,  D1,  D2;

  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj T1T,
          T2B;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  FLA_Obj W1TL,  W1TR,
          W1BL,  W1BR;

  dim_t   b_alg, b;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_TOP );

  while ( FLA_Obj_width( DL ) < FLA_Obj_width( D ) ){

    b = fla_min( b_alg, FLA_Obj_width( DR ) );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &T2,
                           b, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                        /* ** */            /* ** */
                                              &C1, 
                           CB,                &C2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,    &T1T, 
                         &T2B,     b, FLA_TOP );

    FLA_Part_2x2( W1,    &W1TL, &W1TR,
                         &W1BL, &W1BR,     b, FLA_Obj_width( C1 ), FLA_TL );

    // W1TL = C1;

    FLA_Copyt_internal( FLA_NO_TRANSPOSE, C1, W1TL,
                        FLA_Cntl_sub_copyt( cntl ) );

    // W1TL = inv( triu( T1T ) )' * ( C1 + D1' * E );

    FLA_Gemm_internal( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_ONE, D1, E, FLA_ONE, W1TL,
                       FLA_Cntl_sub_gemm1( cntl ) );

    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T1T, W1TL,
                       FLA_Cntl_sub_trsm( cntl ) );

    // C1 = C1 - W1TL;
    // E  = E  - D1 * W1TL;

    FLA_Axpyt_internal( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, W1TL, C1,
                        FLA_Cntl_sub_axpyt( cntl ) );

    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, D1, W1TL, FLA_ONE, E,
                       FLA_Cntl_sub_gemm2( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ T2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                                                  C1, 
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

