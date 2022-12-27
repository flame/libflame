/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_CAQ2_UT_lhfc_blk_var1( FLA_Obj D, FLA_Obj T, FLA_Obj W1, FLA_Obj C, 
                                                                             FLA_Obj E, fla_apcaq2ut_t* cntl )
{
  FLA_Obj DTL,   DTR,      D00, D01, D02, 
          DBL,   DBR,      D10, D11, D12,
                           D20, D21, D22;

  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj T1T,
          T2B;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  FLA_Obj ET,              E0,
          EB,              E1,
                           E2;

  FLA_Obj W1TL,  W1TR,
          W1BL,  W1BR;

  dim_t   b_alg, b;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  FLA_Part_2x2( D,    &DTL, &DTR,
                      &DBL, &DBR,     0, 0, FLA_TL );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_TOP );

  FLA_Part_2x1( E,    &ET, 
                      &EB,            0, FLA_TOP );

  while ( FLA_Obj_width( DBR ) > 0 ){

    b = fla_min( b_alg, FLA_Obj_width( DBR ) );

    FLA_Repart_2x2_to_3x3( DTL, /**/ DTR,       &D00, /**/ &D01, &D02,
                        /* ************* */   /* ******************** */
                                                &D10, /**/ &D11, &D12,
                           DBL, /**/ DBR,       &D20, /**/ &D21, &D22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &T2,
                           b, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                        /* ** */            /* ** */
                                              &C1, 
                           CB,                &C2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( ET,                &E0, 
                        /* ** */            /* ** */
                                              &E1, 
                           EB,                &E2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,    &T1T, 
                         &T2B,     b, FLA_TOP );

    FLA_Part_2x2( W1,    &W1TL, &W1TR,
                         &W1BL, &W1BR,     b, FLA_Obj_width( C1 ), FLA_TL );

    // W1TL = inv( triu( T1T ) )' * ( C1 + D1' * E );
    //      = inv( triu( T1T ) )' * ( C1 + D01' * E0 + D11' * E1 );

    FLA_Copy_internal( E1, W1TL,
                       FLA_Cntl_sub_copy( cntl ) );

    FLA_Trmm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, D11, W1TL,
                       FLA_Cntl_sub_trmm1( cntl ) );

    FLA_Gemm_internal( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_ONE, D01, E0, FLA_ONE, W1TL,
                       FLA_Cntl_sub_gemm1( cntl ) );

    FLA_Axpy_internal( FLA_ONE, C1, W1TL,
                       FLA_Cntl_sub_axpy1( cntl ) );

    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T1T, W1TL,
                       FLA_Cntl_sub_trsm( cntl ) );

    // C1 = C1 - W1TL;
    // E  = E  - D1 * W1TL;
	//  => E0 = E0 - D01 * W1TL;
	//     E1 = E1 - D11 * W1TL;

    FLA_Axpy_internal( FLA_MINUS_ONE, W1TL, C1,
                       FLA_Cntl_sub_axpy2( cntl ) );

    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, D01, W1TL, FLA_ONE, E0,
                       FLA_Cntl_sub_gemm2( cntl ) );

    FLA_Trmm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, D11, W1TL,
                       FLA_Cntl_sub_trmm2( cntl ) );

    FLA_Axpy_internal( FLA_MINUS_ONE, W1TL, E1,
                       FLA_Cntl_sub_axpy3( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &DTL, /**/ &DTR,       D00, D01, /**/ D02,
                                                     D10, D11, /**/ D12,
                            /* ************** */  /* ****************** */
                              &DBL, /**/ &DBR,       D20, D21, /**/ D22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ T2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                                                  C1, 
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &ET,                E0, 
                                                  E1, 
                            /* ** */           /* ** */
                              &EB,                E2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

