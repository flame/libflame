/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_QUD_UT_lhfc_blk_var1( FLA_Obj T, FLA_Obj W,
                                                     FLA_Obj R,
                                          FLA_Obj U, FLA_Obj C,
                                          FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl )
{
  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj UL,    UR,       U0,  U1,  U2;

  FLA_Obj VL,    VR,       V0,  V1,  V2;

  FLA_Obj RT,              R0,
          RB,              R1,
                           R2;

  FLA_Obj T1T,
          T1B;

  FLA_Obj W1TL,  W1TR,
          W1BL,  W1BR;

  dim_t   b_alg, b;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  FLA_Part_1x2( U,    &UL,  &UR,      0, FLA_LEFT );

  FLA_Part_1x2( V,    &VL,  &VR,      0, FLA_LEFT );

  FLA_Part_2x1( R,    &RT, 
                      &RB,            0, FLA_TOP );

  while ( FLA_Obj_width( UL ) < FLA_Obj_width( U ) ){

    b = fla_min( b_alg, FLA_Obj_width( UR ) );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &T2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( UL,  /**/ UR,        &U0, /**/ &U1, &U2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( VL,  /**/ VR,        &V0, /**/ &V1, &V2,
                           b, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( RT,                &R0, 
                        /* ** */            /* ** */
                                              &R1, 
                           RB,                &R2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,    &T1T,
                         &T1B,    b, FLA_TOP );

    FLA_Part_2x2( W,     &W1TL, &W1TR,
                         &W1BL, &W1BR,     b, FLA_Obj_width( R1 ), FLA_TL );

    // W1TL = R1;

    FLA_Copyt_internal( FLA_NO_TRANSPOSE, R1, W1TL,
                        FLA_Cntl_sub_copyt( cntl ) );

    // W1TL = inv( triu( T1T ) )' * ( R1 + U1' * C + V1' * D );

    FLA_Gemm_internal( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_ONE, U1, C, FLA_ONE, W1TL,
                       FLA_Cntl_sub_gemm1( cntl ) );

    FLA_Gemm_internal( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_ONE, V1, D, FLA_ONE, W1TL,
                       FLA_Cntl_sub_gemm2( cntl ) );

    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T1T, W1TL,
                       FLA_Cntl_sub_trsm( cntl ) );

    // R1 = R1 - W1TL;
    // C  = C  - U1 * W1TL;
    // D  = D  + V1 * W1TL;

    FLA_Axpyt_internal( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, W1TL, R1,
                        FLA_Cntl_sub_axpyt( cntl ) );

    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, U1, W1TL, FLA_ONE, C,
                       FLA_Cntl_sub_gemm3( cntl ) );

    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_ONE, V1, W1TL, FLA_ONE, D,
                       FLA_Cntl_sub_gemm4( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ T2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &UL,  /**/ &UR,        U0, U1, /**/ U2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &VL,  /**/ &VR,        V0, V1, /**/ V2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &RT,                R0, 
                                                  R1, 
                            /* ** */           /* ** */
                              &RB,                R2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

