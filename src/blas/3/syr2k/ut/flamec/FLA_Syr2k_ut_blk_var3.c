
#include "FLAME.h"

FLA_Error FLA_Syr2k_ut_blk_var3( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_syr2k_t* cntl )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;

  FLA_Obj BL,    BR,       B0,  B1,  B2;

  FLA_Obj CTL,   CTR,      C00, C01, C02, 
          CBL,   CBR,      C10, C11, C12,
                           C20, C21, C22;

  dim_t b;

  FLA_Scalr_internal( FLA_UPPER_TRIANGULAR, beta, C,
                      FLA_Cntl_sub_scalr( cntl ) );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_TL );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    b = FLA_Determine_blocksize( AR, FLA_RIGHT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b, FLA_RIGHT );

    FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00, /**/ &C01, &C02,
                        /* ************* */   /* ******************** */
                                                &C10, /**/ &C11, &C12,
                           CBL, /**/ CBR,       &C20, /**/ &C21, &C22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    /* C01 = C01 + A0' * B1 */
    FLA_Gemm_internal( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       alpha, A0, B1, FLA_ONE, C01,
                       FLA_Cntl_sub_gemm1( cntl ) );

    /* C12 = C12 + B1' * A2 */
    FLA_Gemm_internal( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       alpha, B1, A2, FLA_ONE, C12,
                       FLA_Cntl_sub_gemm2( cntl ) );

    /* C11 = C11 + A1' * B1 + B1' * A1 */
    FLA_Syr2k_internal( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, 
                        alpha, A1, B1, FLA_ONE, C11,
                        FLA_Cntl_sub_syr2k( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );

    FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00, C01, /**/ C02,
                                                     C10, C11, /**/ C12,
                            /* ************** */  /* ****************** */
                              &CBL, /**/ &CBR,       C20, C21, /**/ C22,
                              FLA_TL );

  }

  return FLA_SUCCESS;
}
