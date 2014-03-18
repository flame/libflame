
#include "FLAME.h"

FLA_Error FLA_Gemm_cn_blk_var3( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  FLA_Obj BL,    BR,       B0,  B1,  B2;

  FLA_Obj CL,    CR,       C0,  C1,  C2;

  dim_t b;

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  while ( FLA_Obj_width( BL ) < FLA_Obj_width( B ) ){

    b = FLA_Determine_blocksize( BR, FLA_RIGHT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    /* C1 = alpha * A * B1 + C1; */
    FLA_Gemm_internal( FLA_CONJ_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       alpha, A, B1, beta, C1, 
                       FLA_Cntl_sub_gemm( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}
