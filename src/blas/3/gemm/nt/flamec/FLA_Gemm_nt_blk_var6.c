
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Gemm_nt_blk_var6( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;

  FLA_Obj BL,    BR,       B0,  B1,  B2;

  dim_t b;

  FLA_Scal_internal( beta, C,
                     FLA_Cntl_sub_scal( cntl ) );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_RIGHT );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_RIGHT );

  while ( FLA_Obj_width( AR ) < FLA_Obj_width( A ) ){

    b = FLA_Determine_blocksize( AL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, &A1, /**/ &A2,
                           b, FLA_LEFT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, &B1, /**/ &B2,
                           b, FLA_LEFT );

    /*------------------------------------------------------------*/

    /* C = alpha * A1 * B1' + C; */
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
                       alpha, A1, B1, FLA_ONE, C,
                       FLA_Cntl_sub_gemm( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, /**/ A1, A2,
                              FLA_RIGHT );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, /**/ B1, B2,
                              FLA_RIGHT );

  }

  return FLA_SUCCESS;
}

#endif
