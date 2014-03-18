
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Gemm_hh_blk_var4( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  FLA_Obj CL,    CR,       C0,  C1,  C2;

  dim_t b;

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_BOTTOM );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_RIGHT );

  while ( FLA_Obj_length( BB ) < FLA_Obj_length( B ) ){

    b = FLA_Determine_blocksize( BT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                                              &B1, 
                        /* ** */            /* ** */
                           BB,                &B2,        b, FLA_TOP );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, &C1, /**/ &C2,
                           b, FLA_LEFT );

    /*------------------------------------------------------------*/

    /* C1 = alpha * A' * B1' + C1; */
    FLA_Gemm_internal( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       alpha, A, B1, beta, C1,
                       FLA_Cntl_sub_gemm( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                            /* ** */           /* ** */
                                                  B1, 
                              &BB,                B2,     FLA_BOTTOM );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, /**/ C1, C2,
                              FLA_RIGHT );

  }

  return FLA_SUCCESS;
}

#endif
