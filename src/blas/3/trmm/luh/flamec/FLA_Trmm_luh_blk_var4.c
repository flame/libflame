
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Trmm_luh_blk_var4( FLA_Diag diagA, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  FLA_Obj BL,    BR,       B0,  B1,  B2;

  dim_t b;

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_RIGHT );

  while ( FLA_Obj_width( BR ) < FLA_Obj_width( B ) ){

    b = FLA_Determine_blocksize( BL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, &B1, /**/ &B2,
                           b, FLA_LEFT );

    /*------------------------------------------------------------*/

    /* B1 = triu( A )' * B1; */
    FLA_Trmm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diagA,
                       alpha, A, B1,
                       FLA_Cntl_sub_trmm( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, /**/ B1, B2,
                              FLA_RIGHT );

  }

  return FLA_SUCCESS;
}

#endif
