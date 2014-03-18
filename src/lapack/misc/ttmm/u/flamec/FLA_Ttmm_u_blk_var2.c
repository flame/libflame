
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Ttmm_u_blk_var2( FLA_Obj A, fla_ttmm_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    // A01 = A01 * triu( A11 )'
    FLA_Trmm_internal( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, A11, A01,
                       FLA_Cntl_sub_trmm( cntl ) );

    // A01 = A01 + A02 * A12'
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_ONE, A02, A12, FLA_ONE, A01,
                       FLA_Cntl_sub_gemm( cntl ) );

    // A11 = triu( A11 ) * triu( A11 )'
    FLA_Ttmm_internal( FLA_UPPER_TRIANGULAR, A11,
                       FLA_Cntl_sub_ttmm( cntl ) );

    // A11 = A11 + A12 * A12'
    FLA_Herk_internal( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_ONE, A12, FLA_ONE, A11,
                       FLA_Cntl_sub_herk( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

  }

  return FLA_SUCCESS;
}

#endif
