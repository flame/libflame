
#include "FLAME.h"

FLA_Error FLA_Trinv_uu_blk_var3( FLA_Obj A, fla_trinv_t* cntl )
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

    // A12 = -triuu( A11 ) \ A12;
    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, 
                       FLA_MINUS_ONE, A11, A12, 
                       FLA_Cntl_sub_trsm1( cntl ) );

    // A02 = A01 * A12 + A02;
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_ONE, A01, A12, FLA_ONE, A02, 
                       FLA_Cntl_sub_gemm( cntl ) );

    // A01 = A01 / triuu( A11 );
    FLA_Trsm_internal( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, 
                       FLA_ONE, A11, A01, 
                       FLA_Cntl_sub_trsm2( cntl ) );

    // A11 = inv( A11 );
    FLA_Trinv_internal( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A11,
                        FLA_Cntl_sub_trinv( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

  }

  return FLA_SUCCESS;
}

