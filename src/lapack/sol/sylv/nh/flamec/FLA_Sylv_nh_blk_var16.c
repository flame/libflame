
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Sylv_nh_blk_var16( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_BR );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( ABR ) < FLA_Obj_length( A ) ){

    b = FLA_Determine_blocksize( CT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, &A01, /**/ &A02,
                                                &A10, &A11, /**/ &A12,
                        /* ************* */   /* ******************** */
                           ABL, /**/ ABR,       &A20, &A21, /**/ &A22,
                           b, b, FLA_TL );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                                              &C1, 
                        /* ** */            /* ** */
                           CB,                &C2,        b, FLA_TOP );

    // Loop Invariant:
    // CT = CT - ATR * sylv( ABR, B', CB )
    // CB = sylv( ABR, B', CB )

    /*------------------------------------------------------------*/

    // C1 = sylv( A11, B', C1 );
    FLA_Sylv_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, 
                       isgn, A11, B, C1, scale,
                       FLA_Cntl_sub_sylv1( cntl ) );

    // C0 = C0 - A01 * C1;
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A01, C1, FLA_ONE, C0,
                       FLA_Cntl_sub_gemm1( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                            /* ************** */  /* ****************** */
                                                     A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                              FLA_BR );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                            /* ** */           /* ** */
                                                  C1, 
                              &CB,                C2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
