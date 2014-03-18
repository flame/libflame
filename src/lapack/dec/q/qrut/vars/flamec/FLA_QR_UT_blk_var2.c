
#include "FLAME.h"

FLA_Error FLA_QR_UT_blk_var2( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TTL,   TTR,      T00, T01, T02, 
          TBL,   TBR,      T10, T11, W12,
                           T20, T21, T22;

  FLA_Obj AB1,   AB2;

  dim_t   b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00, /**/ &T01, &T02,
                        /* ************* */   /* ******************** */
                                                &T10, /**/ &T11, &W12,
                           TBL, /**/ TBR,       &T20, /**/ &T21, &T22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    FLA_Merge_2x1( A11,
                   A21,   &AB1 );

    // Perform a QR factorization via the UT transform on AB1:
    //
    //   / A11 \ -> QB1 R11
    //   \ A21 /
    //
    // where:
    //  - QB1 is formed from UB1 (which is stored column-wise below the
    //    diagonal of AB1) and T11 (which is stored to the upper triangle
    //    of T11).
    //  - R11 is stored to the upper triangle of AB1.
  
    FLA_QR_UT_internal( AB1, T11, 
                        FLA_Cntl_sub_qrut( cntl ) );


    if ( FLA_Obj_width( A12 ) > 0 )
    {
      FLA_Merge_2x1( A12,
                     A22,   &AB2 );

      // Apply the Householder transforms associated with UB1 and T11 to 
      // AB2:
      //
      //   / A12 \ := QB1' / A12 \
      //   \ A22 /         \ A22 /
      //
      // where QB1 is formed from UB1 and T11.

      FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                               AB1, T11, W12, AB2,
                               FLA_Cntl_sub_apqut( cntl ) );
    }


    // Update T
    //
    // T01 = A10' * U11 + A20' * U21;
    //
    // Recall: U11 = trilu( A11 );
    //         U21 = A21;

    FLA_Copyt_external( FLA_CONJ_TRANSPOSE, A10, T01 );
    FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A11, T01 );
    FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_ONE, A20, A21, FLA_ONE, T01 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00, T01, /**/ T02,
                                                     T10, T11, /**/ W12,
                            /* ************** */  /* ****************** */
                              &TBL, /**/ &TBR,       T20, T21, /**/ T22,
                              FLA_TL );
  }

  return FLA_SUCCESS;
}

