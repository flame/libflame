
#include "FLAME.h"

FLA_Error FLA_CAQR_UT_inc_blk_var1( FLA_Obj A, FLA_Obj TW, fla_caqrutinc_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TTL,   WTR,      T00, W01, W02, 
          TBL,   TBR,      T10, T11, W12,
                           T20, T21, T22;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( TW,   &TTL, &WTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x2_to_3x3( TTL, /**/ WTR,       &T00, /**/ &W01, &W02,
                        /* ************* */   /* ******************** */
                                                &T10, /**/ &T11, &W12,
                           TBL, /**/ TBR,       &T20, /**/ &T21, &T22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    FLA_CAQR2_UT_internal( A11,
                           A21, T21, 
                           FLA_Cntl_sub_caqr2ut( cntl ) );


    if ( FLA_Obj_width( A12 ) > 0 )
    {
      FLA_Apply_CAQ2_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                  A21, T21, W12, A12,
                                                 A22,
                                  FLA_Cntl_sub_apcaq2ut( cntl ) );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &WTR,       T00, W01, /**/ W02,
                                                     T10, T11, /**/ W12,
                            /* ************** */  /* ****************** */
                              &TBL, /**/ &TBR,       T20, T21, /**/ T22,
                              FLA_TL );

  }

  return FLA_SUCCESS;
}

