
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_LU_piv_blk_var3( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl )
{
  FLA_Error r_val = FLA_SUCCESS, r_val_sub = FLA_SUCCESS;
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj AL,    AR,       A0,  A1,  A2;

  FLA_Obj pT,              p0,
          pB,              p1,
                           p2;

  FLA_Obj AB0, AB1;

  dim_t b;


  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_2x1( p,    &pT, 
                      &pB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) &&
          FLA_Obj_width( ATL ) < FLA_Obj_width( A )){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( pT,                &p0, 
                        /* ** */            /* ** */
                                              &p1, 
                           pB,                &p2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // Apply previously computed pivots
    FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, p0, A1,
                               FLA_Cntl_sub_appiv1( cntl ) );

    // A01 = trilu( A00 ) \ A10 
    FLA_Trsm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR, 
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A00, A01,
                       FLA_Cntl_sub_trsm1( cntl ) );

    // A11 = A11 - A10 * A01
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A10, A01, FLA_ONE, A11,
                       FLA_Cntl_sub_gemm1( cntl ) );

    // A21 = A21 - A20 * A01
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A20, A01, FLA_ONE, A21,
                       FLA_Cntl_sub_gemm2( cntl ) );

    // AB1 = / A11 \
    //       \ A21 /
    FLA_Merge_2x1( A11,
                   A21,      &AB1 );

    // AB1, p1 = LU_piv( AB1 )
    r_val = FLA_LU_piv_internal( AB1, p1, 
                                 FLA_Cntl_sub_lu( cntl ) );

    // If the unblocked algorithm returns a null pivot, 
    // update the pivot index and return it.
    if ( r_val == FLA_SUCCESS && r_val_sub >= 0 )
    {
        r_val = FLA_Obj_length( A01 ) + r_val_sub;
    }

    // AB0 = / A10 \
    //       \ A20 /
    FLA_Merge_2x1( A10,
                   A20,      &AB0 );

    // Apply pivots to previous columns
    FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, p1, AB0,
                               FLA_Cntl_sub_appiv2( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &pT,                p0, 
                                                  p1, 
                            /* ** */           /* ** */
                              &pB,                p2,     FLA_TOP );

  }

  if ( FLA_Obj_width( ATR ) > 0 )
  {
    /* Apply pivots to untouched columns */
    FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, p, ATR,
                               FLA_Cntl_sub_appiv1( cntl ) );

    /* ATR = trilu( ATL ) \ ATR */ 
    FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, 
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, ATL, ATR );
  }

  return r_val;
}

#endif
