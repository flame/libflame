
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Chol_l_blk_var1( FLA_Obj A, fla_chol_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02,
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  dim_t b;

  int r_val = FLA_SUCCESS;

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

    // A10 = A10 * inv( tril( A00 )' )
    FLA_Trsm_internal( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, A00, A10,
                       FLA_Cntl_sub_trsm( cntl ) );

    // A11 = A11 - A10 * A10'
    FLA_Herk_internal( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A10, FLA_ONE, A11,
                       FLA_Cntl_sub_herk( cntl ) );

    // A11 = chol( A11 )
    r_val = FLA_Chol_internal( FLA_LOWER_TRIANGULAR, A11,
                               FLA_Cntl_sub_chol( cntl ) );

    if ( r_val != FLA_SUCCESS )
      return ( FLA_Obj_length( A00 ) + r_val );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
  }

  return r_val;
}

#endif
