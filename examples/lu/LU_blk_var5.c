#include "FLAME.h"
#include "LU_prototypes.h"

FLA_Error LU_blk_var5( FLA_Obj A, int nb_alg )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    b = min( FLA_Obj_length( ABR ), nb_alg );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    /* A11 = LU( A11 ); */
    LU_unb_var5( A11 );

    //Alternative 1: Call FLA_LU_nopiv instead:
    // FLA_LU_nopiv( A11 );

    /* A12 = inv( trilu( A11 ) ) * A12; */
    FLA_Trsm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
              FLA_ONE, A11, A12 );

    /* A21 = A21 * inv( triu( A11 ) ); */
    FLA_Trsm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, A11, A21 );

    /* A22 = A22 - A21 * A12; */
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
              FLA_MINUS_ONE, A21, A12, FLA_ONE, A22 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
  }

  return FLA_SUCCESS;
}



