
#include "FLAME.h"

FLA_Error FLA_Trsv_lt_blk_var2( FLA_Diag diagA, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj xT,              x0,
          xB,              x1,
                           x2;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_BR );

  FLA_Part_2x1( x,    &xT, 
                      &xB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( ABR ) < FLA_Obj_length( A ) ){

    b = FLA_Determine_blocksize( ATL, FLA_TL, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, &A01, /**/ &A02,
                                                &A10, &A11, /**/ &A12,
                        /* ************* */   /* ******************** */
                           ABL, /**/ ABR,       &A20, &A21, /**/ &A22,
                           b, b, FLA_TL );

    FLA_Repart_2x1_to_3x1( xT,                &x0, 
                                              &x1, 
                        /* ** */            /* ** */
                           xB,                &x2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    /* x1 = tril( A11' ) \ x1; */
    FLA_Trsv_internal( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, diagA,
                       A11, x1,
                       FLA_Cntl_sub_trsv( cntl ) );

    /* x0 = x0 - A10' * x1; */
    FLA_Gemv_internal( FLA_TRANSPOSE,
                       FLA_MINUS_ONE, A10, x1, FLA_ONE, x0,
                       FLA_Cntl_sub_gemv( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                            /* ************** */  /* ****************** */
                                                     A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                              FLA_BR );

    FLA_Cont_with_3x1_to_2x1( &xT,                x0, 
                            /* ** */           /* ** */
                                                  x1, 
                              &xB,                x2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

