#include "FLAME.h"
#include "Chol_prototypes.h"

FLA_Error Chol_unb_var2( FLA_Obj A )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* *************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/

    /* alpha11 = alpha11 - a10t * a10; */
    FLA_Dots( FLA_MINUS_ONE, a10t, a10t, FLA_ONE, alpha11 );

    /* alpha11 = sqrt( alpha11 ); */
    FLA_Sqrt( alpha11 );

    /* a21 = a21 - A20 * a10; */
    FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, a10t, FLA_ONE, a21 );

    /* a21 = a21 / alpha11; */
    FLA_Inv_scal( alpha11, a21 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************* */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );
  }

  return FLA_SUCCESS;
}
