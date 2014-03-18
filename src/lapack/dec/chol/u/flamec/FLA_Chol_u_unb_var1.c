
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Chol_u_unb_var1( FLA_Obj A )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02,
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  int r_val = FLA_SUCCESS;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/

    // a01 = inv( tril( A00 )' ) * a01
    FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a01 );

    // alpha11 = alpha11 - a01' * a01
    FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, a01, FLA_ONE, alpha11 );

    // alpha11 = sqrt( alpha11 )
    r_val = FLA_Sqrt( alpha11 );

    if ( r_val != FLA_SUCCESS )
      return ( FLA_Obj_length( A00 ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );
  }

  return r_val;
}

#endif
