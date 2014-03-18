
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Chol_u_unb_var2( FLA_Obj A )
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

    // alpha11 = alpha11 - a01' * a01
    FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, a01, FLA_ONE, alpha11 );

    // a12t = a12t - a01' * A02
    // a12t' = a12t' - A02' * conj(a01)
    FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, a01, FLA_ONE, a12t );

    // alpha11 = sqrt( alpha11 )
    r_val = FLA_Sqrt( alpha11 );

    if ( r_val != FLA_SUCCESS )
      return ( FLA_Obj_length( A00 ) );

    // a12t = a12t / alpha11
    FLA_Inv_scal_external( alpha11, a12t );

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
