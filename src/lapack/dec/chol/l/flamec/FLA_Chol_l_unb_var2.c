/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Chol_l_unb_var2( FLA_Obj A )
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

    // alpha11 = alpha11 - a10t * a10t'
    FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a10t, a10t, FLA_ONE, alpha11 );

    // a21 = a21 - A20 * a10t'
    FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A20, a10t, FLA_ONE, a21 );

    // alpha11 = sqrt( alpha11 )
    r_val = FLA_Sqrt( alpha11 );

    if ( r_val != FLA_SUCCESS )
      return ( FLA_Obj_length( A00 ) );

    // a21 = a21 / alpha11
    FLA_Inv_scal_external( alpha11, a21 );

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
