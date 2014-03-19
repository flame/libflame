/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_QR_UT_piv_unb_var2( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj TTL,   TTR,      T00,  t01,   T02, 
          TBL,   TBR,      t10t, tau11, t12t,
                           T20,  t21,   T22;

  FLA_Obj pT,              p0,
          pB,              pi1,
                           p2;

  FLA_Obj wT,              w0,
          wB,              omega1,
                           w2;

  FLA_Obj ab1, v;

  // Create workspace
  FLA_Obj_create( FLA_Obj_datatype( T ), 1, FLA_Obj_width( T ), 0, 0, &v );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  FLA_Part_2x1( p,    &pT,
                      &pB,            0, FLA_TOP );

  FLA_Part_2x1( w,    &wT,
                      &wB,            0, FLA_TOP );

  while ( FLA_Obj_min_dim( pB ) > 0 ) {

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************ */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( pT,                &p0,
                        /* ** */            /* *** */
                                              &pi1,
                           pB,                &p2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( wT,                &w0,
                        /* ** */            /* *** */
                                              &omega1,
                           wB,                &w2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/


    //  ** Ignore minus inputs for LAPACK compatability.
    if ( FLA_Obj_lt( pi1, FLA_ZERO ) == FALSE )
    {
      // ** Determine pivot index
      FLA_Amax_external( wB, pi1 );

      // ** BLIS returns -1 if it fails to search the maximum value
      if ( FLA_Obj_lt( pi1, FLA_ZERO ) == TRUE )
        FLA_Set( FLA_ZERO, pi1 );

      // ** Apply a pivot on column norms
      FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, wB );

      // ** Apply a pivot on ABR
      FLA_Apply_pivots( FLA_RIGHT, FLA_TRANSPOSE, pi1, ABR );

      // ** Apply a pivot on TTR
      FLA_Apply_pivots( FLA_RIGHT, FLA_TRANSPOSE, pi1, TTR );
    }
    else
    {
      // ** Do not pivot.
      FLA_Set( FLA_ZERO, pi1 );
    }

    // ** Update the pivot column
    FLA_Merge_2x1( alpha11,
                   a21, &ab1 );

    // ab1 = ab1 - ABL t01
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, ABL, t01, FLA_ONE, ab1 );

    // ** Find the householder reflector on that column
    FLA_Househ2_UT( FLA_LEFT, alpha11,
                              a21,     tau11 );

    // ** Update the pivot row
    FLA_Apply_H2_UT_piv_row( tau11, a12t, a10t, T02,
                             a21,   A22,  A20,  t12t,
                             v );

    // ** Apply pivots on ATR 
    FLA_Apply_pivots( FLA_RIGHT, FLA_TRANSPOSE, pi1, ATR );

    // ** Norm downdate w2 = w2 - columnwisenorm2(a12t)
    FLA_QR_UT_piv_colnorm( FLA_MINUS_ONE, a12t, w2 );

    // ** Update T matrix
    // t01 = a10t' + A20' * u21; 
    FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ********************** */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &pT,                p0,
                                                  pi1,
                            /* ** */           /* *** */
                              &pB,                p2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &wT,                w0,
                                                  omega1,
                            /* ** */           /* *** */
                              &wB,                w2,     FLA_TOP );
  }

  // Free the workspace
  FLA_Obj_free( &v);

  return FLA_SUCCESS;
}

