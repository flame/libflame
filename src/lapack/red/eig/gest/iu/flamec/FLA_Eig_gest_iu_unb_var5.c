
#include "FLAME.h"

FLA_Error FLA_Eig_gest_iu_unb_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj BTL,   BTR,      B00,  b01,    B02, 
          BBL,   BBR,      b10t, beta11, b12t,
                           B20,  b21,    B22;

  //FLA_Obj yL,    yR,       y10t, psi11,  y12t;

  //FLA_Obj y12t_t,
  //        y12t_b;

  FLA_Obj psi11, y12t,
          y21,   Y22;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  //FLA_Part_1x2( Y,    &yL,  &yR,      0, FLA_LEFT );

  FLA_Part_2x2( Y,    &psi11, &y12t,
                      &y21,   &Y22,     1, 1, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x2_to_3x3( BTL, /**/ BTR,       &B00,  /**/ &b01,    &B02,
                        /* ************* */   /* ************************* */
                                                &b10t, /**/ &beta11, &b12t,
                           BBL, /**/ BBR,       &B20,  /**/ &b21,    &B22,
                           1, 1, FLA_BR );

    //FLA_Repart_1x2_to_1x3( yL,  /**/ yR,        &y10t, /**/ &psi11,  &y12t,
    //                       1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    //FLA_Part_2x1( y12t,   &y12t_t,
    //                      &y12t_b,    1, FLA_TOP );

    // alpha11 = inv(beta11) * alpha11 * inv(conj(beta11));
    //         = inv(beta11) * alpha11 * inv(beta11);
    FLA_Inv_scal_external( beta11, alpha11 );
    FLA_Inv_scal_external( beta11, alpha11 );

    //// y12t = alpha11 * b12t;
    //FLA_Copy_external( b12t, y12t_t );
    //FLA_Scal_external( alpha11, y12t_t );
    // psi11 = - 1/2 * alpha11;
    FLA_Copy_external( alpha11, psi11 );
    FLA_Scal_external( FLA_MINUS_ONE_HALF, psi11 );

    // a12t = inv(conj(beta11)) * a12t;
    //      = inv(beta11) * a12t;
    FLA_Inv_scal_external( beta11, a12t );

    //// a12t = a12t - 1/2 * y12t;
    //FLA_Axpy_external( FLA_MINUS_ONE_HALF, y12t_t, a12t );
    // a12t = a12t - 1/2 * alpha11 * b12t;
    FLA_Axpy_external( psi11, b12t, a12t );

    // A22 = A22 - a12t' * b12t - b12t' * a12t;
    FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE,
                        FLA_MINUS_ONE, a12t, b12t, A22 );

    //// a12t = a12t - 1/2 * y12t;
    //FLA_Axpy_external( FLA_MINUS_ONE_HALF, y12t_t, a12t );
    // a12t = a12t - 1/2 * alpha11 * b12t;
    FLA_Axpy_external( psi11, b12t, a12t );

    // a12t   = a12t * inv( tril( B22 ) );
    // a12t^T = inv( triu( B22 )^T ) * a12t^T;
    FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
                       B22, a12t );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &BTL, /**/ &BTR,       B00,  b01,    /**/ B02,
                                                     b10t, beta11, /**/ b12t,
                            /* ************** */  /* *********************** */
                              &BBL, /**/ &BBR,       B20,  b21,    /**/ B22,
                              FLA_TL );

    //FLA_Cont_with_1x3_to_1x2( &yL,  /**/ &yR,        y10t, psi11, /**/  y12t,
    //                          FLA_LEFT );
  }

  return FLA_SUCCESS;
}

