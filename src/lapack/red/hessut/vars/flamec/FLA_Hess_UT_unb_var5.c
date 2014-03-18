
#include "FLAME.h"

FLA_Error FLA_Hess_UT_unb_var5( FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;
  FLA_Obj   U, Z;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &U );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Z );

  r_val = FLA_Hess_UT_step_unb_var5( A, U, Z, T );

  FLA_Obj_free( &U );
  FLA_Obj_free( &Z );

  return r_val;
}

FLA_Error FLA_Hess_UT_step_unb_var5( FLA_Obj A, FLA_Obj U, FLA_Obj Z, FLA_Obj T )
{
  FLA_Obj  ATL,   ATR,      A00,  a01,     A02, 
           ABL,   ABR,      a10t, alpha11, a12t,
                            A20,  a21,     A22;
  FLA_Obj  UTL,   UTR,      U00,  u01,       U02, 
           UBL,   UBR,      u10t, upsilon11, u12t,
                            U20,  u21,       U22;
  FLA_Obj  ZTL,   ZTR,      Z00,  z01,    Z02, 
           ZBL,   ZBR,      z10t, zeta11, z12t,
                            Z20,  z21,    Z22;
  FLA_Obj  TTL,   TTR,      T00,  t01,   T02, 
           TBL,   TBR,      t10t, tau11, t12t,
                            T20,  t21,   T22;
  FLA_Obj  wT,              w0,
           wB,              omega1,
                            w2;
  FLA_Obj  w;

  FLA_Obj  a21_t,
           a21_b;
  FLA_Obj  u21_t,
           u21_b;

  FLA_Datatype datatype_A;
  dim_t        m_A;
  dim_t        b_alg;


  b_alg      = FLA_Obj_length( T );

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );

  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );

  FLA_Set( FLA_ZERO, U );
  FLA_Set( FLA_ZERO, Z );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x2( U,    &UTL, &UTR,
                      &UBL, &UBR,     0, 0, FLA_TL );
  FLA_Part_2x2( Z,    &ZTL, &ZTR,
                      &ZBL, &ZBR,     0, 0, FLA_TL );
  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );
  FLA_Part_2x1( w,    &wT, 
                      &wB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < b_alg )
  {
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00,  /**/ &u01,       &U02,
                        /* ************* */   /* **************************** */
                                                &u10t, /**/ &upsilon11, &u12t,
                           UBL, /**/ UBR,       &U20,  /**/ &u21,       &U22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( ZTL, /**/ ZTR,       &Z00,  /**/ &z01,    &Z02,
                        /* ************* */   /* ************************* */
                                                &z10t, /**/ &zeta11, &z12t,
                           ZBL, /**/ ZBR,       &Z20,  /**/ &z21,    &Z22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************** */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );
    FLA_Repart_2x1_to_3x1( wT,                &w0, 
                        /* ** */            /* ****** */
                                              &omega1, 
                           wB,                &w2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    if ( FLA_Obj_length( ATL ) > 0 )
    {
      // w0 = inv( triu( T00 ) ) * u10t';
      FLA_Copyt( FLA_CONJ_TRANSPOSE, u10t, w0 );
      FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                T00, w0 );

      // a01     = a01     - Z00  * w0;
      // alpha11 = alpha11 - z10t * w0;
      // a21     = a21     - Z20  * w0;
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z00, w0, FLA_ONE, a01 );
      FLA_Dots( FLA_MINUS_ONE, z10t, w0, FLA_ONE, alpha11 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, w0, FLA_ONE, a21 );
 
      // w0 = inv( triu( T00 ) )' * ( U00' * a01 + u10t' * alpha11 + U20' * a21 );
      FLA_Trmvsx( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                  FLA_ONE, U00, a01, FLA_ZERO, w0 );
      FLA_Axpyt( FLA_CONJ_TRANSPOSE, alpha11, u10t, w0 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, U20, a21, FLA_ONE, w0 );
      
      FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                T00, w0 );
  
      // a01     = a01     - U00  * w0;
      // alpha11 = alpha11 - u10t * w0;
      // a21     = a21     - U20  * w0;
      FLA_Trmvsx( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                  FLA_MINUS_ONE, U00, w0, FLA_ONE, a01 );
      FLA_Dots( FLA_MINUS_ONE, u10t, w0, FLA_ONE, alpha11 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, U20, w0, FLA_ONE, a21 );
    }
 
    if ( FLA_Obj_length( a21 ) > 0 )
    {
      FLA_Part_2x1( a21,    &a21_t,
                            &a21_b,            1, FLA_TOP );

      // [ u21, tau11, a21 ] = House( a21 );
      FLA_Househ2_UT( FLA_LEFT,
                      a21_t,
                      a21_b, tau11 );

      // u21 := a21;
      FLA_Copy( a21, u21 );

      // Explicitly set the first element of the Householder vector so we
      // can use it in regular computations.
      FLA_Part_2x1( u21,    &u21_t,
                            &u21_b,            1, FLA_TOP );
      FLA_Set( FLA_ONE, u21_t );
  
      // z01    = A02  * u21;
      // zeta11 = a12t * u21;
      // z21    = A22  * u21;
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, u21, FLA_ZERO, z01 );
      FLA_Dot( a12t, u21, zeta11 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, u21, FLA_ZERO, z21 );
  
      // t01 = U20' * u21;
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, U20, u21, FLA_ZERO, t01 );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00,  u01,       /**/ U02,
                                                     u10t, upsilon11, /**/ u12t,
                            /* ************** */  /* ************************** */
                              &UBL, /**/ &UBR,       U20,  u21,       /**/ U22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &ZTL, /**/ &ZTR,       Z00,  z01,    /**/ Z02,
                                                     z10t, zeta11, /**/ z12t,
                            /* ************** */  /* *********************** */
                              &ZBL, /**/ &ZBR,       Z20,  z21,    /**/ Z22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ************************ */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &wT,                w0, 
                                                  omega1, 
                            /* ** */           /* ****** */
                              &wB,                w2,     FLA_TOP );
  }

  FLA_Obj_free( &w );

  return FLA_SUCCESS;
}

