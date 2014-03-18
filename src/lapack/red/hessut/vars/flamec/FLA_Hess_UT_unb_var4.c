
#include "FLAME.h"

FLA_Error FLA_Hess_UT_unb_var4( FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;
  FLA_Obj   Y, Z;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Y );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Z );

  r_val = FLA_Hess_UT_step_unb_var4( A, Y, Z, T );

  FLA_Obj_free( &Y );
  FLA_Obj_free( &Z );

  return r_val;
}

FLA_Error FLA_Hess_UT_step_unb_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T )
{
  FLA_Obj  ATL,   ATR,      A00,  a01,     A02, 
           ABL,   ABR,      a10t, alpha11, a12t,
                            A20,  a21,     A22;
  FLA_Obj  YTL,   YTR,      Y00,  y01,   Y02, 
           YBL,   YBR,      y10t, psi11, y12t,
                            Y20,  y21,   Y22;
  FLA_Obj  ZTL,   ZTR,      Z00,  z01,    Z02, 
           ZBL,   ZBR,      z10t, zeta11, z12t,
                            Z20,  z21,    Z22;
  FLA_Obj  TTL,   TTR,      T00,  t01,   T02, 
           TBL,   TBR,      t10t, tau11, t12t,
                            T20,  t21,   T22;
  FLA_Obj  dT,              d0,
           dB,              delta1,
                            d2;
  FLA_Obj  eT,              e0,
           eB,              epsilon1,
                            e2;
  FLA_Obj  fT,              f0,
           fB,              phi1,
                            f2;
  FLA_Obj  d, e, f;

  FLA_Obj  inv_tau11;
  FLA_Obj  minus_inv_tau11;
  FLA_Obj  first_elem;
  FLA_Obj  last_elem;
  FLA_Obj  beta;
  FLA_Obj  conj_beta;
  FLA_Obj  dot_product;

  FLA_Obj  a10t_l, a10t_r;
  FLA_Obj  a21_t,
           a21_b;
  FLA_Obj  a2;

  FLA_Datatype datatype_A;
  dim_t        m_A;
  dim_t        b_alg;


  b_alg      = FLA_Obj_length( T );

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );

  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &inv_tau11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &minus_inv_tau11 );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &first_elem );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &last_elem );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &beta );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &conj_beta );
  FLA_Obj_create( datatype_A, 1,   1, 0, 0, &dot_product );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &d );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &e );
  FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );

  FLA_Set( FLA_ZERO, Y );
  FLA_Set( FLA_ZERO, Z );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x2( Y,    &YTL, &YTR,
                      &YBL, &YBR,     0, 0, FLA_TL );
  FLA_Part_2x2( Z,    &ZTL, &ZTR,
                      &ZBL, &ZBR,     0, 0, FLA_TL );
  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );
  FLA_Part_2x1( d,    &dT, 
                      &dB,            0, FLA_TOP );
  FLA_Part_2x1( e,    &eT, 
                      &eB,            0, FLA_TOP );
  FLA_Part_2x1( f,    &fT, 
                      &fB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < b_alg )
  {
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( YTL, /**/ YTR,       &Y00,  /**/ &y01,   &Y02,
                        /* ************* */   /* ************************ */
                                                &y10t, /**/ &psi11, &y12t,
                           YBL, /**/ YBR,       &Y20,  /**/ &y21,   &Y22,
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
    FLA_Repart_2x1_to_3x1( dT,                &d0, 
                        /* ** */            /* ****** */
                                              &delta1, 
                           dB,                &d2,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( eT,                &e0, 
                        /* ** */            /* ******** */
                                              &epsilon1, 
                           eB,                &e2,        1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( fT,                &f0, 
                        /* ** */            /* **** */
                                              &phi1, 
                           fB,                &f2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // Save first element of a10_r and set it to one so we can use a10t as
    // u10t in subsequent computations. We will restore a10_r later on.
    if ( FLA_Obj_length( ATL ) > 0 )
    {
      FLA_Part_1x2( a10t,   &a10t_l, &a10t_r,     1, FLA_RIGHT );
      FLA_Copy( a10t_r, last_elem );
      FLA_Set( FLA_ONE, a10t_r );
    }
 
    FLA_Merge_2x1( alpha11,
                   a21,      &a2 );

    // alpha11 = alpha11 - u10t * y10t' - z10t * u10t';
    // a21     = a21     - U20  * y10t' - Z20  * u10t';
    FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a10t, FLA_ONE, a2 );

    // a12t = a12t - u10t * Y20' - z10t * U20';
    FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, A20, z10t, FLA_ONE, a12t );

    // Restore last element of a10t.
    if ( FLA_Obj_length( ATL ) > 0 )
    {
      FLA_Copy( last_elem, a10t_r );
    }

    if ( FLA_Obj_length( A22 ) > 0 )
    {
      FLA_Part_2x1( a21,    &a21_t,
                            &a21_b,            1, FLA_TOP );

      // [ u21, tau11, a21 ] = House( a21 );
      FLA_Househ2_UT( FLA_LEFT,
                      a21_t,
                      a21_b, tau11 );
  
      // inv_tau11            =  1 / tau11;
      // minus_inv_tau11      = -1 / tau11;
      FLA_Set( FLA_ONE, inv_tau11 );
      FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      FLA_Copy( inv_tau11, minus_inv_tau11 );
      FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );

      // Save first element of a21_t and set it to one.
      FLA_Copy( a21_t, first_elem );
      FLA_Set( FLA_ONE, a21_t );
  
      // y21 = A22' * u21;
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, y21 );
  
      // z21 = A22 * u21;
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, z21 );
 
      // y21 = y21 - Y20 * ( U20' * u21 ) - U20 * ( Z20' * u21 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, d0 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Y20, a21, FLA_ZERO, e0 );
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, a21, FLA_ZERO, f0 );
  
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, y21 );
  
      // t01 = U20' * u21;
      FLA_Copy( d0, t01 );

      // z21 = z21 - U20 * ( Y20' * u21 ) - Z20 * ( U20' * u21 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, e0, FLA_ONE, z21 );
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, d0, FLA_ONE, z21 );
  
      // beta      = u21' * z21 / 2;
      // conj_beta = conj(beta);
      FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      FLA_Inv_scal( FLA_TWO, beta );
      FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
  
      // y21' = ( y21' - beta / tau * u21' ) / tau;
      // y21  = ( y21 - conj(beta) / tau * u21 ) / tau;
      FLA_Scal( minus_inv_tau11, conj_beta );
      FLA_Axpy( conj_beta, a21, y21 );
      FLA_Scal( inv_tau11, y21 );
  
      // z21 = ( z21 - beta / tau * u21 ) / tau;
      FLA_Scal( minus_inv_tau11, beta );
      FLA_Axpy( beta, a21, z21 );
      FLA_Scal( inv_tau11, z21 );
  
      // a12t = a12t * ( I - u21 * u21' / tau );
      //      = a12t - ( a12t * u21 ) * u21' / tau;
      FLA_Dot( a12t, a21, dot_product );
      FLA_Scal( minus_inv_tau11, dot_product );
      FLA_Axpyt( FLA_CONJ_TRANSPOSE, dot_product, a21, a12t );
  
      // A02 = A02 * ( I - u21 * u21' / tau );
      //     = A02 - ( A02 * u21 ) * u21' / tau;
      FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, e0 );
      FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, e0, a21, A02 );
  
      // Restore first element of a21.
      FLA_Copy( first_elem, a21_t );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &YTL, /**/ &YTR,       Y00,  y01,   /**/ Y02,
                                                     y10t, psi11, /**/ y12t,
                            /* ************** */  /* ********************** */
                              &YBL, /**/ &YBR,       Y20,  y21,   /**/ Y22,
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
    FLA_Cont_with_3x1_to_2x1( &dT,                d0, 
                                                  delta1, 
                            /* ** */           /* ****** */
                              &dB,                d2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &eT,                e0, 
                                                  epsilon1, 
                            /* ** */           /* ******** */
                              &eB,                e2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &fT,                f0, 
                                                  phi1, 
                            /* ** */           /* **** */
                              &fB,                f2,     FLA_TOP );
  }

  FLA_Obj_free( &inv_tau11 );
  FLA_Obj_free( &minus_inv_tau11 );
  FLA_Obj_free( &first_elem );
  FLA_Obj_free( &last_elem );
  FLA_Obj_free( &beta );
  FLA_Obj_free( &conj_beta );
  FLA_Obj_free( &dot_product );
  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
  FLA_Obj_free( &f );

  return FLA_SUCCESS;
}

