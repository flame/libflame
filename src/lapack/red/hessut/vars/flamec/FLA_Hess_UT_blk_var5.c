
#include "FLAME.h"

FLA_Error FLA_Hess_UT_blk_var5( FLA_Obj A, FLA_Obj T )
{
  FLA_Obj  ATL,   ATR,      A00, A01, A02, 
           ABL,   ABR,      A10, A11, A12,
                            A20, A21, A22;
  FLA_Obj  UT,              U0,
           UB,              U1,
                            U2;
  FLA_Obj  ZT,              Z0,
           ZB,              Z1,
                            Z2;
  FLA_Obj  TL,    TR,       T0, T1, W12; 

  FLA_Obj  U, Z;
  FLA_Obj  UB_l;
  FLA_Obj  ZB_l;
  FLA_Obj  WT_l;
  FLA_Obj  T1_tl;
  FLA_Obj  none, none2, none3;
  FLA_Datatype datatype_A;
  dim_t        m_A;
  dim_t        b_alg, b, bb;

  b_alg      = FLA_Obj_length( T );

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );

  FLA_Obj_create( datatype_A, m_A, b_alg, 0, 0, &U );
  FLA_Obj_create( datatype_A, m_A, b_alg, 0, 0, &Z );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x1( U,    &UT, 
                      &UB,            0, FLA_TOP );
  FLA_Part_2x1( Z,    &ZT, 
                      &ZB,            0, FLA_TOP );
  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT ); 

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) )
  {
    b = min( FLA_Obj_length( ABR ), b_alg );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );
    FLA_Repart_2x1_to_3x1( UT,                &U0, 
                        /* ** */            /* ** */
                                              &U1, 
                           UB,                &U2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( ZT,                &Z0, 
                        /* ** */            /* ** */
                                              &Z1, 
                           ZB,                &Z2,        b, FLA_BOTTOM );
    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &W12,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x2( T1,     &T1_tl, &none,   
                          &none2, &none3,    b, b, FLA_TL ); 

    bb = min( FLA_Obj_length( ABR ) - 1, b_alg );

    FLA_Part_1x2( UB,     &UB_l,  &none,    bb, FLA_LEFT ); 
    FLA_Part_1x2( ZB,     &ZB_l,  &none,    bb, FLA_LEFT ); 

    // [ ABR, UB, ZB, T1 ] = FLA_Hess_UT_step_unb_var5( ABR, UB, ZB, T1, b );
    //FLA_Hess_UT_step_unb_var5( ABR, UB, ZB, T1_tl );
    FLA_Hess_UT_step_opt_var5( ABR, UB, ZB, T1_tl );

    // ATR = ATR - ATR * UB * inv( triu ( T1 ) ) * UB' );
    if ( FLA_Obj_length( ATR ) > 0 )
    {
      // NOTE: We use ZT as temporary workspace.
      FLA_Part_1x2( ZT,     &WT_l,  &none,    bb, FLA_LEFT ); 
      FLA_Part_2x2( T1,     &T1_tl, &none,   
                            &none2, &none3,   bb, bb, FLA_TL ); 

      // WT_l = ATR * UB_l * inv( triu( T1 ) ).
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                         FLA_ONE, ATR, UB_l, FLA_ZERO, WT_l );
      FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, 
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, T1_tl, WT_l );

      // ATR = ATR - WT_l * UB_l'
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                         FLA_MINUS_ONE, WT_l, UB_l, FLA_ONE, ATR );
    }

    //  / A12 \  =  Q11' * / / A12 \ - / Z1 \ * inv( triu( T1 ) ) * U2' \
    //  \ A22 /            \ \ A22 /   \ Z2 /                           /
    //
    //   where Q11 corresponds to the block Householder transformation
    //   associated with UB and T1.
    if ( FLA_Obj_width( A12 ) > 0 )
    {
      FLA_Obj ABR2, ABR2_b;
      FLA_Obj UB_b;

      // NOTE: Since A12.n > 0, we are guaranteed to not be at an edge case,
      // namely the case where bb = b - 1 = ABR.m - 1, thus we are free to use
      // the "full" matrix partitions in this scope block (ie: ZB instead of
      // ZB_l).

      // W12 = U2'
      // W12 = inv( triu( T1 ) ) * W12;
      FLA_Copyt_external( FLA_CONJ_TRANSPOSE, U2, W12 );
      FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                         FLA_NONUNIT_DIAG, FLA_ONE, T1_tl, W12 );

      FLA_Merge_2x1( A12,
                     A22,   &ABR2 );

      //  / A12 \  =  / A12 \  -  / Z1 \ * W12
      //  \ A22 /     \ A22 /     \ Z2 / 
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                         FLA_MINUS_ONE, ZB, W12, FLA_ONE, ABR2 );

      // Omit the top row of UB so it has [implicit] unit diagonal, allowing us
      // to use FLA_Apply_Q_UT() to apply the block Householder transformation
      // corresponding to UB and T1. This trick is valid since the top row of
      // ABR2 would normally be unchanged by the transformation (ie: multiplied
      // by identity).
      FLA_Part_2x1( UB,     &none,   
                            &UB_b,     1, FLA_TOP );
      FLA_Part_2x1( ABR2,   &none,   
                            &ABR2_b,   1, FLA_TOP );

      // Apply Q11' to A12 and A22 from the left:
      //
      //   / A12 \  =  / I - / U1 \ * inv( triu( T1 ) ) * / U1 \' \' / A12 \
      //   \ A22 /     \     \ U2 /                       \ U2 /  /  \ A22 /
      //
      FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                      UB_b, T1_tl, W12, ABR2_b );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &UT,                U0, 
                                                  U1, 
                            /* ** */           /* ** */
                              &UB,                U2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &ZT,                Z0, 
                                                  Z1, 
                            /* ** */           /* ** */
                              &ZB,                Z2,     FLA_TOP );
    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ W12,
                              FLA_LEFT );
  }

  FLA_Obj_free( &U );
  FLA_Obj_free( &Z );

  return FLA_SUCCESS;
}

