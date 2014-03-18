
#include "FLAME.h"

FLA_Error FLA_Hess_UT_blk_var1( FLA_Obj A, FLA_Obj T )
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
  FLA_Obj  TL,    TR,       T0, T1, T2; 

  FLA_Obj  U, Z;
  FLA_Obj  ABR_l;
  FLA_Obj  UB_l;
  FLA_Obj  T1_tl;
  FLA_Obj  WT_l;
  FLA_Obj  none, none2, none3;
  FLA_Obj  UB_tl,
           UB_bl;
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
    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &T2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x2( T1,     &T1_tl, &none,   
                          &none2, &none3,   b, b, FLA_TL ); 

    bb = min( FLA_Obj_length( ABR ) - 1, b_alg );

    FLA_Part_1x2( ABR,    &ABR_l, &none,    bb, FLA_LEFT ); 
    FLA_Part_1x2( UB,     &UB_l,  &none,    bb, FLA_LEFT ); 

    // [ ABR, T1 ] = FLA_Hess_UT_step_unb_var1( ABR, T1, b );
    //FLA_Hess_UT_step_unb_var1( ABR, T1_tl );
    FLA_Hess_UT_step_opt_var1( ABR, T1_tl );

    // Build UB from ABR, with explicit unit subdiagonal and zeros.
    FLA_Copy_external( ABR_l, UB_l );
    FLA_Part_2x1( UB_l,   &UB_tl, 
                          &UB_bl,            1, FLA_TOP );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_UNIT_DIAG, UB_bl );
    FLA_Set( FLA_ZERO, UB_tl );

    // ATR = ATR - ATR * UB * inv( triu( T1 ) ) * UB' );
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
    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ T2,
                              FLA_LEFT );
  }

  FLA_Obj_free( &U );
  FLA_Obj_free( &Z );

  return FLA_SUCCESS;
}

