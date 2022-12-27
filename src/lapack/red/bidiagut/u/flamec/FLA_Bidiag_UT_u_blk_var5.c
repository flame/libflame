/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u_blk_var5( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  FLA_Obj  ATL,   ATR,      A00, A01, A02, 
           ABL,   ABR,      A10, A11, A12,
                            A20, A21, A22;
  FLA_Obj  UT,              U0,
           UB,              U1,
                            U2;
  FLA_Obj  VT,              V0,
           VB,              V1,
                            V2;
  FLA_Obj  YT,              Y0,
           YB,              Y1,
                            Y2; 
  FLA_Obj  ZT,              Z0,
           ZB,              Z1,
                            Z2;
  FLA_Obj  TUL,   TUR,      TU0, TU1, TU2;
  FLA_Obj  TVL,   TVR,      TV0, TV1, TV2;

  FLA_Obj  U, V, Y, Z;
  FLA_Obj  ABR_l, ABR_t;
  FLA_Obj  UB_l, U2_l;
  FLA_Obj  VB_l, V2_l;
  FLA_Obj  YB_l, Y2_l;
  FLA_Obj  ZB_l, Z2_l;
  FLA_Obj  TU1_tl;
  FLA_Obj  TV1_tl;
  FLA_Obj  none, none2, none3;
  FLA_Obj  VB_tl,
           VB_bl;
  FLA_Datatype datatype_A;
  dim_t        m_A, n_A;
  dim_t        b_alg, b;

  b_alg      = FLA_Obj_length( TU );

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );
  n_A        = FLA_Obj_width( A );

  FLA_Obj_create( datatype_A, m_A, b_alg, 0, 0, &U );
  FLA_Obj_create( datatype_A, n_A, b_alg, 0, 0, &V );
  FLA_Obj_create( datatype_A, n_A, b_alg, 0, 0, &Y );
  FLA_Obj_create( datatype_A, m_A, b_alg, 0, 0, &Z );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,   0, 0, FLA_TL );
  FLA_Part_2x1( U,    &UT,
                      &UB,            0, FLA_TOP );
  FLA_Part_2x1( V,    &VT,
                      &VB,            0, FLA_TOP );
  FLA_Part_2x1( Y,    &YT,
                      &YB,            0, FLA_TOP );
  FLA_Part_2x1( Z,    &ZT,
                      &ZB,            0, FLA_TOP );
  FLA_Part_1x2( TU,   &TUL, &TUR,      0, FLA_LEFT ); 
  FLA_Part_1x2( TV,   &TVL, &TVR,      0, FLA_LEFT ); 

  while ( FLA_Obj_min_dim( ABR ) > 0 )
  {
    b = fla_min( FLA_Obj_min_dim( ABR ), b_alg );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );
    FLA_Repart_2x1_to_3x1( UT,                &U0,
                        /* ** */            /* ** */
                                              &U1,
                           UB,                &U2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( VT,                &V0,
                        /* ** */            /* ** */
                                              &V1,
                           VB,                &V2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( YT,                &Y0,
                        /* ** */            /* ** */
                                              &Y1,
                           YB,                &Y2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( ZT,                &Z0,
                        /* ** */            /* ** */
                                              &Z1,
                           ZB,                &Z2,        b, FLA_BOTTOM );
    FLA_Repart_1x2_to_1x3( TUL, /**/ TUR,       &TU0, /**/ &TU1, &TU2,
                           b, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( TVL, /**/ TVR,       &TV0, /**/ &TV1, &TV2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x2( TU1,     &TU1_tl, &none,   
                           &none2,  &none3,   b, b, FLA_TL ); 

    FLA_Part_2x2( TV1,     &TV1_tl, &none,   
                           &none2,  &none3,   b, b, FLA_TL ); 

    FLA_Part_1x2( ABR,    &ABR_l, &none,    b, FLA_LEFT );
    FLA_Part_2x1( ABR,    &ABR_t,
                          &none,            b, FLA_TOP );

    FLA_Part_1x2( UB,     &UB_l,  &none,    b, FLA_LEFT );
    FLA_Part_1x2( VB,     &VB_l,  &none,    b, FLA_LEFT );
    FLA_Part_1x2( YB,     &YB_l,  &none,    b, FLA_LEFT );
    FLA_Part_1x2( ZB,     &ZB_l,  &none,    b, FLA_LEFT );

    FLA_Part_2x1( UB_l,   &none,
                          &U2_l,            b, FLA_TOP );
    FLA_Part_2x1( VB_l,   &none,
                          &V2_l,            b, FLA_TOP );
    FLA_Part_2x1( YB_l,   &none, 
                          &Y2_l,            b, FLA_TOP );
    FLA_Part_2x1( ZB_l,   &none,
                          &Z2_l,            b, FLA_TOP );

    // [ ABR, YB, ZB, TU1, TV1 ] = FLA_Bidiag_UT_u_step_unb_var5( ABR, TU1, TV1, b );
    //FLA_Bidiag_UT_u_step_unb_var5( ABR, YB, ZB, TU1_tl, TV1_tl );
    FLA_Bidiag_UT_u_step_opt_var5( ABR, YB, ZB, TU1_tl, TV1_tl );

    if ( FLA_Obj_length( A22 ) > 0 )
    {
      // Build UB from ABR, with explicit unit subdiagonal and zeros.
      FLA_Copy( ABR_l, UB_l );
      FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_UNIT_DIAG, UB_l );

      // Build VB from ABR, with explicit unit subdiagonal and zeros.
      FLA_Copyt( FLA_TRANSPOSE, ABR_t, VB_l );
      FLA_Part_2x1( VB_l,   &VB_tl,
                            &VB_bl,            1, FLA_TOP );
      FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_UNIT_DIAG, VB_bl );
      FLA_Set( FLA_ZERO, VB_tl );

      // A22 = A22 - U2 * Y2' - Z2 * V2';
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                         FLA_MINUS_ONE, U2_l, Y2_l, FLA_ONE, A22 );
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                         FLA_MINUS_ONE, Z2_l, V2_l, FLA_ONE, A22 );
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
    FLA_Cont_with_3x1_to_2x1( &VT,                V0,
                                                  V1,
                            /* ** */           /* ** */
                              &VB,                V2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &YT,                Y0,
                                                  Y1,
                            /* ** */           /* ** */
                              &YB,                Y2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &ZT,                Z0,
                                                  Z1,
                            /* ** */           /* ** */
                              &ZB,                Z2,     FLA_TOP );
    FLA_Cont_with_1x3_to_1x2( &TUL, /**/ &TUR,       TU0, TU1, /**/ TU2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &TVL, /**/ &TVR,       TV0, TV1, /**/ TV2,
                              FLA_LEFT );
  }

  FLA_Obj_free( &U );
  FLA_Obj_free( &V );
  FLA_Obj_free( &Y );
  FLA_Obj_free( &Z );

  return FLA_SUCCESS;
}

