
#include "FLAME.h"

FLA_Error FLA_Apply_QUD_UT_inc_lhfc_blk_var1( FLA_Obj T, FLA_Obj W,
                                                         FLA_Obj B,
                                              FLA_Obj U, FLA_Obj C,
                                              FLA_Obj V, FLA_Obj D, fla_apqudutinc_t* cntl )
{
  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj UL,    UR,       U0,  U1,  U2;

  FLA_Obj VL,    VR,       V0,  V1,  V2;

  FLA_Obj WT,              W0,
          WB,              W1,
                           W2;

  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  dim_t   b;

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  FLA_Part_1x2( U,    &UL,  &UR,      0, FLA_LEFT );

  FLA_Part_1x2( V,    &VL,  &VR,      0, FLA_LEFT );

  FLA_Part_2x1( W,    &WT,
                      &WB,            0, FLA_TOP );

  FLA_Part_2x1( B,    &BT,
                      &BB,            0, FLA_TOP );

  while ( FLA_Obj_width( UL ) < FLA_Obj_width( U ) ){

    b = FLA_Determine_blocksize( UR, FLA_RIGHT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &T2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( UL,  /**/ UR,        &U0, /**/ &U1, &U2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( VL,  /**/ VR,        &V0, /**/ &V1, &V2,
                           b, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( WT,                &W0,
                        /* ** */            /* ** */
                                              &W1,
                           WB,                &W2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( BT,                &B0,
                        /* ** */            /* ** */
                                              &B1,
                           BB,                &B2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    /*
         Apply Q' to B1, C, and D from the left:

           / B1 \          / B1 \
           | C  |  =  Q' * | C  |
           \ D  /          \ D  /

         where Q is formed from U1, V1, and T1.
    */

    FLA_Apply_QUD_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                               T1, W1,
                                   B1,
                               U1, C,
                               V1, D, FLA_Cntl_sub_apqudut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ T2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &UL,  /**/ &UR,        U0, U1, /**/ U2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &VL,  /**/ &VR,        V0, V1, /**/ V2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &WT,                W0,
                                                  W1,
                            /* ** */           /* ** */
                              &WB,                W2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                                                  B1,
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

