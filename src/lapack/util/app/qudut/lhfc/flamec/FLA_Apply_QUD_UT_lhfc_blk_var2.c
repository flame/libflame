
#include "FLAME.h"

FLA_Error FLA_Apply_QUD_UT_lhfc_blk_var2( FLA_Obj T, FLA_Obj W,
                                                     FLA_Obj R,
                                          FLA_Obj U, FLA_Obj C,
                                          FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl )
{
  FLA_Obj WL,    WR,       W0,  W1,  W2;

  FLA_Obj RL,    RR,       R0,  R1,  R2;

  FLA_Obj CL,    CR,       C0,  C1,  C2;

  FLA_Obj DL,    DR,       D0,  D1,  D2;

  dim_t b;

  FLA_Part_1x2( W,    &WL,  &WR,      0, FLA_LEFT );

  FLA_Part_1x2( R,    &RL,  &RR,      0, FLA_LEFT );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );

  while ( FLA_Obj_width( RL ) < FLA_Obj_width( R ) ){

    b = FLA_Determine_blocksize( RR, FLA_RIGHT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( WL,  /**/ WR,        &W0, /**/ &W1, &W2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( RL,  /**/ RR,        &R0, /**/ &R1, &R2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    //  Apply Q' to R1, C1, and D1 from the left:
    //
    //  / R1 \       / R1 \
    //  | C1 | =  Q' | C1 |
    //  \ D1 /       \ D1 /
    //
    //  where Q is formed from U, V, and T.

    FLA_Apply_QUD_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                T, W1,
                                   R1,
                                U, C1,
                                V, D1, FLA_Cntl_sub_apqudut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &WL,  /**/ &WR,        W0, W1, /**/ W2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &RL,  /**/ &RR,        R0, R1, /**/ R2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

