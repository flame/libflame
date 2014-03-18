
#include "FLAME.h"

FLA_Error FLA_Syrk_ln_blk_var2_ht( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, FLA_Syrk_t* cntl )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj CTL,   CTR,      C00, C01, C02, 
          CBL,   CBR,      C10, C11, C12,
                           C20, C21, C22;

  int b;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    b = 1;

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* ** */
                                              &A1, 
                           AB,                &A2,        b, FLA_BOTTOM );

    FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00, /**/ &C01, &C02,
                        /* ************* */   /* ******************** */
                                                &C10, /**/ &C11, &C12,
                           CBL, /**/ CBR,       &C20, /**/ &C21, &C22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    /* C21 = C21 + A2 * A1' */
    FLA_Gemm_nt_blk_var1_ht( alpha, A2, A1, beta, C21,
                          NULL );

    /* C11 = C11 + A1 * A1' */
    FLA_Syrk_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, 
                       alpha, *FLASH_OBJ_PTR_AT( A1 ), beta, *FLASH_OBJ_PTR_AT( C11 ),
                       NULL );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  A1, 
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00, C01, /**/ C02,
                                                     C10, C11, /**/ C12,
                            /* ************** */  /* ****************** */
                              &CBL, /**/ &CBR,       C20, C21, /**/ C22,
                              FLA_TL );

  }

  return FLA_SUCCESS;
}

