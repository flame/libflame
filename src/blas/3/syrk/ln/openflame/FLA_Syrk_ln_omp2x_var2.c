
#include "FLAME.h"
#include "FLA_Syrk_ln_omp.h"

FLA_Error FLA_Syrk_ln_omp2x_var2( FLA_Obj A, FLA_Obj C )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj CTL,   CTR,      C00, C01, C02, 
          CBL,   CBR,      C10, C11, C12,
                           C20, C21, C22;

  FLA_Obj A2_T,            A2_0,
          A2_B,            A2_1,
                           A2_2;

  FLA_Obj C21_T,           C21_0,
          C21_B,           C21_1,
                           C21_2;
  int b, b2;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_TL );

  #pragma intel omp parallel taskq
  {
  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    b = FLA_Task_compute_blocksize( 0, A, AT, FLA_TOP );

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

    FLA_Part_2x1( A2,     &A2_T, 
                          &A2_B,            0, FLA_TOP );

    FLA_Part_2x1( C21,    &C21_T, 
                          &C21_B,            0, FLA_TOP );

    while ( FLA_Obj_length( C21_T ) < FLA_Obj_length( C21 ) ){

      b2 = FLA_Task_compute_blocksize( 1, C21, C21_T, FLA_TOP );

      FLA_Repart_2x1_to_3x1( A2_T,              &A2_0, 
                           /* ** */           /* ** */
                                                &A2_1, 
                             A2_B,              &A2_2,        b2, FLA_BOTTOM );

      FLA_Repart_2x1_to_3x1( C21_T,             &C21_0, 
                           /* ** */            /* ** */
                                                &C21_1, 
                             C21_B,             &C21_2,       b2, FLA_BOTTOM );

      /*------------------------------------------------------------*/

      #pragma intel omp task captureprivate(A2_1, A1, C21_1)
      {
      /*    C1 = alpha * A1 * B' + C1;   */
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A2_1, A1, FLA_ONE, C21_1 );
      }

      /*------------------------------------------------------------*/

      FLA_Cont_with_3x1_to_2x1( &A2_T,              A2_0, 
                                                    A2_1, 
                               /* ** */            /* ** */
                                &A2_B,              A2_2,     FLA_TOP );

      FLA_Cont_with_3x1_to_2x1( &C21_T,             C21_0, 
                                                    C21_1, 
                               /* ** */            /* ** */
                                &C21_B,             C21_2,    FLA_TOP );

    }
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

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    b = FLA_Task_compute_blocksize( 0, A, AT, FLA_TOP );

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

    #pragma intel omp task captureprivate(C11, A1)
    {
    /* C11 = C11 + A1 * A1' */
    FLA_Syrk_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_ONE, A1, FLA_ONE, C11 );
    }

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
  }

  return FLA_SUCCESS;
}

