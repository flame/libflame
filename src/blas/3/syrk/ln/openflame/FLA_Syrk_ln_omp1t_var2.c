
#include "FLAME.h"
#include "FLA_Syrk_ln_omp.h"

FLA_Error FLA_Syrk_ln_omp1t_var2( FLA_Obj A, FLA_Obj C )
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

    #pragma intel omp task captureprivate(C11, A1, A2, C21)
    {
    /* C21 = C21 + A2 * A1' */
    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A2, A1, FLA_ONE, C21 );

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

