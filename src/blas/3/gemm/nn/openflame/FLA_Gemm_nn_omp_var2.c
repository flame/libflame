
#include "FLAME.h"

FLA_Error FLA_Gemm_nn_omp_var2( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj C, fla_gemm_t* cntl )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  int b;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_BOTTOM );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_BOTTOM );

  #pragma intel omp parallel taskq
  {
  while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){

    b = FLA_Determine_blocksize( A, AB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                              &A1, 
                        /* ** */            /* ** */
                           AB,                &A2,        b, FLA_TOP );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                                              &C1, 
                        /* ** */            /* ** */
                           CB,                &C2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    #pragma intel omp task captureprivate( A1, C1 )
    {
    /* C1 = alpha * A1 * B + C1; */
    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       alpha, A1, B, FLA_ONE, C1 );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                            /* ** */           /* ** */
                                                  A1, 
                              &AB,                A2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                            /* ** */           /* ** */
                                                  C1, 
                              &CB,                C2,     FLA_BOTTOM );

  }
  }

  return FLA_SUCCESS;
}

