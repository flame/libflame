#include "FLAME.h"
#include "Gemm_prototypes.h"

int Gemm_unb_var3( FLA_Obj A, FLA_Obj B, FLA_Obj C )
{
  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;

  FLA_Obj CT,              C0,
          CB,              c1t,
                           C2;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* *** */
                                              &a1t, 
                           AB,                &A2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                        /* ** */            /* *** */
                                              &c1t, 
                           CB,                &C2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Gemv( FLA_TRANSPOSE, FLA_ONE, B, a1t, FLA_ONE, c1t );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  a1t, 
                            /* ** */           /* *** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                                                  c1t, 
                            /* ** */           /* *** */
                              &CB,                C2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}

