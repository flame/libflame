/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Gemm_nn_omp_var13( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj C, fla_gemm_t* cntl )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  FLA_Obj BL,    BR,       B0,  B1,  B2;

  FLA_Obj C1L,   C1R,      C10,  C11,  C12;

  int b_m, b_n;
 

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_TOP );


  #pragma intel omp parallel taskq
  {
  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) )
  {
    b_m = FLA_Determine_blocksize( A, AT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* ** */
                                              &A1, 
                           AB,                &A2,        b_m, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                        /* ** */            /* ** */
                                              &C1, 
                           CB,                &C2,        b_m, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

    FLA_Part_1x2( C1,    &C1L,  &C1R,      0, FLA_LEFT );

    while ( FLA_Obj_width( BL ) < FLA_Obj_width( B ) ){

      b_n = FLA_Determine_blocksize( B, BL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

      FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                             b_n, FLA_RIGHT );

      FLA_Repart_1x2_to_1x3( C1L,  /**/ C1R,        &C10, /**/ &C11, &C12,
                             b_n, FLA_RIGHT );

      /*------------------------------------------------------------*/

      #pragma intel omp task captureprivate( A1, B1, C11 )
      {
      /*    C1 = alpha * A * B1 + C1; */
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                         alpha, A1, B1, FLA_ONE, C11 );
      }

      /*------------------------------------------------------------*/

      FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                                FLA_LEFT );

      FLA_Cont_with_1x3_to_1x2( &C1L,  /**/ &C1R,        C10, C11, /**/ C12,
                                FLA_LEFT );

    }


    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  A1, 
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                                                  C1, 
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );

  }
  }

  return FLA_SUCCESS;
}

