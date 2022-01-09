/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Gemm_nn_omp_var51( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj C, fla_gemm_t* cntl )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;

  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  FLA_Obj AT,              A01,
          AB,              A11,
                           A21;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;
  FLA_Obj C1_local;

  integer b_m, b_k;
  integer i;


  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_TOP );

  #pragma intel omp parallel taskq
  {
  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) )
  {
    b_k = FLA_Determine_blocksize( A, AL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );
    //b_c = min( FLA_Obj_width( AR ), nb_alg_c );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b_k, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                        /* ** */            /* ** */
                                              &B1, 
                           BB,                &B2,        b_k, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( A1,   &AT, 
                        &AB,            0, FLA_TOP );

    FLA_Part_2x1( C,    &CT, 
                        &CB,            0, FLA_TOP );
  
    while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) )
    {
      b_m = FLA_Determine_blocksize( A, AT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );
      //b_r = min( FLA_Obj_length( AB ), nb_alg_r );

      // Get the index of the current partition.
      i = FLA_Obj_length( CT ) / b_m;
  
      FLA_Repart_2x1_to_3x1( AT,                &A01, 
                          /* ** */            /* ** */
                                                &A11, 
                             AB,                &A21,      b_m, FLA_BOTTOM );

      FLA_Repart_2x1_to_3x1( CT,                &C0, 
                          /* ** */            /* ** */
                                                &C1,
                             CB,                &C2,       b_m, FLA_BOTTOM );
  
      /*------------------------------------------------------------*/
  
      /* C1 = alpha * A11 * B1 + C1; */
      //// FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
      ////          alpha, A11, B1, FLA_ONE, C1 );
  
      #pragma intel omp task captureprivate( i, A11, B1, C1 ) private( C1_local )
      {
      FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C1, &C1_local );
      FLA_Obj_set_to_zero( C1_local );

      /*    C1_local = alpha * A11 * B1 + C1_local; */
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                         alpha, A11, B1, FLA_ONE, C1_local );

      // Acquire lock[i] (the lock for C1).
      omp_set_lock( &fla_omp_lock[i] );

      /* C1 += C1_local */
      FLA_Axpy_external( FLA_ONE, C1_local, C1 );

      // Release lock[i] (the lock for C1).
      omp_unset_lock( &fla_omp_lock[i] );

      FLA_Obj_free( &C1_local );
      }
  
      /*------------------------------------------------------------*/
  
      FLA_Cont_with_3x1_to_2x1( &AT,                A01, 
                                                    A11, 
                              /* ** */           /* ** */
                                &AB,                A21,    FLA_TOP );

      FLA_Cont_with_3x1_to_2x1( &CT,                C0,
                                                    C1,
                              /* ** */           /* ** */
                                &CB,                C2,     FLA_TOP );
    }
  
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                                                  B1, 
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );
  }
  }

  return FLA_SUCCESS;
}

