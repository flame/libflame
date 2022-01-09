/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Gemm_nn_omp_var35( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj C, fla_gemm_t* cntl )
{
  FLA_Obj BL,    BR,       B0,  B1,  B2;

  FLA_Obj CL,    CR,       C0,  C1,  C2;

  FLA_Obj AL,    AR,       A0,  A1,  A2;

  FLA_Obj BT,              B01,
          BB,              B11,
                           B21;
  FLA_Obj C1_local;

  integer i, j, lock_ldim, lock_i;
  integer b_n, b_k;


  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  #pragma intel omp parallel taskq
  {
  while ( FLA_Obj_width( BL ) < FLA_Obj_width( B ) )
  {
    b_n = FLA_Determine_blocksize( B, BL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b_n, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                           b_n, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

    FLA_Part_2x1( B1,   &BT, 
                        &BB,            0, FLA_TOP );
  
    while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) )
    {
      b_k = FLA_Determine_blocksize( A, AL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

      // Get the index of the current partition.
      //j = FLA_Obj_width( CL ) / b_n;
      //i = FLA_Obj_width( AL ) / b_k;
      //lock_ldim = FLA_get_num_threads_in_n_dim(omp_get_num_threads());
      lock_i = FLA_Obj_width( CL ) / b_n;
  
      FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                             b_k, FLA_RIGHT );

      FLA_Repart_2x1_to_3x1( BT,                &B01,
                          /* ** */            /* ** */
                                                &B11,
                             BB,                &B21,       b_k, FLA_BOTTOM );
  
      /*------------------------------------------------------------*/
  
      /*    C1 = alpha * A1 * B11 + C1; */
      // FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
      //           alpha, A1, B11, FLA_ONE, C1 );

      #pragma intel omp task captureprivate( lock_i, A1, B11, C1 ) private( C1_local )
      {
      FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C1, &C1_local );
      FLA_Obj_set_to_zero( C1_local );

      // C1_local = alpha * A1 * B11 + C1_local;
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                         alpha, A1, B11, FLA_ONE, C1_local );

      // Acquire lock[i] (the lock for C1).
      omp_set_lock( &fla_omp_lock[lock_i] );

      // C1 += C1_local
      FLA_Axpy_external( FLA_ONE, C1_local, C1 );
      //FLA_Axpy_sync_pipeline2( j*lock_ldim, FLA_ONE, C1_local, C1 );

      // Release lock[i] (the lock for C1).
      omp_unset_lock( &fla_omp_lock[lock_i] );

      FLA_Obj_free( &C1_local );
      }

      /*------------------------------------------------------------*/
  
      FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                                FLA_LEFT );

      FLA_Cont_with_3x1_to_2x1( &BT,                B01, 
                                                    B11, 
                              /* ** */           /* ** */
                                &BB,                B21,    FLA_TOP );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                              FLA_LEFT );
  }
  }

  return FLA_SUCCESS;
}

