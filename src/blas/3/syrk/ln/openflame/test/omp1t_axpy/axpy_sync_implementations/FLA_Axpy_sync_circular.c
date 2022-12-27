/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef OPENMP
FLA_Error FLA_Axpy_sync_circular( FLA_Obj alpha, FLA_Obj X, FLA_Obj B )
{
  FLA_Obj XL,    XR,       X0,  X1,  X2;
  FLA_Obj BL,    BR,       B0,  B1,  B2;

  int n_stages    = FLA_omp_get_num_stages();
  int stage_width = FLA_omp_compute_stage_width( X );
  int thread_num  = omp_get_thread_num();
  int n_done      = 0;
  int b, i;

  // Start thread i on the ith panel partition of B.
  FLA_Part_1x2( X,    &XL,  &XR,    stage_width*thread_num, FLA_LEFT );
  FLA_Part_1x2( B,    &BL,  &BR,    stage_width*thread_num, FLA_LEFT );

  while ( n_done++ < n_stages ){

    // The last lockable partition may be smaller than the others.
    b = fla_min( FLA_Obj_width( XR ), stage_width );
    
    FLA_Repart_1x2_to_1x3( XL,  /**/ XR,        &X0, /**/ &X1, &X2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    // Get the index of the current partition.
    i = FLA_Obj_width(XL)/stage_width;

    // Acquire lock[i] (the lock for X1 and B1).
    omp_set_lock( &fla_omp_lock[i] );

    // B1 := alpha * X1 + B1
    FLA_Axpy_external( alpha, X1, B1 );

    // Release lock[i] (the lock for X1 and B1).
    omp_unset_lock( &fla_omp_lock[i] );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &XL,  /**/ &XR,        X0, X1, /**/ X2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );

    // If this thread reaches the last partition, wrap back around to 
    // the first partition for the next iteration.
    if( FLA_Obj_width( XL ) == FLA_Obj_width( X ) )
    {
      FLA_Part_1x2( X,    &XL,  &XR,      0, FLA_LEFT );
      FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );
    }

  }

  return FLA_SUCCESS;
}
#endif
