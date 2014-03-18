
#include "FLAME.h"

#ifdef OPENMP
FLA_Error FLA_Axpy_sync_pipeline( FLA_Obj alpha, FLA_Obj X, FLA_Obj B )
{
  FLA_Obj XL,    XR,       X0,  X1,  X2;
  FLA_Obj BL,    BR,       B0,  B1,  B2;

  int b, i, nb_alg;

  FLA_Part_1x2( X,    &XL,  &XR,      0, FLA_LEFT );
  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  // Compute the width of one lockable partition.
  nb_alg = FLA_omp_compute_stage_width( X );

  while ( FLA_Obj_width( XL ) < FLA_Obj_width( X ) ){

    b = min( FLA_Obj_width( XR ), nb_alg );

    FLA_Repart_1x2_to_1x3( XL,  /**/ XR,        &X0, /**/ &X1, &X2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    // Get the index of the current partition.
    i = FLA_Obj_width(XL)/nb_alg;

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

  }

  return FLA_SUCCESS;
}
#endif
