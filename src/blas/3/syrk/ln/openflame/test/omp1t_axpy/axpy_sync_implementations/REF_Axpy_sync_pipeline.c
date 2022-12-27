/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef OPENMP
FLA_Error REF_Axpy_sync_pipeline( FLA_Obj alpha, FLA_Obj X, FLA_Obj B )
{
  double* x_buf, *b_buf;
  integer     x_m, x_n, x_ldim, b_ldim;
  integer     b, j, j2, j_part, nb_alg;
  int     i_one = 1;
  double  alpha_value;

  // Compute the width of one lockable partition.
  nb_alg    = FLA_omp_compute_stage_width( X );
  x_n       = FLA_Obj_width( X );
  x_m       = FLA_Obj_length( X );
  x_buf     = FLA_Obj_buffer_at_view( X );
  b_buf     = FLA_Obj_buffer_at_view( B );
  x_ldim    = FLA_Obj_ldim( X );
  b_ldim    = FLA_Obj_ldim( B );
  alpha_value = FLA_DOUBLE_VALUE( alpha );

  for( j = 0; j < x_n; j += nb_alg )
  {
    b = fla_min( x_n-j, nb_alg );

    /*------------------------------------------------------------*/

    // Get the index of the current partition.
    j_part = j/nb_alg;

    // Acquire lock[j_part] (the lock for X1 and B1).
    omp_set_lock( &fla_omp_lock[j_part] );

    // B1 := alpha * X1 + B1
    for( j2 = 0; j2 < b; ++j2 )
      FLA_C2F( daxpy) ( &x_m,
                        &alpha_value,
                        x_buf + (j+j2)*x_ldim, &i_one,
                        b_buf + (j+j2)*b_ldim, &i_one );

    // Release lock[j_part] (the lock for X1 and B1).
    omp_unset_lock( &fla_omp_lock[j_part] );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}
#endif
