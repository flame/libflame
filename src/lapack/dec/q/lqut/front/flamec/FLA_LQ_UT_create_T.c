/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_LQ_UT_create_T( FLA_Obj A, FLA_Obj* T )
{
  FLA_Datatype datatype;
  dim_t        b_alg, k;
  dim_t        rs_T, cs_T;

  // Query the datatype of A.
  datatype = FLA_Obj_datatype( A );

  // Query the blocksize from the library.
  b_alg = FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  // Scale the blocksize by a pre-set global constant.
  b_alg = ( dim_t )( ( ( double ) b_alg ) * FLA_LQ_INNER_TO_OUTER_B_RATIO );

  // Adjust the blocksize with respect to the min-dim of A.
  b_alg = fla_min(b_alg, FLA_Obj_min_dim( A ));

  // Query the length of A.
  k = FLA_Obj_length( A );

  // Figure out whether T should be row-major or column-major.
  if ( FLA_Obj_row_stride( A ) == 1 )
  {
    rs_T = 1;
    cs_T = b_alg;
  }
  else // if ( FLA_Obj_col_stride( A ) == 1 )
  {
    rs_T = k;
    cs_T = 1;
  }

  // Create a b_alg x k matrix to hold the block Householder transforms that
  // will be accumulated within the LQ factorization algorithm.
  FLA_Obj_create( datatype, b_alg, k, rs_T, cs_T, T );

  return FLA_SUCCESS;
}

