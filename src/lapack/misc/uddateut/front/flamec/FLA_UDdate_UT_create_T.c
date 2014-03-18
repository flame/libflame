
#include "FLAME.h"

FLA_Error FLA_UDdate_UT_create_T( FLA_Obj R, FLA_Obj* T )
{
  FLA_Datatype datatype;
  dim_t        b_alg, k;
  dim_t        rs_T, cs_T;

  // Query the datatype of R.
  datatype = FLA_Obj_datatype( R );

  // Query the blocksize from the library.
  b_alg = FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  // We want the width of R, which is the same as that of C and D. Also,
  // R should be square, so we could grab either dimension.
  k = FLA_Obj_width( R );

  // Figure out whether T should be row-major or column-major.
  if ( FLA_Obj_row_stride( R ) == 1 )
  {
    rs_T = 1;
    cs_T = b_alg;
  }
  else // if ( FLA_Obj_col_stride( R ) == 1 )
  {
    rs_T = k;
    cs_T = 1;
  }

  // Create a b_alg x k matrix to hold the block Householder transforms that
  // will be accumulated within the UDdate operation algorithm.
  FLA_Obj_create( datatype, b_alg, k, rs_T, cs_T, T );

  return FLA_SUCCESS;
}

