/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_trinv_t*     fla_trinv_cntl_leaf;
extern fla_trinv_t*     fla_trinv_cntl;
extern fla_blocksize_t* fla_trinv_var3_bsize;

FLA_Error FLA_Trinv( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A )
{
  FLA_Datatype datatype;
  int          m_A, r_val = 0;
  int          FLA_TRINV_VAR3_BLOCKSIZE;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Trinv_check( uplo, diag, A );

  // Determine the datatype of the operation.
  datatype = FLA_Obj_datatype( A );

  // Extract the appropriate blocksize for the given datatype.
  FLA_TRINV_VAR3_BLOCKSIZE = FLA_Blocksize_extract( datatype, fla_trinv_var3_bsize );

  // Determine the dimension of A.
  m_A = FLA_Obj_length( A );

  // Invoke FLA_Trinv_internal() with the appropriate control tree.
  if      ( m_A <= FLA_TRINV_VAR3_BLOCKSIZE )
  {
    r_val = FLA_Trinv_internal( uplo, diag, A, fla_trinv_cntl_leaf );
  }
  else if ( FLA_TRINV_VAR3_BLOCKSIZE < m_A )
  {
    r_val = FLA_Trinv_internal( uplo, diag, A, fla_trinv_cntl );
  }

  return r_val;
}

