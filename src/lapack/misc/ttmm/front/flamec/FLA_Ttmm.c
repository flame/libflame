/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_ttmm_t*      fla_ttmm_cntl_leaf;
extern __thread fla_ttmm_t*      fla_ttmm_cntl;
extern __thread fla_blocksize_t* fla_ttmm_var1_bsize;

FLA_Error FLA_Ttmm( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Datatype datatype;
  int          m_A, r_val = 0;
  int          FLA_TTMM_VAR1_BLOCKSIZE;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Ttmm_check( uplo, A );

  // Determine the datatype of the operation.
  datatype = FLA_Obj_datatype( A );

  // Extract the appropriate blocksize for the given datatype.
  FLA_TTMM_VAR1_BLOCKSIZE = FLA_Blocksize_extract( datatype, fla_ttmm_var1_bsize );

  // Determine the dimension of A.
  m_A = FLA_Obj_length( A );

  // Invoke FLA_Ttmm_internal() with the appropriate control tree.
  if      ( m_A <= FLA_TTMM_VAR1_BLOCKSIZE )
  {
    r_val = FLA_Ttmm_internal( uplo, A, fla_ttmm_cntl_leaf );
  }
  else if ( FLA_TTMM_VAR1_BLOCKSIZE < m_A )
  {
    r_val = FLA_Ttmm_internal( uplo, A, fla_ttmm_cntl );
  }

  return r_val;
}

