/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_gemm_t* fla_gemm_cntl_mm_op;

FLA_Error FLA_Gemm( FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Error r_val = FLA_SUCCESS;
  
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Gemm_check( transa, transb, alpha, A, B, beta, C );

#ifdef FLA_ENABLE_BLAS3_FRONT_END_CNTL_TREES
  r_val = FLA_Gemm_internal( transa, transb, alpha, A, B, beta, C, fla_gemm_cntl_mm_op );
#else
  r_val = FLA_Gemm_external( transa, transb, alpha, A, B, beta, C );
#endif

  return r_val;
}

