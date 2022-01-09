/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_apqut_t* flash_apqut_cntl_blas;
extern TLS_CLASS_SPEC fla_apqut_t* fla_apqut_cntl_leaf;

FLA_Error FLASH_Apply_Q_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B )
{
  FLA_Error r_val;
  dim_t     b_alg;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Apply_Q_UT_check( side, trans, direct, storev, A, T, W, B );

  // Inspect the length of TTL to get the blocksize used by the QR/LQ
  // factorization, which will be our inner blocksize for Apply_Q_UT.
  b_alg = FLASH_Obj_scalar_length_tl( T );

  // The traditional (non-incremental) Apply_Q_UT algorithm-by-blocks
  // requires that the algorithmic blocksize be equal to the storage
  // blocksize.
  if ( b_alg != FLASH_Obj_scalar_width_tl( T ) )
  {
    FLA_Print_message( "FLASH_Apply_Q_UT() requires that b_alg == b_store",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  // Adjust the blocksize of the control tree node for the flat subproblem.
  if ( FLA_Cntl_blocksize( fla_apqut_cntl_leaf ) != NULL )
    FLA_Blocksize_set( FLA_Cntl_blocksize( fla_apqut_cntl_leaf ),
                       b_alg, b_alg, b_alg, b_alg );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Invoke FLA_Apply_Q_UT_internal() with the standard control tree.
  r_val = FLA_Apply_Q_UT_internal( side, trans, direct, storev, A, T, W, B,
                                   flash_apqut_cntl_blas );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

