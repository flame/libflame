/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_qrut_t*   flash_qrut_cntl;
extern TLS_CLASS_SPEC fla_qrut_t*   fla_qrut_cntl_leaf;

FLA_Error FLASH_QR_UT( FLA_Obj A, FLA_Obj TW )
{
  FLA_Error r_val;
  dim_t     b_alg, b_flash;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_QR_UT_check( A, TW );

  // *** The current hierarchical QR_UT algorithm assumes that the matrix
  // has a hierarchical depth of 1. We check for that here, because we
  //  anticipate that we'll use a more general algorithm in the future, and
  // we don't want to forget to remove the constraint. ***
  if ( FLASH_Obj_depth( A ) != 1 )
  {
    FLA_Print_message( "FLASH_QR_UT() currently only supports matrices of depth 1",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  // Inspect the length of TTL to get the blocksize used by the QR
  // factorization, which will be our inner blocksize for Apply_Q_UT.
  b_alg   = FLASH_Obj_scalar_length_tl( TW );
  b_flash = FLASH_Obj_scalar_width_tl( TW );

  // The traditional (non-incremental) QR_UT algorithm-by-blocks requires
  // that the algorithmic blocksize be equal to the storage blocksize.
  if ( b_alg != b_flash )
  {
    FLA_Print_message( "FLASH_QR_UT() requires that b_alg == b_store",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  // The traditional (non-incremental) QR_UT algorithm-by-blocks requires
  // that min_dim(A) % b_flash == 0.
  if ( FLASH_Obj_scalar_min_dim( A ) % b_flash != 0 )
  {
    FLA_Print_message( "FLASH_QR_UT() requires that min_dim( A ) %% b_store == 0",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  // Begin a parallel region.
  FLASH_Queue_begin();

  // Invoke FLA_QR_UT_internal() with hierarchical control tree.
  r_val = FLA_QR_UT_internal( A, TW, flash_qrut_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

