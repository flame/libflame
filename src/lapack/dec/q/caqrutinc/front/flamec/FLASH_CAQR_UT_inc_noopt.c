/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_caqrutinc_t* flash_caqrutinc_cntl;

FLA_Error FLASH_CAQR_UT_inc_noopt( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW )
{
  FLA_Error r_val = FLA_SUCCESS;
  dim_t     nb_part;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_CAQR_UT_inc_check( p, A, ATW, R, RTW );

  // Compute the partition length from the number of partitions.
  nb_part = FLA_CAQR_UT_inc_compute_blocks_per_part( p, A );

  // Begin a parallel region.
  FLASH_Queue_begin();

  // Perform incremental QR's on each of the p partitions.
  FLA_CAQR_UT_inc_factorize_panels( nb_part, A, ATW );

  // Copy the triangles of A into R.
  FLA_CAQR_UT_inc_copy_triangles( nb_part, A, R );

  // Perform an incremental CAQR on the resulting upper triangular R's in A.
  FLA_CAQR_UT_inc_blk_var1( R, RTW, flash_caqrutinc_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

