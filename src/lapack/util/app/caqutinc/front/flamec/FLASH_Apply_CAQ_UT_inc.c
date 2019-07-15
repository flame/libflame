/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_apcaqutinc_t* flash_apcaqutinc_cntl;

FLA_Error FLASH_Apply_CAQ_UT_inc( dim_t p,
                                  FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                  FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW, FLA_Obj W, FLA_Obj B )
{
  FLA_Error r_val;
  dim_t     nb_part;
  FLA_Obj   WT, WB;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Apply_CAQ_UT_inc_check( side, trans, direct, storev, A, ATW, R, RTW, W, B );

  // Compute the partition length from the number of partitions.
  nb_part = FLA_CAQR_UT_inc_compute_blocks_per_part( p, R );

  // Begin a parallel region.
  FLASH_Queue_begin();

  // Apply the individual Q's from the incremental QR factorizations.
  FLA_Apply_CAQ_UT_inc_apply_panels( nb_part, A, ATW, W, B );

  FLA_Part_2x1( W,   &WT,
                     &WB,    1, FLA_TOP );

  // Apply the Q from the factorization of the upper triangular R's.
  r_val = FLA_Apply_CAQ_UT_inc_internal( side, trans, direct, storev,
                                         R, RTW, WT, B, flash_apcaqutinc_cntl );


  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

