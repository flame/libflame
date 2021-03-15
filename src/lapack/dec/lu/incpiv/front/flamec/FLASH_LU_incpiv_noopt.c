/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_lu_t* flash_lu_incpiv_cntl;

FLA_Error FLASH_LU_incpiv_noopt( FLA_Obj A, FLA_Obj p, FLA_Obj L )
{
  dim_t     nb_alg;
  FLA_Error r_val;
  
  // Inspect the width of a the top-left element of L to get the algorithmic
  // blocksize we'll use throughout the LU_incpiv algorithm.
  nb_alg = FLASH_Obj_scalar_width_tl( L );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLASH_LU_incpiv_var1( A, p, L, nb_alg, flash_lu_incpiv_cntl );
  
  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}
