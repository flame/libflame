/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_FS_incpiv( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj b )
{
  dim_t     nb_alg;
  FLA_Error r_val;
  FLA_Bool  enable_supermatrix;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_FS_incpiv_check( A, p, L, b );

  // *** The current forward substitution algorithm implemented assumes that 
  // the matrix has a hierarchical depth of 1. We check for that here, because
  // we anticipate that we'll use a more general algorithm in the future, and 
  // we don't want to forget to remove the constraint. ***
  if ( FLASH_Obj_depth( A ) != 1 )
  {
    FLA_Print_message( "FLASH_FS_incpiv() currently only supports matrices of depth 1",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }
  
  // Inspect the width of a the top-left element of L to get the algorithmic
  // blocksize we'll use throughout the LU_incpiv algorithm.
  nb_alg = FLASH_Obj_scalar_width_tl( L );

  // Find the status of SuperMatrix.
  enable_supermatrix = FLASH_Queue_get_enabled();

  // Temporarily disable SuperMatrix.
  FLASH_Queue_disable();
  
  // Execute tasks.
  r_val = FLASH_FS_incpiv_aux1( A, p, L, b, nb_alg );

  // Restore SuperMatrix to its previous status.
  if ( enable_supermatrix )
     FLASH_Queue_enable();

  return r_val;
}

