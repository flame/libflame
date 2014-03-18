
#include "FLAME.h"

extern fla_lu_t* flash_lu_piv_cntl;

FLA_Error FLASH_LU_piv( FLA_Obj A, FLA_Obj p )
{
  FLA_Error r_val = FLA_SUCCESS;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LU_piv_check( A, p );

  // *** The current LU_piv algorithm implemented assumes that
  // the matrix has a hierarchical depth of 1. We check for that here, because
  // we anticipate that we'll use a more general algorithm in the future, and
  // we don't want to forget to remove the constraint. ***
  if ( FLASH_Obj_depth( A ) != 1 )
  {
    FLA_Print_message( "FLASH_LU_piv() currently only supports matrices of depth 1",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  // Begin a parallel region.
  FLASH_Queue_begin();

  // Invoke FLA_LU_piv_internal() with large control tree.
  FLA_LU_piv_internal( A, p, flash_lu_piv_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  // Check for singularity.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    r_val = FLASH_LU_find_zero_on_diagonal( A );

  return r_val;
}

