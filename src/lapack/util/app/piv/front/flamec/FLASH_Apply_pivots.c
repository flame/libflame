
#include "FLAME.h"

extern fla_appiv_t* flash_appiv_cntl;

FLA_Error FLASH_Apply_pivots( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A )

{
  FLA_Error r_val;
  FLA_Bool  enable_supermatrix;

  // Check parameters.

  // *** The current Apply_pivots algorithm implemented assumes that
  // the matrix has a hierarchical depth of 1. We check for that here, because
  // we anticipate that we'll use a more general algorithm in the future, and
  // we don't want to forget to remove the constraint. ***
  if ( FLASH_Obj_depth( A ) != 1 )
  {
    FLA_Print_message( "FLASH_Apply_pivots() currently only supports matrices of depth 1",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  // Find the status of SuperMatrix.
  enable_supermatrix = FLASH_Queue_get_enabled();

  // Temporarily disable SuperMatrix.
  FLASH_Queue_disable();

  // Invoke FLA_Apply_pivots_internal() with large control tree.
  r_val = FLA_Apply_pivots_internal( side, trans, p, A, flash_appiv_cntl );

  // Restore SuperMatrix to its previous status.
  if ( enable_supermatrix )
     FLASH_Queue_enable();

  return r_val;
}

