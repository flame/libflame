
#include "FLAME.h"

extern fla_apqudut_t*    fla_apqudut_cntl_leaf;
extern fla_apqudutinc_t* flash_apqudutinc_cntl;

FLA_Error FLASH_Apply_QUD_UT_inc( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                  FLA_Obj T, FLA_Obj W,
                                             FLA_Obj R,
                                  FLA_Obj U, FLA_Obj C,
                                  FLA_Obj V, FLA_Obj D )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Apply_QUD_UT_inc_check( side, trans, direct, storev, T, W, R, U, C, V, D );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Invoke _internal() back-end with the standard control tree.
  r_val = FLA_Apply_QUD_UT_inc_internal( side, trans, direct, storev,
                                         T, W, R, U, C, V, D, flash_apqudutinc_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

