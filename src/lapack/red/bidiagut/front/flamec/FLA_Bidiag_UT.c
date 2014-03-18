
#include "FLAME.h"

extern fla_bidiagut_t* fla_bidiagut_cntl_fused;
extern fla_bidiagut_t* fla_bidiagut_cntl_nofus;
extern fla_bidiagut_t* fla_bidiagut_cntl_plain;

FLA_Error FLA_Bidiag_UT( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Bidiag_UT_check( A, TU, TV );

  if ( FLA_Obj_row_stride( A ) == 1 &&
       FLA_Obj_row_stride( TU ) == 1 &&
       FLA_Obj_row_stride( TV ) == 1 &&
       FLA_Obj_is_double_precision( A ) )
    // Temporary modification to "nofus"; 
    // fused operations are not working for row-major, ex) bl1_ddotsv2 
    r_val = FLA_Bidiag_UT_internal( A, TU, TV, fla_bidiagut_cntl_plain );
  else
    r_val = FLA_Bidiag_UT_internal( A, TU, TV, fla_bidiagut_cntl_plain );

  return r_val;
}

