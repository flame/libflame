
#include "FLAME.h"

extern fla_appiv_t* fla_appiv_cntl_leaf;

FLA_Error FLA_Apply_pivots( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  
  r_val = FLA_Apply_pivots_internal( side, trans, p, A, fla_appiv_cntl_leaf );
  
  return r_val;
}

