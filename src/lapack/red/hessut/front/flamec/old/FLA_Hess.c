
#include "FLAME.h"

extern fla_hess_t*      fla_hess_cntl;
extern fla_blocksize_t* fla_hess_blocksize;

FLA_Error FLA_Hess( FLA_Obj A, FLA_Obj t, int ilo, int ihi )
{
  FLA_Error r_val;

  // Note the following constraints:  0 <= ilo <= n-1,  0 <= ihi <= n-1.
  r_val = FLA_Hess_internal( A, t, ilo, ihi, fla_hess_cntl );

  return r_val;
}

