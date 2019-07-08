/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_LU_piv_task_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj A, FLA_Obj p, fla_lu_t* cntl )
{
  return FLA_LU_piv_internal_ts( FLA_cntl_init_i, A, p,
                              FLA_cntl_init_i->FLA_Cntl_init_flamec_i->fla_lu_piv_cntl_leaf );
}
#endif

extern fla_lu_t* fla_lu_piv_cntl_leaf;

FLA_Error FLA_LU_piv_task( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl )
{
  return FLA_LU_piv_internal( A, p,
                              fla_lu_piv_cntl_leaf );
}

