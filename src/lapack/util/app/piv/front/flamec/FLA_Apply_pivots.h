/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Apply_pivots_ln.h"
#include "FLA_Apply_pivots_lt.h"
#include "FLA_Apply_pivots_rn.h"
#include "FLA_Apply_pivots_rt.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Apply_pivots_internal_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
#endif
FLA_Error FLA_Apply_pivots_internal( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Apply_pivots_ln_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
#endif
FLA_Error FLA_Apply_pivots_ln( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Apply_pivots_lt_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
#endif
FLA_Error FLA_Apply_pivots_lt( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Apply_pivots_rn_ts(FLA_cntl_init_s *FLA_cntl_init_i,  FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
#endif
FLA_Error FLA_Apply_pivots_rn( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Apply_pivots_rt_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
#endif
FLA_Error FLA_Apply_pivots_rt( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );

