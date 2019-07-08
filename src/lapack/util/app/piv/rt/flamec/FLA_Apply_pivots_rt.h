/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Apply_pivots_rt_opt_var1_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj p, FLA_Obj A );
#endif
FLA_Error FLA_Apply_pivots_rt_opt_var1( FLA_Obj p, FLA_Obj A );

