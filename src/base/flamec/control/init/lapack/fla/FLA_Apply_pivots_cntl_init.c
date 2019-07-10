/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

fla_appiv_t* fla_appiv_cntl_leaf = NULL;

void FLA_Apply_pivots_cntl_init()
{
	// Create a control tree to invoke LAPACK.
	fla_appiv_cntl_leaf = FLA_Cntl_appiv_obj_create( FLA_FLAT,
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS
	                                                 FLA_UNBLOCKED_EXTERN,
#else
	                                                 FLA_UNB_OPT_VARIANT1,
#endif
	                                                 NULL,
	                                                 NULL );
}

void FLA_Apply_pivots_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_appiv_cntl_leaf );
}

