/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Hess_UT_vars.h"

FLA_Error FLA_Hess_UT( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Hess_UT_internal( FLA_Obj A, FLA_Obj T, fla_hessut_t* cntl );

FLA_Error FLA_Hess_UT_create_T( FLA_Obj A, FLA_Obj* T );

FLA_Error FLA_Hess_UT_recover_tau( FLA_Obj T, FLA_Obj t );
