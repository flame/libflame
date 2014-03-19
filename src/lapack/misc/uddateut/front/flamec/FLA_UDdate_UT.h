/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_UDdate_UT_vars.h"

FLA_Error FLA_UDdate_UT( FLA_Obj R,
                         FLA_Obj C,
                         FLA_Obj D, FLA_Obj T );

FLA_Error FLA_UDdate_UT_internal( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl );

FLA_Error FLA_UDdate_UT_create_T( FLA_Obj R, FLA_Obj* T );

FLA_Error FLA_UDdate_UT_update_rhs( FLA_Obj T, FLA_Obj bR,
                                    FLA_Obj C, FLA_Obj bC,
                                    FLA_Obj D, FLA_Obj bD );

FLA_Error FLA_UDdate_UT_solve( FLA_Obj R, FLA_Obj bR, FLA_Obj x );
