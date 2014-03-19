/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Apply_QUD_UT_inc_lhfc.h"

FLA_Error FLASH_Apply_QUD_UT_inc( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                  FLA_Obj T, FLA_Obj W,
                                             FLA_Obj R,
                                  FLA_Obj U, FLA_Obj C,
                                  FLA_Obj V, FLA_Obj D );

FLA_Error FLA_Apply_QUD_UT_inc_internal( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                         FLA_Obj T, FLA_Obj W,
                                                    FLA_Obj R,
                                         FLA_Obj U, FLA_Obj C,
                                         FLA_Obj V, FLA_Obj D, fla_apqudutinc_t* cntl );

FLA_Error FLA_Apply_QUD_UT_inc_lhfc( FLA_Obj T, FLA_Obj W,
                                                FLA_Obj R,
                                     FLA_Obj U, FLA_Obj C,
                                     FLA_Obj V, FLA_Obj D, fla_apqudutinc_t* cntl );

FLA_Error FLASH_Apply_QUD_UT_inc_create_workspace( FLA_Obj T, FLA_Obj R, FLA_Obj* W );

