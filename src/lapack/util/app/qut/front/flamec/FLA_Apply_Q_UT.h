/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Apply_Q_UT_lnfc.h"
#include "FLA_Apply_Q_UT_lnfr.h"
#include "FLA_Apply_Q_UT_lnbc.h"
#include "FLA_Apply_Q_UT_lnbr.h"
#include "FLA_Apply_Q_UT_lhfc.h"
#include "FLA_Apply_Q_UT_lhfr.h"
#include "FLA_Apply_Q_UT_lhbc.h"
#include "FLA_Apply_Q_UT_lhbr.h"

#include "FLA_Apply_Q_UT_rhbc.h"
#include "FLA_Apply_Q_UT_rhbr.h"
#include "FLA_Apply_Q_UT_rhfc.h"
#include "FLA_Apply_Q_UT_rhfr.h"
#include "FLA_Apply_Q_UT_rnbc.h"
#include "FLA_Apply_Q_UT_rnbr.h"
#include "FLA_Apply_Q_UT_rnfc.h"
#include "FLA_Apply_Q_UT_rnfr.h"

FLA_Error FLA_Apply_Q_UT_internal( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );

FLA_Error FLA_Apply_Q_UT_lnfc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lnfr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lnbc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lnbr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhfc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhfr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhbc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhbr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );

FLA_Error FLA_Apply_Q_UT_rhbc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rhbr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rhfc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rhfr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnbc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnbr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnfc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnfr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );

FLA_Error FLA_Apply_Q_UT_create_workspace( FLA_Obj T, FLA_Obj B, FLA_Obj* W );
FLA_Error FLA_Apply_Q_UT_create_workspace_side( FLA_Side side, FLA_Obj T, FLA_Obj B, FLA_Obj* W );

FLA_Error FLASH_Apply_Q_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B );
FLA_Error FLASH_Apply_Q_UT_create_workspace( FLA_Obj TW, FLA_Obj B, FLA_Obj* W );

