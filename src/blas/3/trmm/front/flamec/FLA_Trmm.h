/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Trmm_llc.h"
#include "FLA_Trmm_llh.h"
#include "FLA_Trmm_lln.h"
#include "FLA_Trmm_llt.h"
#include "FLA_Trmm_luc.h"
#include "FLA_Trmm_luh.h"
#include "FLA_Trmm_lun.h"
#include "FLA_Trmm_lut.h"
#include "FLA_Trmm_rlc.h"
#include "FLA_Trmm_rlh.h"
#include "FLA_Trmm_rln.h"
#include "FLA_Trmm_rlt.h"
#include "FLA_Trmm_ruc.h"
#include "FLA_Trmm_ruh.h"
#include "FLA_Trmm_run.h"
#include "FLA_Trmm_rut.h"

FLA_Error FLA_Trmm_internal( FLA_Side side, FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );

FLA_Error FLA_Trmm_llc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_llh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_lln( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_llt( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_luc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_luh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_lun( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_lut( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_rlc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_rlh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_rln( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_rlt( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_ruc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_ruh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_run( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );
FLA_Error FLA_Trmm_rut( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl );

