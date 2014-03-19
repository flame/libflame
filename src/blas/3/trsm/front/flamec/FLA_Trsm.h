/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Trsm_llc.h"
#include "FLA_Trsm_llh.h"
#include "FLA_Trsm_lln.h"
#include "FLA_Trsm_llt.h"
#include "FLA_Trsm_luc.h"
#include "FLA_Trsm_luh.h"
#include "FLA_Trsm_lun.h"
#include "FLA_Trsm_lut.h"
#include "FLA_Trsm_rlc.h"
#include "FLA_Trsm_rlh.h"
#include "FLA_Trsm_rln.h"
#include "FLA_Trsm_rlt.h"
#include "FLA_Trsm_ruc.h"
#include "FLA_Trsm_ruh.h"
#include "FLA_Trsm_run.h"
#include "FLA_Trsm_rut.h"

FLA_Error FLA_Trsm_internal( FLA_Side side, FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );

FLA_Error FLA_Trsm_llc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_llh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_lln( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_llt( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_luc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_luh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_lun( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_lut( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_rlc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_rlh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_rln( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_rlt( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_ruc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_ruh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_run( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_rut( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );

