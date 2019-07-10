/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Gemm_cc.h"
#include "FLA_Gemm_ch.h"
#include "FLA_Gemm_cn.h"
#include "FLA_Gemm_ct.h"
#include "FLA_Gemm_hc.h"
#include "FLA_Gemm_hh.h"
#include "FLA_Gemm_hn.h"
#include "FLA_Gemm_ht.h"
#include "FLA_Gemm_nc.h"
#include "FLA_Gemm_nh.h"
#include "FLA_Gemm_nn.h"
#include "FLA_Gemm_nt.h"
#include "FLA_Gemm_tc.h"
#include "FLA_Gemm_th.h"
#include "FLA_Gemm_tn.h"
#include "FLA_Gemm_tt.h"

FLA_Error FLA_Gemm_internal( FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );

FLA_Error FLA_Gemm_cc( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_ch( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_cn( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_ct( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_hc( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_hh( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_hn( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_ht( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_nc( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_nh( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_nn( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_nt( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_tc( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_th( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_tn( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );
FLA_Error FLA_Gemm_tt( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl );

