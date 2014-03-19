/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Apply_G_lf.h"
#include "FLA_Apply_G_lb.h"
#include "FLA_Apply_G_rf.h"
#include "FLA_Apply_G_rb.h"

FLA_Error FLA_Apply_G( FLA_Side side, FLA_Direct direct, FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_internal( FLA_Side side, FLA_Direct direct, FLA_Obj G, FLA_Obj A );

#include "FLA_Givens2.h"

#include "FLA_Apply_GTG.h"

#include "FLA_Apply_GT_2x2.h"
#include "FLA_Apply_G_2x2.h"
#include "FLA_Apply_G_1x2.h"

#include "FLA_Apply_G_mx2_opt.h"
#include "FLA_Apply_G_mx2_asm.h"

#include "FLA_Apply_G_mx3_opt.h"
#include "FLA_Apply_G_mx3b_opt.h"
#include "FLA_Apply_G_mx3_asm.h"
#include "FLA_Apply_G_mx3b_asm.h"

#include "FLA_Apply_G_mx4s_opt.h"
#include "FLA_Apply_G_mx4s_asm.h"

