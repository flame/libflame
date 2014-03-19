/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Ttmm_l.h"
#include "FLA_Ttmm_u.h"

FLA_Error FLA_Ttmm_internal( FLA_Uplo uplo, FLA_Obj A, fla_ttmm_t* cntl );
FLA_Error FLA_Ttmm_l( FLA_Obj A, fla_ttmm_t* cntl );
FLA_Error FLA_Ttmm_u( FLA_Obj A, fla_ttmm_t* cntl );
