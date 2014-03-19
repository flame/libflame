/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Copyt_n.h"
#include "FLA_Copyt_t.h"
#include "FLA_Copyt_c.h"
#include "FLA_Copyt_h.h"

FLA_Error FLA_Copyt_internal( FLA_Trans trans, FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl );
FLA_Error FLA_Copyt_n( FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl );
FLA_Error FLA_Copyt_t( FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl );
FLA_Error FLA_Copyt_c( FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl );
FLA_Error FLA_Copyt_h( FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl );

