/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Trinv_ln.h"
#include "FLA_Trinv_lu.h"
#include "FLA_Trinv_un.h"
#include "FLA_Trinv_uu.h"

FLA_Error FLA_Trinv_internal( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A, fla_trinv_t* cntl );

FLA_Error FLA_Trinv_ln( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_lu( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_un( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_uu( FLA_Obj A, fla_trinv_t* cntl );

