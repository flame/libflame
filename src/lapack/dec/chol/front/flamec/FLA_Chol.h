/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Chol_l.h"
#include "FLA_Chol_u.h"

FLA_Error FLA_Chol_internal( FLA_Uplo uplo, FLA_Obj A, fla_chol_t* cntl );
FLA_Error FLA_Chol_l( FLA_Obj A, fla_chol_t* cntl );
FLA_Error FLA_Chol_u( FLA_Obj A, fla_chol_t* cntl );

FLA_Error FLA_Chol_solve( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_Chol_solve( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj X );
