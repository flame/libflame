/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Copyr_l.h"
#include "FLA_Copyr_u.h"

FLA_Error FLASH_Copyr( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B );

FLA_Error FLA_Copyr_internal( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl );
FLA_Error FLA_Copyr_l( FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl );
FLA_Error FLA_Copyr_u( FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl );

