/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- top-level front-end prototypes ------------------------------------------

FLA_Error FLASH_Axpy( FLA_Obj alpha, FLA_Obj A, FLA_Obj B );
FLA_Error FLASH_Axpyt( FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B );
FLA_Error FLASH_Copy( FLA_Obj A, FLA_Obj B );
FLA_Error FLASH_Copyt( FLA_Trans trans, FLA_Obj A, FLA_Obj B );
FLA_Error FLASH_Scal( FLA_Obj alpha, FLA_Obj A );
FLA_Error FLASH_Scalr( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A );

