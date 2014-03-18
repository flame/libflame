
// --- top-level front-end prototypes ------------------------------------------

FLA_Error FLASH_Axpy( FLA_Obj alpha, FLA_Obj A, FLA_Obj B );
FLA_Error FLASH_Axpyt( FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B );
FLA_Error FLASH_Copy( FLA_Obj A, FLA_Obj B );
FLA_Error FLASH_Copyt( FLA_Trans trans, FLA_Obj A, FLA_Obj B );
FLA_Error FLASH_Scal( FLA_Obj alpha, FLA_Obj A );
FLA_Error FLASH_Scalr( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A );

