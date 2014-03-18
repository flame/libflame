
#include "FLA_Chol_l.h"
#include "FLA_Chol_u.h"

FLA_Error FLA_Chol_internal( FLA_Uplo uplo, FLA_Obj A, fla_chol_t* cntl );
FLA_Error FLA_Chol_l( FLA_Obj A, fla_chol_t* cntl );
FLA_Error FLA_Chol_u( FLA_Obj A, fla_chol_t* cntl );

FLA_Error FLA_Chol_solve( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_Chol_solve( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj X );
