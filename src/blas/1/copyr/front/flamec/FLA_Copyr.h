
#include "FLA_Copyr_l.h"
#include "FLA_Copyr_u.h"

FLA_Error FLASH_Copyr( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B );

FLA_Error FLA_Copyr_internal( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl );
FLA_Error FLA_Copyr_l( FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl );
FLA_Error FLA_Copyr_u( FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl );

