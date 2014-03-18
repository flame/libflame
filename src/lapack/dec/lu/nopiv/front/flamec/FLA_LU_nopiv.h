
#include "FLA_LU_nopiv_vars.h"

FLA_Error FLA_LU_nopiv_internal( FLA_Obj A, fla_lu_t* cntl );

FLA_Error FLA_LU_nopiv_solve( FLA_Obj A, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_LU_nopiv_solve( FLA_Obj A, FLA_Obj B, FLA_Obj X );
