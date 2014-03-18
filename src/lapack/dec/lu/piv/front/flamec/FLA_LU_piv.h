
#include "FLA_LU_piv_vars.h"

FLA_Error FLA_LU_piv_internal( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl );

FLA_Error FLA_LU_piv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_LU_piv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj B, FLA_Obj X );
