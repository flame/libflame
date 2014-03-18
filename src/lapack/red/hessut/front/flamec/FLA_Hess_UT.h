
#include "FLA_Hess_UT_vars.h"

FLA_Error FLA_Hess_UT( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Hess_UT_internal( FLA_Obj A, FLA_Obj T, fla_hessut_t* cntl );

FLA_Error FLA_Hess_UT_create_T( FLA_Obj A, FLA_Obj* T );

FLA_Error FLA_Hess_UT_recover_tau( FLA_Obj T, FLA_Obj t );
