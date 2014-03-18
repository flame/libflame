
#include "FLA_Ttmm_l.h"
#include "FLA_Ttmm_u.h"

FLA_Error FLA_Ttmm_internal( FLA_Uplo uplo, FLA_Obj A, fla_ttmm_t* cntl );
FLA_Error FLA_Ttmm_l( FLA_Obj A, fla_ttmm_t* cntl );
FLA_Error FLA_Ttmm_u( FLA_Obj A, fla_ttmm_t* cntl );
