
#include "FLA_Scalr_l.h"
#include "FLA_Scalr_u.h"

FLA_Error FLA_Scalr_internal( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, fla_scalr_t* cntl );

FLA_Error FLA_Scalr_l( FLA_Obj alpha, FLA_Obj A, fla_scalr_t* cntl );
FLA_Error FLA_Scalr_u( FLA_Obj alpha, FLA_Obj A, fla_scalr_t* cntl );

