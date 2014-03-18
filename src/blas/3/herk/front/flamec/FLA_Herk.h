
#include "FLA_Herk_lh.h"
#include "FLA_Herk_ln.h"
#include "FLA_Herk_uh.h"
#include "FLA_Herk_un.h"

FLA_Error FLA_Herk_internal( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_herk_t* cntl );

FLA_Error FLA_Herk_lh( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_herk_t* cntl );
FLA_Error FLA_Herk_ln( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_herk_t* cntl );
FLA_Error FLA_Herk_uh( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_herk_t* cntl );
FLA_Error FLA_Herk_un( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_herk_t* cntl );

