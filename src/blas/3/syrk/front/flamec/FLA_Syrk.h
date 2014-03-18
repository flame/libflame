
#include "FLA_Syrk_ln.h"
#include "FLA_Syrk_lt.h"
#include "FLA_Syrk_un.h"
#include "FLA_Syrk_ut.h"

FLA_Error FLA_Syrk_internal( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl );

FLA_Error FLA_Syrk_ln( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl );
FLA_Error FLA_Syrk_lt( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl );
FLA_Error FLA_Syrk_un( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl );
FLA_Error FLA_Syrk_ut( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl );

