
#include "FLA_Symm_ll.h"
#include "FLA_Symm_lu.h"
#include "FLA_Symm_rl.h"
#include "FLA_Symm_ru.h"

FLA_Error FLA_Symm_internal( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl );

FLA_Error FLA_Symm_ll( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl );
FLA_Error FLA_Symm_lu( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl );
FLA_Error FLA_Symm_rl( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl );
FLA_Error FLA_Symm_ru( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl );

