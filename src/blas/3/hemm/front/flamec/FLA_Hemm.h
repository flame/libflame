
#include "FLA_Hemm_ll.h"
#include "FLA_Hemm_lu.h"
#include "FLA_Hemm_rl.h"
#include "FLA_Hemm_ru.h"

FLA_Error FLA_Hemm_internal( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_hemm_t* cntl );

FLA_Error FLA_Hemm_ll( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_hemm_t* cntl );
FLA_Error FLA_Hemm_lu( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_hemm_t* cntl );
FLA_Error FLA_Hemm_rl( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_hemm_t* cntl );
FLA_Error FLA_Hemm_ru( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_hemm_t* cntl );

