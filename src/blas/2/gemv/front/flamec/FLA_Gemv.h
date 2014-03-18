
#include "FLA_Gemv_h.h"
#include "FLA_Gemv_n.h"
#include "FLA_Gemv_t.h"

FLA_Error FLA_Gemv_internal( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );

FLA_Error FLA_Gemv_h( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );
FLA_Error FLA_Gemv_n( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );
FLA_Error FLA_Gemv_t( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );

