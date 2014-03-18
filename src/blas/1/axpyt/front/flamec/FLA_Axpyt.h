
#include "FLA_Axpyt_n.h"
#include "FLA_Axpyt_t.h"
#include "FLA_Axpyt_c.h"
#include "FLA_Axpyt_h.h"

FLA_Error FLA_Axpyt_internal( FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl );
FLA_Error FLA_Axpyt_n( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl );
FLA_Error FLA_Axpyt_t( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl );
FLA_Error FLA_Axpyt_c( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl );
FLA_Error FLA_Axpyt_h( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl );

