
#include "FLA_Apply_pivots_ln.h"
#include "FLA_Apply_pivots_lt.h"
#include "FLA_Apply_pivots_rn.h"
#include "FLA_Apply_pivots_rt.h"

FLA_Error FLA_Apply_pivots_internal( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );

FLA_Error FLA_Apply_pivots_ln( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
FLA_Error FLA_Apply_pivots_lt( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
FLA_Error FLA_Apply_pivots_rn( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
FLA_Error FLA_Apply_pivots_rt( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );

