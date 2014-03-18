
#include "FLA_Trsv_lc.h"
#include "FLA_Trsv_ln.h"
#include "FLA_Trsv_lt.h"
#include "FLA_Trsv_uc.h"
#include "FLA_Trsv_un.h"
#include "FLA_Trsv_ut.h"

FLA_Error FLA_Trsv_internal( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );

FLA_Error FLA_Trsv_lc( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_ln( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_lt( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_uc( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_un( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_ut( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );

