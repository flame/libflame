
#include "FLA_Trinv_ln.h"
#include "FLA_Trinv_lu.h"
#include "FLA_Trinv_un.h"
#include "FLA_Trinv_uu.h"

FLA_Error FLA_Trinv_internal( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A, fla_trinv_t* cntl );

FLA_Error FLA_Trinv_ln( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_lu( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_un( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_uu( FLA_Obj A, fla_trinv_t* cntl );

