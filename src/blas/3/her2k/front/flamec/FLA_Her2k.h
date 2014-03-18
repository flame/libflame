
#include "FLA_Her2k_lh.h"
#include "FLA_Her2k_ln.h"
#include "FLA_Her2k_uh.h"
#include "FLA_Her2k_un.h"

FLA_Error FLA_Her2k_internal( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl );

FLA_Error FLA_Her2k_lh( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl );
FLA_Error FLA_Her2k_ln( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl );
FLA_Error FLA_Her2k_uh( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl );
FLA_Error FLA_Her2k_un( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl );

