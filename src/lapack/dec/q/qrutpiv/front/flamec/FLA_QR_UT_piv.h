
#include "FLA_QR_UT_piv_vars.h"

FLA_Error FLA_QR_UT_piv( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p );

FLA_Error FLA_QR_UT_piv_internal( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_piv_colnorm( FLA_Obj alpha, FLA_Obj A, FLA_Obj b );

// The source files are located at src/base/flamec/check/lapack
FLA_Error FLA_QR_UT_piv_check( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p );
FLA_Error FLA_QR_UT_piv_internal_check( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_piv_colnorm_check( FLA_Obj alpha, FLA_Obj A, FLA_Obj b );
