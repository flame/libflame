
#include "FLA_Apply_QUD_UT_lhfc.h"

FLA_Error FLA_Apply_QUD_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                            FLA_Obj T, FLA_Obj W,
                                       FLA_Obj R,
                            FLA_Obj U, FLA_Obj C,
                            FLA_Obj V, FLA_Obj D );

FLA_Error FLA_Apply_QUD_UT_internal( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                     FLA_Obj T, FLA_Obj W,
                                                FLA_Obj R,
                                     FLA_Obj U, FLA_Obj C,
                                     FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl );

FLA_Error FLA_Apply_QUD_UT_lhfc( FLA_Obj T, FLA_Obj W,
                                            FLA_Obj R,
                                 FLA_Obj U, FLA_Obj C,
                                 FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl );

FLA_Error FLA_Apply_QUD_UT_create_workspace( FLA_Obj T, FLA_Obj R, FLA_Obj* W );

