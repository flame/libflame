
#include "FLA_Apply_CAQ_UT_inc_lhfc.h"

FLA_Error FLASH_Apply_CAQ_UT_inc( dim_t p,
                                  FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                  FLA_Obj A, FLA_Obj ATW,
                                  FLA_Obj R, FLA_Obj RTW, FLA_Obj W, FLA_Obj B );

FLA_Error FLA_Apply_CAQ_UT_inc_apply_panels( dim_t nb_part, FLA_Obj A, FLA_Obj ATW, FLA_Obj W, FLA_Obj B );

FLA_Error FLASH_Apply_CAQ_UT_inc_create_workspace( dim_t p, FLA_Obj TW, FLA_Obj B, FLA_Obj* W );

FLA_Error FLA_Apply_CAQ_UT_inc_internal( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                         FLA_Obj R, FLA_Obj TW, FLA_Obj W1, FLA_Obj B,
                                         fla_apcaqutinc_t* cntl );

FLA_Error FLA_Apply_CAQ_UT_inc_lhfc( FLA_Obj R, FLA_Obj TW, FLA_Obj W1, FLA_Obj B,
                                     fla_apcaqutinc_t* cntl );

