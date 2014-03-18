
#include "FLAME.h"

extern fla_appiv_t* fla_appiv_cntl_leaf;

FLA_Error FLA_Apply_pivots_task( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
   //return FLA_Apply_pivots_unb_external( side, trans, p, A );
   return FLA_Apply_pivots_internal( side, trans,
                                     p, A,
                                     fla_appiv_cntl_leaf );
}

FLA_Error FLA_Apply_pivots_ln_task( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
   //return FLA_Apply_pivots_unb_external( FLA_LEFT, FLA_NO_TRANSPOSE, p, A );
   return FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE,
                                     p, A,
                                     fla_appiv_cntl_leaf );
}

