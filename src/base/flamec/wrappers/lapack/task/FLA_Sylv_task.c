
#include "FLAME.h"

extern fla_sylv_t* fla_sylv_cntl_leaf;

FLA_Error FLA_Sylv_task( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  return FLA_Sylv_internal( transa, transb,
                            isgn, A, B, C, scale,
                            fla_sylv_cntl_leaf );
}

FLA_Error FLA_Sylv_nn_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  //return FLA_Sylv_unb_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
  return FLA_Sylv_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                            isgn, A, B, C, scale,
                            fla_sylv_cntl_leaf );
}

FLA_Error FLA_Sylv_nh_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  //return FLA_Sylv_unb_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, isgn, A, B, C, scale );
  return FLA_Sylv_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                            isgn, A, B, C, scale,
                            fla_sylv_cntl_leaf );
}

FLA_Error FLA_Sylv_hn_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  //return FLA_Sylv_unb_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
  return FLA_Sylv_internal( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                            isgn, A, B, C, scale,
                            fla_sylv_cntl_leaf );
}

FLA_Error FLA_Sylv_hh_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  //return FLA_Sylv_unb_external( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE, isgn, A, B, C, scale );
  return FLA_Sylv_internal( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                            isgn, A, B, C, scale,
                            fla_sylv_cntl_leaf );
}

