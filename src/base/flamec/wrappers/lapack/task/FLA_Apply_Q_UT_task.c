/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_apqut_t* fla_apqut_cntl_leaf;

FLA_Error FLA_Apply_Q_UT_task( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( side, trans, direct, storev,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_lhbc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_BACKWARD, FLA_COLUMNWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_lhbr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_BACKWARD, FLA_ROWWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_lhfc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_lhfr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_lnbc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_BACKWARD, FLA_COLUMNWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_lnbr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_BACKWARD, FLA_ROWWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_lnfc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_lnfr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_rhbc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_RIGHT, FLA_CONJ_TRANSPOSE, FLA_BACKWARD, FLA_COLUMNWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_rhbr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_RIGHT, FLA_CONJ_TRANSPOSE, FLA_BACKWARD, FLA_ROWWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_rhfc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_RIGHT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_rhfr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_RIGHT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_rnbc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_BACKWARD, FLA_COLUMNWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_rnbr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_BACKWARD, FLA_ROWWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_rnfc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

FLA_Error FLA_Apply_Q_UT_rnfr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  return FLA_Apply_Q_UT_internal( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                                  A, T, W, B,
                                  fla_apqut_cntl_leaf );
}

