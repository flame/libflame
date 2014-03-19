/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Sylv_blk_external( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  FLA_Error e_val = FLA_FAILURE;

#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  e_val = FLA_Sylv_unb_external( transa, transb, isgn, A, B, C, scale );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return e_val;
}

FLA_Error FLA_Sylv_nn_blk_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  return FLA_Sylv_blk_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
}

FLA_Error FLA_Sylv_nh_blk_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  return FLA_Sylv_blk_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, isgn, A, B, C, scale );
}

FLA_Error FLA_Sylv_hn_blk_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  return FLA_Sylv_blk_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
}

FLA_Error FLA_Sylv_hh_blk_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  return FLA_Sylv_blk_external( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE, isgn, A, B, C, scale );
}

