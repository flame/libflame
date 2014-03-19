/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error REF_Lyap_n( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale )
{
  FLA_Trans transB;

  if ( FLA_Obj_is( isgn, FLA_MINUS_ONE ) )
    FLA_Scal( FLA_MINUS_ONE, C );

  if ( FLA_Obj_is_complex( A ) ) transB = FLA_CONJ_TRANSPOSE;
  else                           transB = FLA_TRANSPOSE;

  FLA_Sylv_blk_external( FLA_NO_TRANSPOSE, transB, FLA_ONE, A, A, C, scale );
  //FLA_Sylv( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, FLA_ONE, A, A, C, scale );

  return FLA_SUCCESS;
}
