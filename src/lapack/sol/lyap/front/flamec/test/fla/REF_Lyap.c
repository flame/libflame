
#include "FLAME.h"

FLA_Error REF_Lyap( FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale )
{
  if ( FLA_Obj_is( isgn, FLA_MINUS_ONE ) )
    FLA_Scal( FLA_MINUS_ONE, C );

  if ( FLA_Obj_is_complex( A ) )
  {
    if ( trans == FLA_NO_TRANSPOSE )
      FLA_Sylv_blk_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, FLA_ONE, A, A, C, scale );
    else
      FLA_Sylv_blk_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, A, C, scale );
  }
  else // if ( FLA_Obj_is_real( A ) )
  {
    if ( trans == FLA_NO_TRANSPOSE )
      FLA_Sylv_blk_external( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A, A, C, scale );
    else
      FLA_Sylv_blk_external( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, A, C, scale );
  }

  return FLA_SUCCESS;
}

