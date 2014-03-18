
#include "FLAME.h"

FLA_Error REF_Lyap_h( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale )
{
  FLA_Trans transA;

  if ( FLA_Obj_is_complex( A ) ) transA = FLA_CONJ_TRANSPOSE;
  else                           transA = FLA_TRANSPOSE;

  if ( FLA_Obj_is( isgn, FLA_MINUS_ONE ) )
    FLA_Scal( FLA_MINUS_ONE, C );

  FLA_Sylv_blk_external( transA, FLA_NO_TRANSPOSE, FLA_ONE, A, A, C, scale );
  //FLA_Sylv( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, A, C, scale );

  return FLA_SUCCESS;
}
