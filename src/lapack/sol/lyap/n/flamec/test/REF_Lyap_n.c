
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
