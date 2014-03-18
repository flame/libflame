
#include "FLAME.h"

FLA_Error REF_Eig_gest( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_blk_external( inv, uplo, A, B );
}

