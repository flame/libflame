
#include "FLAME.h"

extern fla_trinv_t* fla_trinv_cntl_leaf;

FLA_Error FLA_Trinv_task( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A, fla_trinv_t* cntl )
{
  return FLA_Trinv_internal( uplo, diag, A,
                             fla_trinv_cntl_leaf );
}

FLA_Error FLA_Trinv_ln_task( FLA_Obj A, fla_trinv_t* cntl )
{
  //return FLA_Trinv_unb_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
  return FLA_Trinv_internal( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A,
                             fla_trinv_cntl_leaf );
}

FLA_Error FLA_Trinv_lu_task( FLA_Obj A, fla_trinv_t* cntl )
{
  //return FLA_Trinv_unb_external( FLA_LOWER_TRIANGULAR, FLA_UNIT_DIAG, A );
  return FLA_Trinv_internal( FLA_LOWER_TRIANGULAR, FLA_UNIT_DIAG, A,
                             fla_trinv_cntl_leaf );
}

FLA_Error FLA_Trinv_un_task( FLA_Obj A, fla_trinv_t* cntl )
{
  //return FLA_Trinv_unb_external( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
  return FLA_Trinv_internal( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A,
                             fla_trinv_cntl_leaf );
}

FLA_Error FLA_Trinv_uu_task( FLA_Obj A, fla_trinv_t* cntl )
{
  //return FLA_Trinv_unb_external( FLA_UPPER_TRIANGULAR, FLA_UNIT_DIAG, A );
  return FLA_Trinv_internal( FLA_UPPER_TRIANGULAR, FLA_UNIT_DIAG, A,
                             fla_trinv_cntl_leaf );
}

