/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_trinv_t* fla_trinv_cntl_leaf;

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

