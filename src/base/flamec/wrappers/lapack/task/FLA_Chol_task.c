/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_chol_t* fla_chol_cntl_leaf;

FLA_Error FLA_Chol_task( FLA_Uplo uplo, FLA_Obj A, fla_chol_t* cntl )
{
  return FLA_Chol_internal( uplo, A,
                            fla_chol_cntl_leaf );
}

FLA_Error FLA_Chol_l_task( FLA_Obj A, fla_chol_t* cntl )
{
  //return FLA_Chol_unb_external( FLA_LOWER_TRIANGULAR, A );
  return FLA_Chol_internal( FLA_LOWER_TRIANGULAR, A,
                            fla_chol_cntl_leaf );
}

FLA_Error FLA_Chol_u_task( FLA_Obj A, fla_chol_t* cntl )
{
  //return FLA_Chol_unb_external( FLA_UPPER_TRIANGULAR, A );
  return FLA_Chol_internal( FLA_UPPER_TRIANGULAR, A,
                            fla_chol_cntl_leaf );
}

