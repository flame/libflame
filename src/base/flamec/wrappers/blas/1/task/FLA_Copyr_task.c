/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Copyr_task( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl )
{
  return FLA_Copyr_external( uplo, A, B );
}

FLA_Error FLA_Copyr_l_task( FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl )
{
  return FLA_Copyr_external( FLA_LOWER_TRIANGULAR, A, B );
}

FLA_Error FLA_Copyr_u_task( FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl )
{
  return FLA_Copyr_external( FLA_UPPER_TRIANGULAR, A, B );
}

