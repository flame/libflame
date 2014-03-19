/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Copyt_task( FLA_Trans trans, FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl )
{
  return FLA_Copyt_external( trans, A, B );
}

FLA_Error FLA_Copyt_n_task( FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl )
{
  return FLA_Copyt_external( FLA_NO_TRANSPOSE, A, B );
}

FLA_Error FLA_Copyt_t_task( FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl )
{
  return FLA_Copyt_external( FLA_TRANSPOSE, A, B );
}

FLA_Error FLA_Copyt_c_task( FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl )
{
  return FLA_Copyt_external( FLA_CONJ_NO_TRANSPOSE, A, B );
}

FLA_Error FLA_Copyt_h_task( FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl )
{
  return FLA_Copyt_external( FLA_CONJ_TRANSPOSE, A, B );
}
