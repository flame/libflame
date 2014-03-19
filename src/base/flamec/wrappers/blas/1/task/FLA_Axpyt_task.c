/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Axpyt_task( FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl )
{
  return FLA_Axpyt_external( trans, alpha, A, B );
}

FLA_Error FLA_Axpyt_n_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl )
{
  return FLA_Axpyt_external( FLA_NO_TRANSPOSE, alpha, A, B );
}

FLA_Error FLA_Axpyt_t_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl )
{
  return FLA_Axpyt_external( FLA_TRANSPOSE, alpha, A, B );
}

FLA_Error FLA_Axpyt_c_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl )
{
  return FLA_Axpyt_external( FLA_CONJ_NO_TRANSPOSE, alpha, A, B );
}

FLA_Error FLA_Axpyt_h_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl )
{
  return FLA_Axpyt_external( FLA_CONJ_TRANSPOSE, alpha, A, B );
}
