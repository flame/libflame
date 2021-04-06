/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Hess( FLA_Obj A, FLA_Obj t, integer ilo, integer ihi );
FLA_Error FLA_Hess_internal( FLA_Obj A, FLA_Obj t, integer ilo, integer ihi, fla_hess_t* cntl );
FLA_Error FLA_Hess_blk_var1( FLA_Obj A, FLA_Obj t, integer ilo, integer ihi, fla_hess_t* cntl );
FLA_Error FLA_Hess_blk_var2( FLA_Obj A, FLA_Obj t, integer ilo, integer ihi, fla_hess_t* cntl );
FLA_Error FLA_Hess_step_unb_var1( FLA_Obj A, FLA_Obj t, FLA_Obj B, FLA_Obj C, integer nb_alg );
FLA_Error FLA_Hess_step_unb_var2( FLA_Obj A, FLA_Obj t, FLA_Obj B, FLA_Obj C, integer nb_alg );

void MyFLA_extract_column_block( FLA_Obj A, integer ini, integer n, FLA_Obj * B );
void MyFLA_extract_row_block( FLA_Obj A, integer ini, integer n, FLA_Obj * B );
void MyFLA_int_swap( integer * a, integer * b );
void MyFLA_v_int_swap( integer n, integer * a, integer * b );

