/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Svd_ext.h"
#include "FLA_Svd_uv.h"

FLA_Error FLA_Svd_compute_scaling( FLA_Obj A, FLA_Obj sigma );

FLA_Error FLA_Svd( FLA_Svd_type jobu, FLA_Svd_type jobv, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V );
FLA_Error FLA_Svd_ext( FLA_Svd_type jobu, FLA_Trans transu,
                       FLA_Svd_type jobv, FLA_Trans transv,
                       FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V );
