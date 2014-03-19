/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Hevd_ln.h"
#include "FLA_Hevd_lv.h"
//#include "FLA_Hevd_un.h"
//#include "FLA_Hevd_uv.h"

FLA_Error FLA_Hevd_compute_scaling( FLA_Uplo uplo, FLA_Obj A, FLA_Obj sigma );

FLA_Error FLA_Hevd( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l );

