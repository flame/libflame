/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Bsvd_n.h"   // No U and V are requested.
#include "FLA_Bsvd_v.h"   // Vectors of U and V are requested.
#include "FLA_Bsvd_ext.h" // Extension to LAPACK-like interface.

FLA_Error FLA_Bsvd_create_workspace( FLA_Obj d, FLA_Obj *G, FLA_Obj *H );
FLA_Error FLA_Bsvd( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H, 
                    FLA_Svd_type jobu, FLA_Obj U, 
                    FLA_Svd_type jobv, FLA_Obj V );
FLA_Error FLA_Bsvd_ext( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H,
                        FLA_Svd_type jobu, FLA_Obj U,
                        FLA_Svd_type jobv, FLA_Obj V,
                        FLA_Bool apply_Uh2C, FLA_Obj C );
