/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

//#include "FLA_Bidiag_UT_l.h"
#include "FLA_Bidiag_UT_u.h"

FLA_Error FLA_Bidiag_UT( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_internal( FLA_Obj A, FLA_Obj TU, FLA_Obj TV, fla_bidiagut_t* cntl );

FLA_Error FLA_Bidiag_UT_l( FLA_Obj A, FLA_Obj TU, FLA_Obj TV, fla_bidiagut_t* cntl );
FLA_Error FLA_Bidiag_UT_u( FLA_Obj A, FLA_Obj TU, FLA_Obj TV, fla_bidiagut_t* cntl );

FLA_Error FLA_Bidiag_UT_create_T( FLA_Obj A, FLA_Obj* TU, FLA_Obj* TV );

FLA_Error FLA_Bidiag_UT_recover_tau( FLA_Obj TU, FLA_Obj TV, FLA_Obj tu, FLA_Obj tv );

FLA_Error FLA_Bidiag_UT_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_u_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_l_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e );

FLA_Error FLA_Bidiag_UT_extract_real_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_u_extract_real_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_l_extract_real_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e );

FLA_Error FLA_Bidiag_UT_scale_diagonals( FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Bidiag_UT_u_scale_diagonals( FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Bidiag_UT_l_scale_diagonals( FLA_Obj alpha, FLA_Obj A );

FLA_Error FLA_Bidiag_UT_realify( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_l_realify_unb( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_l_realify_opt( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_u_realify_unb( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_u_realify_opt( FLA_Obj A, FLA_Obj d, FLA_Obj e );

FLA_Error FLA_Bidiag_UT_realify_diagonals( FLA_Uplo uplo, FLA_Obj a, FLA_Obj b, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_realify_diagonals_opt( FLA_Obj a, FLA_Obj b, FLA_Obj d, FLA_Obj e );

FLA_Error FLA_Bidiag_UT_form_U( FLA_Obj A, FLA_Obj T, FLA_Obj U );
FLA_Error FLA_Bidiag_UT_form_V( FLA_Obj A, FLA_Obj S, FLA_Obj V );

FLA_Error FLA_Bidiag_UT_form_U_ext( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T, FLA_Trans transu, FLA_Obj U );
FLA_Error FLA_Bidiag_UT_form_V_ext( FLA_Uplo uplo, FLA_Obj A, FLA_Obj S, FLA_Trans transv, FLA_Obj V );

