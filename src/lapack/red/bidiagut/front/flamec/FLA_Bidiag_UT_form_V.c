/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_form_V( FLA_Obj A, FLA_Obj S, FLA_Obj V )
{
    FLA_Uplo uplo;

    if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
        FLA_Bidiag_UT_form_V_check( A, S, V );

    uplo = ( FLA_Obj_length( A ) >= FLA_Obj_width( A ) ?
             FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR );

    FLA_Bidiag_UT_form_V_ext( uplo, A, S, 
                              FLA_NO_TRANSPOSE, V );
    
    return FLA_SUCCESS;
}

