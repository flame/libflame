
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

