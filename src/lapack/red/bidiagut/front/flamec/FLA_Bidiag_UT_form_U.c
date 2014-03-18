
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_form_U( FLA_Obj A, FLA_Obj T, FLA_Obj U )
{
    FLA_Uplo uplo;

    if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
        FLA_Bidiag_UT_form_U_check( A, T, U );

    uplo = ( FLA_Obj_length( A ) >= FLA_Obj_width( A ) ?
             FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR );

    FLA_Bidiag_UT_form_U_ext( uplo, A, T, 
                              FLA_NO_TRANSPOSE, U );
    
    return FLA_SUCCESS;
}

