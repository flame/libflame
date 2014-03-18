
#include "FLAME.h"

// The resulting Q is same regardless of the given uplo flag (use the symmetry).
FLA_Error FLA_Tridiag_UT_form_Q( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T, FLA_Obj Q )
{
    FLA_Error r_val = FLA_SUCCESS;
    FLA_Obj   ATL, ATR,
              ABL, ABR;
    FLA_Obj   QTL, QTR,
              QBL, QBR;
    FLA_Obj   TL,  TR;

    if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
        FLA_Tridiag_UT_form_Q_check( uplo, A, T, Q );

    // Adjust T.
    FLA_Part_1x2( T,    &TL, &TR,     1, FLA_RIGHT );

    if ( FLA_Obj_is( A, Q ) == FALSE )
    {
        FLA_Part_2x2( Q,  &QTL, &QTR,
                          &QBL, &QBR,   1, 1, FLA_TL );

        FLA_Set( FLA_ONE,  QTL );
        FLA_Set( FLA_ZERO, QTR );
        FLA_Set( FLA_ZERO, QBL );

        if ( uplo == FLA_LOWER_TRIANGULAR )
        {
            FLA_Part_2x2( A,  &ATL, &ATR,
                              &ABL, &ABR,   1, 1, FLA_TR );
            FLA_QR_UT_form_Q( ABL, TL, QBR );
        }
        else // ( uplo == FLA_UPPER_TRIANGULAR )
        {
            FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
        }
    }
    else
    {
        // Shift the Householder vectors one row/column towards the diagonal.
        FLA_Tridiag_UT_shift_U( uplo, A );

        FLA_Part_2x2( A,  &ATL, &ATR,
                          &ABL, &ABR,   1, 1, FLA_TL );

        if ( uplo == FLA_LOWER_TRIANGULAR )
        {
            FLA_QR_UT_form_Q( ABR, TL, ABR );
        }
        else // ( uplo == FLA_UPPER_TRIANGULAR )
        {
            FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
        }
    }
    return r_val;
}

