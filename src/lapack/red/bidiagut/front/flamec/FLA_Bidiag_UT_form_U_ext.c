/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_form_U_ext( FLA_Uplo  uplo,  FLA_Obj A, FLA_Obj T,
                                    FLA_Trans transu, FLA_Obj U )
{
    //if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    //  FLA_Bidiag_UT_form_U_ext_check( uplo, A, T, transu, U );

    if ( transu == FLA_NO_TRANSPOSE ||
         transu == FLA_CONJ_NO_TRANSPOSE )
    {

        if ( uplo == FLA_UPPER_TRIANGULAR )
        {
            FLA_QR_UT_form_Q( A, T, U );
        }
        else // if ( uplo == FLA_LOWER_TRIANGULAR )
        {
            FLA_Obj ATL, ATR,
                    ABL, ABR;

            FLA_Obj UTL, UTR,
                    UBL, UBR;

            FLA_Obj TL,  TR;

            dim_t   b = ( FLA_Obj_length( A ) - 1 );

            FLA_Part_1x2( T,    &TL,  &TR,     1, FLA_RIGHT );
            FLA_Part_2x2( U,    &UTL, &UTR,
                                &UBL, &UBR,    1, 1, FLA_TL );

            if ( FLA_Obj_is( A, U ) == FALSE )
            {
                FLA_Set( FLA_ONE,  UTL ); FLA_Set( FLA_ZERO, UTR );
                FLA_Set( FLA_ZERO, UBL );

                FLA_Part_2x2( A,    &ATL, &ATR,
                                    &ABL, &ABR,    1, b, FLA_TL );

                FLA_QR_UT_form_Q( ABL, TL, UBR );
            }
            else
            {
                FLA_Obj p, pt, pb;
                FLA_Part_2x2( A,    &ATL, &ATR,
                                    &ABL, &ABR,    1, b+1, FLA_TL );

                FLA_Obj_create( FLA_INT, b+1,1, 0, 0, &p );
                FLA_Part_2x1( p,    &pt,
                                    &pb, 1, FLA_BOTTOM );
                FLA_Set( FLA_ONE, pt );
                FLA_Set( FLA_ZERO, pb );
                FLA_Apply_pivots ( FLA_RIGHT, FLA_NO_TRANSPOSE, p, ABL );
                FLA_Obj_free(&p );

                FLA_Set( FLA_ONE,  UTL );
                FLA_Set( FLA_ZERO, UBL );
                FLA_Set( FLA_ZERO, UTR );

                FLA_Part_1x2( ABL,  &ABL,
                                    &ABR, 1, FLA_LEFT );

                FLA_QR_UT_form_Q( ABR, TL, UBR );
            }
        }
    }
    else
    {
        FLA_Uplo uplo_flip = ( uplo == FLA_UPPER_TRIANGULAR ?
                               FLA_LOWER_TRIANGULAR : FLA_UPPER_TRIANGULAR );

        FLA_Obj_flip_base( &A );
        FLA_Obj_flip_view( &A );

        // A and U should have different base objects
        FLA_Bidiag_UT_form_V_ext( uplo_flip,  A, T,
                                  FLA_CONJ_TRANSPOSE, U );

        FLA_Obj_flip_base( &A );

        // As we use QR and LQ for constructing U and V, 
        // conjugation naturally fits there. 
        // Never apply conjugation separately here even if flipping trick is applied.
        //FLA_Conjugate( U );
    }

    return FLA_SUCCESS;
}

