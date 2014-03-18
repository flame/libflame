
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_form_V_ext( FLA_Uplo  uplo,  FLA_Obj A, FLA_Obj S,
                                    FLA_Trans transv, FLA_Obj V )
{
  //if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
  //      FLA_Bidiag_UT_form_V_ext_check( uplo, A, S, transv, V );


    if ( transv == FLA_TRANSPOSE ||
         transv == FLA_CONJ_TRANSPOSE )
    {
        if ( uplo == FLA_UPPER_TRIANGULAR )
        {
            FLA_Obj ATL, ATR,
                    ABL, ABR;

            FLA_Obj VTL, VTR,
                    VBL, VBR;

            FLA_Obj SL,  SR;

            dim_t   b = ( FLA_Obj_width( A ) - 1 );

            FLA_Part_1x2( S,    &SL,  &SR,     1, FLA_RIGHT );
            FLA_Part_2x2( V,    &VTL, &VTR,
                                &VBL, &VBR,    1, 1, FLA_TL );

            if ( FLA_Obj_is( A, V ) == FALSE )
            {
                FLA_Set( FLA_ONE,  VTL ); FLA_Set( FLA_ZERO, VTR );
                FLA_Set( FLA_ZERO, VBL );


                FLA_Part_2x2( A,    &ATL, &ATR,
                                    &ABL, &ABR,    b, b, FLA_TR );

                FLA_LQ_UT_form_Q( ATR, SL, VBR );
            }
            else
            {
                FLA_Obj p, pt, pb;
                FLA_Part_2x2( A,    &ATL, &ATR,
                                    &ABL, &ABR,    b+1, b, FLA_TR );

                FLA_Obj_create( FLA_INT, b+1, 1, 0, 0, &p );
                FLA_Part_2x1( p,    &pt,
                                    &pb, 1, FLA_BOTTOM );
                FLA_Set( FLA_ONE, pt );
                FLA_Set( FLA_ZERO, pb );

                FLA_Apply_pivots ( FLA_LEFT, FLA_TRANSPOSE, p, ATR );
                FLA_Obj_free( &p );

                FLA_Set( FLA_ONE,  VTL );
                FLA_Set( FLA_ZERO, VBL );
                FLA_Set( FLA_ZERO, VTR );

                FLA_Part_2x1( ATR,  &ATR,
                                    &ABR, 1, FLA_TOP );

                FLA_LQ_UT_form_Q( ABR, SL, VBR );
            }
        }
        else // if ( uplo == FLA_LOWER_TRIANGULAR )
        {
            FLA_LQ_UT_form_Q( A, S, V );
        }
    }
    else
    {
        FLA_Uplo uplo_flip = ( uplo == FLA_UPPER_TRIANGULAR ?
                               FLA_LOWER_TRIANGULAR : FLA_UPPER_TRIANGULAR );

        FLA_Obj_flip_base( &A );
        FLA_Obj_flip_view( &A );

        // A and U should have different base objects.
        FLA_Bidiag_UT_form_U_ext( uplo_flip, A, S,
                                  FLA_NO_TRANSPOSE, V );

        FLA_Obj_flip_base( &A );

        // As we use QR and LQ for constructing U and V, 
        // conjugation naturally fits there. 
        // Never apply conjugation separately here even if flipping trick is applied.
        // FLA_Conjugate( V );
    }

    return FLA_SUCCESS;
}

