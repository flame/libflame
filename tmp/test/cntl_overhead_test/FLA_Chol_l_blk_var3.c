
#include "FLAME.h"

FLA_Error FLA_Chol_l_blk_var3_ht( FLA_Obj A, FLA_Chol_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  int b, value = 0;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    b = 1;

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    /* ********************************************************************* */

    if ( FLA_Obj_is_real( A ) )
    {
      FLA_Chol_unb_external( FLA_LOWER_TRIANGULAR, *FLASH_OBJ_PTR_AT( A11 ), NULL );

      if ( value != FLA_SUCCESS )
        return ( FLA_Obj_length( A00 ) + value );

      FLA_Trsm_rlt_blk_var3_ht( FLA_NONUNIT_DIAG, FLA_ONE, A11, A21,
                             NULL );

      FLA_Syrk_ln_blk_var2_ht( FLA_MINUS_ONE, A21, FLA_ONE, A22,
                            NULL );
    }
    else
    {
      value = FLA_Chol_internal( FLA_LOWER_TRIANGULAR, A11, 
                                 NULL );

      if ( value != FLA_SUCCESS )
        return ( FLA_Obj_length( A00 ) + value );

      FLA_Trsm_internal( FLA_RIGHT, FLA_LOWER_TRIANGULAR, 
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, 
                         FLA_ONE, A11, A21,
                         NULL );

      FLA_Herk_internal( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, 
                         FLA_MINUS_ONE, A21, FLA_ONE, A22,
                         NULL );
    }

    /* ********************************************************************* */

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
  }

  return value;
}

