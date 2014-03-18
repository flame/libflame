
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Syrk_lt_blk_var6( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  dim_t b;

  FLA_Scalr_internal( FLA_LOWER_TRIANGULAR, beta, C,
                      FLA_Cntl_sub_scalr( cntl ) );

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){

    b = FLA_Determine_blocksize( AT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                              &A1, 
                        /* ** */            /* ** */
                           AB,                &A2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    /* C = C + A1' * A1 */
    FLA_Syrk_internal( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, 
                       alpha, A1, FLA_ONE, C,
                       FLA_Cntl_sub_syrk( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                            /* ** */           /* ** */
                                                  A1, 
                              &AB,                A2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
