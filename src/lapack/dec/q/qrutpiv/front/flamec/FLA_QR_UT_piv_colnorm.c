
#include "FLAME.h"

FLA_Error FLA_QR_UT_piv_colnorm( FLA_Obj alpha, FLA_Obj A, FLA_Obj b )
{
  FLA_Obj AL,   AR,       A0,  a1,  A2;
  FLA_Obj bT,             b0,
          bB,             beta1, 
                          b2;

  FLA_Obj val2_a1, val2_a1_real;

  // A and b has matching dimensions.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_QR_UT_piv_colnorm_check( alpha, A, b );
  
  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, &val2_a1 );
  FLA_Obj_create( FLA_Obj_datatype( b ), 1, 1, 0, 0, &val2_a1_real );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_2x1( b,    &bT,
                      &bB,            0, FLA_TOP );
  
  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( bT,                &b0,
                        /* ** */            /* ** */
                                              &beta1,
                           bB,                &b2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    // Using dot product is a bit dangerous when a1 is close to 
    // under/over flow limits.
    // The matrix should be properly scaled before using QR_UT_piv.
    FLA_Dot( a1, a1, val2_a1 );
    FLA_Obj_extract_real_part( val2_a1, val2_a1_real );
    FLA_Axpy( alpha, val2_a1_real, beta1 );
    /*------------------------------------------------------------*/
    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &bT,                b0,
                                                  beta1,
                            /* ** */           /* ** */
                              &bB,                b2,     FLA_TOP );
  }

  FLA_Obj_free( &val2_a1 );
  FLA_Obj_free( &val2_a1_real );

  return FLA_SUCCESS;
}

