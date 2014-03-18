
#include "FLAME.h"

FLA_Error FLA_Norm1( FLA_Obj A, FLA_Obj norm )
{
  FLA_Obj AL,   AR,       A0,  a1,  A2;

  FLA_Obj b;
  FLA_Obj bL,   bR,       b0,  beta1,  b2;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Norm1_check( A, norm );

  FLA_Obj_create( FLA_Obj_datatype( A ), 1, FLA_Obj_width( A ), 0, 0, &b );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_1x2( b,    &bL,  &bR,      0, FLA_LEFT );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( bL,  /**/ bR,        &b0, /**/ &beta1, &b2,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Asum( a1, beta1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &bL,  /**/ &bR,        b0, beta1, /**/ b2,
                              FLA_LEFT );
  }

  FLA_Max_abs_value( b, norm );

  FLA_Obj_free( &b );

  return FLA_SUCCESS;
}

