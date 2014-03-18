
#include "FLAME.h"

FLA_Error FLA_LQ_UT_solve( FLA_Obj A, FLA_Obj T, FLA_Obj B, FLA_Obj X )
{
  FLA_Obj W;
  FLA_Obj AL, AR;
  FLA_Obj XT, XB;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LQ_UT_solve_check( A, T, B, X );

  FLA_Apply_Q_UT_create_workspace( T, X, &W );

  FLA_Part_1x2( A,   &AL, &AR,    FLA_Obj_length( A ), FLA_LEFT );
  FLA_Part_2x1( X,   &XT,
                     &XB,    FLA_Obj_length( B ), FLA_TOP );

  FLA_Copy_external( B, XT );

  FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                     FLA_NONUNIT_DIAG, FLA_ONE, AL, XT );

  FLA_Set( FLA_ZERO, XB );

  FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                  A, T, W, X );

  FLA_Obj_free( &W );

  return FLA_SUCCESS;
}

