
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Gemv_t_blk_var2( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;

  FLA_Obj yT,              y0,
          yB,              y1,
                           y2;

  dim_t b;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_RIGHT );

  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_BOTTOM );

  while ( FLA_Obj_width( AR ) < FLA_Obj_width( A ) ){

    b = FLA_Determine_blocksize( AL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, &A1, /**/ &A2,
                           b, FLA_LEFT );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                                              &y1, 
                        /* ** */            /* ** */
                           yB,                &y2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    /* y1 = alpha * A1' * x + y1 */
    FLA_Gemv_internal( FLA_TRANSPOSE, 
                       alpha, A1, x, beta, y1,
                       FLA_Cntl_sub_gemv( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, /**/ A1, A2,
                              FLA_RIGHT );

    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                            /* ** */           /* ** */
                                                  y1, 
                              &yB,                y2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
