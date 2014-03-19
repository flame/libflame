/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Gemv_h_blk_var6( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj xT,              x0,
          xB,              x1,
                           x2;

  dim_t b;

  FLA_Scal_internal( beta, y,
                     FLA_Cntl_sub_scal( cntl ) );

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_BOTTOM );

  FLA_Part_2x1( x,    &xT, 
                      &xB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){

    b = FLA_Determine_blocksize( AT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                              &A1, 
                        /* ** */            /* ** */
                           AB,                &A2,        b, FLA_TOP );

    FLA_Repart_2x1_to_3x1( xT,                &x0, 
                                              &x1, 
                        /* ** */            /* ** */
                           xB,                &x2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    /* y = alpha * A1' * x1 + y; */
    FLA_Gemv_internal( FLA_CONJ_TRANSPOSE, 
                       alpha, A1, x1, FLA_ONE, y,
                       FLA_Cntl_sub_gemv( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                            /* ** */           /* ** */
                                                  A1, 
                              &AB,                A2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &xT,                x0, 
                            /* ** */           /* ** */
                                                  x1, 
                              &xB,                x2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
