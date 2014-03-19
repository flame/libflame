/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Norm_inf( FLA_Obj A, FLA_Obj norm )
{
  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;

  FLA_Obj bT,              b0,
          bB,              beta1,
                           b2;
  FLA_Obj b;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Norm_inf_check( A, norm );

  FLA_Obj_create( FLA_Obj_datatype( A ), FLA_Obj_length( A ), 1, 0, 0, &b );

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( b,    &bT, 
                      &bB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* *** */
                                              &a1t, 
                           AB,                &A2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( bT,                &b0, 
                        /* ** */            /* ***** */
                                              &beta1, 
                           bB,                &b2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Asum( a1t, beta1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  a1t, 
                            /* ** */           /* *** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &bT,                b0, 
                                                  beta1, 
                            /* ** */           /* ***** */
                              &bB,                b2,     FLA_TOP );

  }

  FLA_Max_abs_value( b, norm );

  FLA_Obj_free( &b );

  return FLA_SUCCESS;
}

