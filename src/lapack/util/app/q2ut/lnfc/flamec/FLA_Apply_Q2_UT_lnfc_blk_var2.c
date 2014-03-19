/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_Q2_UT_lnfc_blk_var2( FLA_Obj D, FLA_Obj T, FLA_Obj W1, FLA_Obj C, 
                                                                           FLA_Obj E, fla_apq2ut_t* cntl )
{
  FLA_Obj DT,              D0,
          DB,              D1,
                           D2;

  FLA_Obj TT,              T0,
          TB,              T1,
                           T2;

  FLA_Obj ET,              E0,
          EB,              E1,
                           E2;

  dim_t b;

  FLA_Part_2x1( D,    &DT, 
                      &DB,            0, FLA_BOTTOM );

  FLA_Part_2x1( T,    &TT, 
                      &TB,            0, FLA_BOTTOM );

  FLA_Part_2x1( E,    &ET, 
                      &EB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( DB ) < FLA_Obj_length( D ) ){

    b = FLA_Determine_blocksize( DT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( DT,                &D0, 
                                              &D1, 
                        /* ** */            /* ** */
                           DB,                &D2,        b, FLA_TOP );

    FLA_Repart_2x1_to_3x1( TT,                &T0, 
                                              &T1, 
                        /* ** */            /* ** */
                           TB,                &T2,        b, FLA_TOP );

    FLA_Repart_2x1_to_3x1( ET,                &E0, 
                                              &E1, 
                        /* ** */            /* ** */
                           EB,                &E2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    //  / C  \ =  Q / C  \
    //  \ E1 /      \ E1 /
    //
    //  where Q is formed from D1 and T1.

    FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                              D1, T1, W1, C,
                                          E1, FLA_Cntl_sub_apq2ut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &DT,                D0, 
                            /* ** */           /* ** */
                                                  D1, 
                              &DB,                D2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &TT,                T0, 
                            /* ** */           /* ** */
                                                  T1, 
                              &TB,                T2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &ET,                E0, 
                            /* ** */           /* ** */
                                                  E1, 
                              &EB,                E2,     FLA_BOTTOM );
  }

  return FLA_SUCCESS;
}

