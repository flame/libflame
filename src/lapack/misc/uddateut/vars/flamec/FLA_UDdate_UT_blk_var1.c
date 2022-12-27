/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_UDdate_UT_blk_var1( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl )
{
  FLA_Obj RTL,   RTR,      R00, R01, R02, 
          RBL,   RBR,      R10, R11, R12,
                           R20, R21, R22;

  FLA_Obj CL,    CR,       C0,  C1,  C2;

  FLA_Obj DL,    DR,       D0,  D1,  D2;

  FLA_Obj TL,    TR,       T0,  T1,  W12;

  FLA_Obj T1T,  T1B;

  FLA_Obj W12T, W12B;

  dim_t   b_alg, b;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  FLA_Part_2x2( R,    &RTL, &RTR,
                      &RBL, &RBR,     0, 0, FLA_TL );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  while ( FLA_Obj_min_dim( RBR ) > 0 ){

    b = fla_min( b_alg, FLA_Obj_min_dim( RBR ) );

    FLA_Repart_2x2_to_3x3( RTL, /**/ RTR,       &R00, /**/ &R01, &R02,
                        /* ************* */   /* ******************** */
                                                &R10, /**/ &R11, &R12,
                           RBL, /**/ RBR,       &R20, /**/ &R21, &R22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &W12,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,   &T1T,
                        &T1B,   b, FLA_TOP );

    /*
       Perform an up/downdate of the upper triangular factor R11 via 
       up/downdating UT Householder transformations:

         [ R11, ...
           C1,  ...
           D1, T1T ] = FLA_UDdate_UT( R11, ...
                                      C1, ...
                                      D1, T1T );

       by updating R11 in such a way that removes the contributions of the rows
       in D1 while simultaneously adding new contributions to the factorization
       from the rows of C1. Note that C1 and D1 are also updated in the process.
    */

    FLA_UDdate_UT_internal( R11,
                            C1,
                            D1, T1T,
                            FLA_Cntl_sub_uddateut( cntl ) );


    if ( FLA_Obj_width( R12 ) > 0 )
    {
      FLA_Part_2x1( W12,  &W12T, 
                          &W12B,  b, FLA_TOP );

      /*
         Apply Q' to R12, C2, and D2 from the left:

           / R12 \          / R12 \
           | C2  |  =  Q' * | C2  |
           \ D2  /          \ D2  /

         where Q is formed from C1 and D1.
      */

      FLA_Apply_QUD_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                 T1T, W12T,
                                      R12,
                                 C1,  C2,
                                 D1,  D2, FLA_Cntl_sub_apqudut( cntl ) );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &RTL, /**/ &RTR,       R00, R01, /**/ R02,
                                                     R10, R11, /**/ R12,
                            /* ************** */  /* ****************** */
                              &RBL, /**/ &RBR,       R20, R21, /**/ R22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ W12,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

