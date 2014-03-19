/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Copyr_u_blk_var4( FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj BTL,   BTR,      B00, B01, B02, 
          BBL,   BBR,      B10, B11, B12,
                           B20, B21, B22;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_BR );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_BR );

  while ( FLA_Obj_min_dim( ATL ) > 0 ){

    b = FLA_Determine_blocksize( ATL, FLA_TL, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, &A01, /**/ &A02,
                                                &A10, &A11, /**/ &A12,
                        /* ************* */   /* ******************** */
                           ABL, /**/ ABR,       &A20, &A21, /**/ &A22,
                           b, b, FLA_TL );

    FLA_Repart_2x2_to_3x3( BTL, /**/ BTR,       &B00, &B01, /**/ &B02,
                                                &B10, &B11, /**/ &B12,
                        /* ************* */   /* ******************** */
                           BBL, /**/ BBR,       &B20, &B21, /**/ &B22,
                           b, b, FLA_TL );

    /*------------------------------------------------------------*/

    // B11 = triu( A11 );
    FLA_Copyr_internal( FLA_UPPER_TRIANGULAR, A11, B11,
                        FLA_Cntl_sub_copyr( cntl ) );

    // B01 = A01;
    FLA_Copy_internal( A01, B01,
                       FLA_Cntl_sub_copy( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                            /* ************** */  /* ****************** */
                                                     A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                              FLA_BR );

    FLA_Cont_with_3x3_to_2x2( &BTL, /**/ &BTR,       B00, /**/ B01, B02,
                            /* ************** */  /* ****************** */
                                                     B10, /**/ B11, B12,
                              &BBL, /**/ &BBR,       B20, /**/ B21, B22,
                              FLA_BR );
  }

  return FLA_SUCCESS;
}

