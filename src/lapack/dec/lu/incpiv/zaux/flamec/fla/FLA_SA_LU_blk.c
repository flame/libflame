/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_SA_LU_blk( FLA_Obj U, 
                         FLA_Obj D, FLA_Obj p, FLA_Obj L, dim_t nb_alg )
{
  FLA_Obj UTL,   UTR,      U00, U01, U02, 
          UBL,   UBR,      U10, U11, U12,
                           U20, U21, U22;

  FLA_Obj DL,    DR,       D0,  D1,  D2;

  FLA_Obj pT,              p0,
          pB,              p1,
                           p2;

  FLA_Obj LT,              L0,
          LB,              L1,
                           L2;

  FLA_Obj L1_sqr, L1_rest;

  dim_t b;

  FLA_Part_2x2( U,    &UTL, &UTR,
                      &UBL, &UBR,     0, 0, FLA_TL );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );

  FLA_Part_2x1( p,    &pT, 
                      &pB,            0, FLA_TOP );

  FLA_Part_2x1( L,    &LT, 
                      &LB,            0, FLA_TOP );

  while ( FLA_Obj_length( UTL ) < FLA_Obj_length( U ) )
  {
    b = fla_min( FLA_Obj_length( UBR ), nb_alg );

    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00, /**/ &U01, &U02,
                        /* ************* */   /* ******************** */
                                                &U10, /**/ &U11, &U12,
                           UBL, /**/ UBR,       &U20, /**/ &U21, &U22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( pT,                  &p0, 
                        /* ** */              /* ** */
                                                &p1, 
                           pB,                  &p2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( LT,                  &L0, 
                        /* ** */              /* ** */
                                                &L1, 
                           LB,                  &L2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_1x2( L1,    &L1_sqr, &L1_rest,      b, FLA_LEFT );


    FLA_SA_LU_unb( U11,
                   D1, p1, L1_sqr );

    FLA_SA_Apply_pivots( U12,
                         D2, p1 );

    FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, L1_sqr, U12 );

    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_MINUS_ONE, D1, U12, FLA_ONE, D2 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00, U01, /**/ U02,
                                                     U10, U11, /**/ U12,
                            /* ************** */  /* ****************** */
                              &UBL, /**/ &UBR,       U20, U21, /**/ U22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &pT,                   p0, 
                                                     p1, 
                            /* ** */              /* ** */
                              &pB,                   p2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &LT,                   L0, 
                                                     L1, 
                            /* ** */              /* ** */
                              &LB,                   L2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}
