/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_SA_FS_blk( FLA_Obj L, 
                         FLA_Obj D, FLA_Obj p, FLA_Obj C, 
                                               FLA_Obj E, dim_t nb_alg )
{
  FLA_Obj LT,              L0,
          LB,              L1,
                           L2;

  FLA_Obj DL,    DR,       D0,  D1,  D2;

  FLA_Obj pT,              p0,
          pB,              p1,
                           p2;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  FLA_Obj L1_sqr, L1_rest;

  dim_t b;

  FLA_Part_2x1( L,    &LT, 
                      &LB,            0, FLA_TOP );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );

  FLA_Part_2x1( p,    &pT, 
                      &pB,            0, FLA_TOP );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_TOP );

  while ( FLA_Obj_length( LT ) < FLA_Obj_length( L ) )
  {
    b = fla_min( FLA_Obj_length( LB ), nb_alg );

    FLA_Repart_2x1_to_3x1( LT,                &L0, 
                        /* ** */            /* ** */
                                              &L1, 
                           LB,                &L2,        b, FLA_BOTTOM );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,      &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( pT,                &p0, 
                        /* ** */            /* ** */
                                              &p1, 
                           pB,                &p2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                        /* ** */            /* ** */
                                              &C1, 
                           CB,                &C2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_1x2( L1,    &L1_sqr, &L1_rest,      b, FLA_LEFT );


    FLA_SA_Apply_pivots( C1,
                         E, p1 );

    FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, L1_sqr, C1 );

    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_MINUS_ONE, D1, C1, FLA_ONE, E );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &LT,                L0, 
                                                  L1, 
                            /* ** */           /* ** */
                              &LB,                L2,     FLA_TOP );

    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,     D0, D1, /**/ D2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &pT,                p0, 
                                                  p1, 
                            /* ** */           /* ** */
                              &pB,                p2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                                                  C1, 
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}
