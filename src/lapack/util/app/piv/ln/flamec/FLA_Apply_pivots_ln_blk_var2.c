/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Apply_pivots_ln_blk_var2_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj pT,              p0,
          pB,              pi1,
                           p2;

  dim_t b;

  FLA_Part_2x1_ts( FLA_cntl_init_i, A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1_ts( FLA_cntl_init_i, p,    &pT, 
                      &pB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ) {

    b = FLA_Determine_blocksize_ts( FLA_cntl_init_i, AB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1_ts( FLA_cntl_init_i, AT,                &A0, 
                        /* ** */            /* ** */
                                              &A1, 
                           AB,                &A2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1_ts( FLA_cntl_init_i, pT,                &p0, 
                        /* ** */            /* ** */
                                              &pi1, 
                           pB,                &p2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    /* Apply pivots to a block and matrix */
    FLA_Apply_pivots_internal_ts( FLA_cntl_init_i, FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB,
                               FLA_Cntl_sub_appiv( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1_ts( FLA_cntl_init_i, &AT,                A0, 
                                                  A1, 
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1_ts( FLA_cntl_init_i, &pT,                p0, 
                                                  pi1, 
                            /* ** */           /* ** */
                              &pB,                p2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}
#endif

FLA_Error FLA_Apply_pivots_ln_blk_var2( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj pT,              p0,
          pB,              pi1,
                           p2;

  dim_t b;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( p,    &pT, 
                      &pB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ) {

    b = FLA_Determine_blocksize( AB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* ** */
                                              &A1, 
                           AB,                &A2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( pT,                &p0, 
                        /* ** */            /* ** */
                                              &pi1, 
                           pB,                &p2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    /* Apply pivots to a block and matrix */
    FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB,
                               FLA_Cntl_sub_appiv( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  A1, 
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &pT,                p0, 
                                                  pi1, 
                            /* ** */           /* ** */
                              &pB,                p2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

