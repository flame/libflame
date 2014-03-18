
#include "FLAME.h"

FLA_Error FLA_Apply_Q2_UT_lhfc_blk_var2( FLA_Obj D, FLA_Obj T, FLA_Obj W1, FLA_Obj C, 
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
                      &DB,            0, FLA_TOP );

  FLA_Part_2x1( T,    &TT, 
                      &TB,            0, FLA_TOP );

  FLA_Part_2x1( E,    &ET, 
                      &EB,            0, FLA_TOP );

  while ( FLA_Obj_length( DT ) < FLA_Obj_length( D ) ){

    b = FLA_Determine_blocksize( DB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( DT,                &D0, 
                        /* ** */            /* ** */
                                              &D1, 
                           DB,                &D2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( TT,                &T0, 
                        /* ** */            /* ** */
                                              &T1, 
                           TB,                &T2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( ET,                &E0, 
                        /* ** */            /* ** */
                                              &E1, 
                           EB,                &E2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    //  / C  \ =  Q' / C  \
    //  \ E1 /       \ E1 /
    //
    //  where Q is formed from D1 and T1.

    FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                              D1, T1, W1, C,
                                          E1, FLA_Cntl_sub_apq2ut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &DT,                D0, 
                                                  D1, 
                            /* ** */           /* ** */
                              &DB,                D2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &TT,                T0, 
                                                  T1, 
                            /* ** */           /* ** */
                              &TB,                T2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &ET,                E0, 
                                                  E1, 
                            /* ** */           /* ** */
                              &EB,                E2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

