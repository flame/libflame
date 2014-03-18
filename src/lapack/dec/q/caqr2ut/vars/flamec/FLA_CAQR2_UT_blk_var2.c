
#include "FLAME.h"

FLA_Error FLA_CAQR2_UT_blk_var2( FLA_Obj U,
                                 FLA_Obj D, FLA_Obj T, fla_caqr2ut_t* cntl )
{
  FLA_Obj DT,              D0,
          DB,              D1,
                           D2;

  FLA_Obj TT,              T0,
          TB,              T1,
                           T2;

  dim_t b;

  FLA_Part_2x1( D,    &DT, 
                      &DB,            0, FLA_TOP );

  FLA_Part_2x1( T,    &TT, 
                      &TB,            0, FLA_TOP );

  while ( FLA_Obj_length( DT ) < FLA_Obj_length( D ) ){

    b = FLA_Determine_blocksize( DB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( DT,                &D0, 
                        /* ** */            /* ****** */
                                              &D1, 
                           DB,                &D2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( TT,                &T0, 
                        /* ** */            /* ****** */
                                              &T1, 
                           TB,                &T2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // [ U, ...
    //   D1, T ] = FLA_CAQR2_UT( U
    //                           D1, T1 );

    FLA_CAQR2_UT_internal( U,
                           D1, T1, 
                           FLA_Cntl_sub_caqr2ut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &DT,                D0, 
                                                  D1, 
                            /* ** */           /* ****** */
                              &DB,                D2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &TT,                T0, 
                                                  T1, 
                            /* ** */           /* ****** */
                              &TB,                T2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

