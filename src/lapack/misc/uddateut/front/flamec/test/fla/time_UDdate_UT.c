/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


void time_UDdate_UT(
                 integer variant, integer type, integer n_repeats, integer mB, integer mC, integer mD, integer n,
                 FLA_Obj B, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj R, FLA_Obj E,
                 double *dtime, double *diff, double *gflops );


void time_UDdate_UT(
                 integer variant, integer type, integer n_repeats, integer mB, integer mC, integer mD, integer n,
                 FLA_Obj B, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj R, FLA_Obj E,
                 double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    RR, EE, R_save, C_save, D_save;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, R, &RR );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, E, &EE );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, R, &R_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, D, &D_save );

  FLA_Copy_external( R, R_save );
  FLA_Copy_external( C, C_save );
  FLA_Copy_external( D, D_save );

//FLA_Obj_show( "B", B, "%10.3e", "" );
//FLA_Obj_show( "C", C, "%10.3e", "" );
//FLA_Obj_show( "D", D, "%10.3e", "" );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( R_save, R );
    FLA_Copy_external( C_save, C );
    FLA_Copy_external( D_save, D );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:{
      switch( type ){
      case FLA_ALG_FRONT:
        FLA_UDdate_UT( R, C, D, T );
        //FLA_UDdate_UT_unb_var1( R, C, D, T );
        //FLA_UDdate_UT_opt_var1( R, C, D, T );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = fla_min( *dtime, dtime_old );

  }

  {
//FLA_Obj_show( "R", R, "%10.3e", "" );
//FLA_Obj_show( "E", E, "%10.3e", "" );
//FLA_Obj_show( "R", R, "%10.3e + %10.3e", "" );
//FLA_Obj_show( "E", E, "%10.3e + %10.3e", "" );

    FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, R, R, FLA_ZERO, RR );
    FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, E, E, FLA_ZERO, EE );

    *diff = FLA_Max_elemwise_diff( RR, EE );
  }

  *gflops = 2 * ( ( mC + mD ) * n * n +
                  ( mC + mD ) * n * 6 ) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( R ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( R_save, R );
  FLA_Copy_external( C_save, C );
  FLA_Copy_external( D_save, D );

  FLA_Obj_free( &R_save );
  FLA_Obj_free( &C_save );
  FLA_Obj_free( &D_save );

  FLA_Obj_free( &RR );
  FLA_Obj_free( &EE );
}

