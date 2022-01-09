/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


void time_UDdate_UT_inc(
                 integer variant, integer type, integer n_repeats, integer mB, integer mC, integer mD, integer n,
                 FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj E,
                 double *dtime, double *diff, double *gflops );


void time_UDdate_UT_inc(
                 integer variant, integer type, integer n_repeats, integer mB, integer mC, integer mD, integer n,
                 FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj E,
                 double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    RR, EE, R_save, C_save, D_save;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, R, &RR );
  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, E, &EE );

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, R, &R_save );
  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_save );
  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, D, &D_save );

  FLASH_Copy( R, R_save );
  FLASH_Copy( C, C_save );
  FLASH_Copy( D, D_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLASH_Copy( R_save, R );
    FLASH_Copy( C_save, C );
    FLASH_Copy( D_save, D );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:{
      switch( type ){
      case FLA_ALG_FRONT:
        FLASH_UDdate_UT_inc( R, C, D, T, W );

        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );

  }

  {
    FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, R, R, FLA_ZERO, RR );
    FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, E, E, FLA_ZERO, EE );

    *diff = FLASH_Max_elemwise_diff( RR, EE );
  }

  *gflops = 2 * ( ( mC + mD ) * n * n +
                  ( mC + mD ) * n * 6 ) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( R ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( R_save, R );
  FLASH_Copy( C_save, C );
  FLASH_Copy( D_save, D );

  FLASH_Obj_free( &R_save );
  FLASH_Obj_free( &C_save );
  FLASH_Obj_free( &D_save );

  FLASH_Obj_free( &RR );
  FLASH_Obj_free( &EE );
}

