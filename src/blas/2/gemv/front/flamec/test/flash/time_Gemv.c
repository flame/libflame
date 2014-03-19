/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Gemv( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
void time_Gemv(
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj x, FLA_Obj y, FLA_Obj y_ref,
               double *dtime, double *diff, double *gflops );


void time_Gemv( 
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj x, FLA_Obj y, FLA_Obj y_ref,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    y_old, A_flat, x_flat, y_flat;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, y, &y_old );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, A, &A_flat );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, x, &x_flat );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, y, &y_flat );

  FLASH_Copy( y, y_old );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( y_old, y );
    FLASH_Obj_flatten( A, A_flat );
    FLASH_Obj_flatten( x, x_flat );
    FLASH_Obj_flatten( y, y_flat );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A_flat, x_flat, FLA_ONE, y_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A, x, FLA_ONE, y );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 1
    case 1:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A_flat, x_flat, FLA_ONE, y_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A, x, FLA_ONE, y );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 2
    case 2:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Gemv( FLA_TRANSPOSE, FLA_ONE, A_flat, x_flat, FLA_ONE, y_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Gemv( FLA_TRANSPOSE, FLA_ONE, A, x, FLA_ONE, y );
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


  if ( type == FLA_ALG_REFERENCE )
  {
    FLASH_Obj_hierarchify( y_flat, y_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLASH_Max_elemwise_diff( y, y_ref );
  }

  *gflops = 2.0 * m * n / 
            dtime_old / 
            1.0e9;

  if ( param_combo == 0 )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( y_old, y );

  FLASH_Obj_free( &y_old );
  FLASH_Obj_free( &A_flat );
  FLASH_Obj_free( &x_flat );
  FLASH_Obj_free( &y_flat );
}

