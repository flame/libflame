/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"



#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Chol( FLA_Trans trans, FLA_Obj A );

void time_Chol(
                integer param_combo, integer type, integer nrepeats, integer m,
                FLA_Obj A, FLA_Obj A_ref,
                double *dtime, double *diff, double *gflops );


void time_Chol(
                integer param_combo, integer type, integer nrepeats, integer m,
                FLA_Obj A, FLA_Obj A_ref,
                double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_old, A_flat;

  /*if ( type == FLA_ALG_REFERENCE )
  {
    *gflops = *diff = 0.0;
    return;
  }*/

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_old );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, A, &A_flat );

  FLASH_Copy( A, A_old );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( A_old, A );
    FLASH_Obj_flatten( A, A_flat );

    *dtime = FLA_Clock();

    switch( param_combo ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Chol( FLA_LOWER_TRIANGULAR, A_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Chol( FLA_LOWER_TRIANGULAR, A );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 1:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Chol( FLA_UPPER_TRIANGULAR, A_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Chol( FLA_UPPER_TRIANGULAR, A );
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

  if ( type == FLA_ALG_REFERENCE ){
    FLASH_Obj_hierarchify( A_flat, A_ref );
    *diff = 0.0;
  }
  else{
    *diff = FLASH_Max_elemwise_diff( A, A_ref );
  }

  *gflops = 1.0 / 3.0 * m * m * m /
            dtime_old / 1e9;

  *dtime = dtime_old;

  FLASH_Copy( A_old, A );

  FLASH_Obj_free( &A_old );
  FLASH_Obj_free( &A_flat );
}

