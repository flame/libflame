/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"



#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_SPDinv( FLA_Trans trans, FLA_Obj C );

void time_SPDinv(
                integer param_combo, integer type, integer nrepeats, integer m,
                FLA_Obj C, FLA_Obj C_ref,
                double *dtime, double *diff, double *gflops );


void time_SPDinv(
                integer param_combo, integer type, integer nrepeats, integer m,
                FLA_Obj C, FLA_Obj C_ref,
                double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    C_old, C_flat;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, C, &C_flat );

  FLASH_Copy( C, C_old );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( C_old, C );
    FLASH_Obj_flatten( C, C_flat );

    *dtime = FLA_Clock();

    switch( param_combo ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_SPDinv( FLA_LOWER_TRIANGULAR, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_SPDinv( FLA_LOWER_TRIANGULAR, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 1:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_SPDinv( FLA_UPPER_TRIANGULAR, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_SPDinv( FLA_UPPER_TRIANGULAR, C );
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

  if ( type == FLA_ALG_REFERENCE ){
    FLASH_Obj_hierarchify( C_flat, C_ref );
    *diff = 0.0;
  }
  else{
    *diff = FLASH_Max_elemwise_diff( C, C_ref );
  }

  *gflops = 1.0 * 
            FLASH_Obj_scalar_length( C ) * 
            FLASH_Obj_scalar_length( C ) * 
            FLASH_Obj_scalar_length( C ) / 
            dtime_old / 1e9;

  *dtime = dtime_old;

  FLASH_Copy( C_old, C );

  FLASH_Obj_free( &C_old );
  FLASH_Obj_free( &C_flat );
}

