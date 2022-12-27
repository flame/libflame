/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Symm( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C );
void time_Symm(
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Symm( 
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    C_old, A_flat, B_flat, C_flat;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, A, &A_flat );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, B, &B_flat );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, C, &C_flat );

  FLASH_Copy( C, C_old );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( C_old, C );
    FLASH_Obj_flatten( A, A_flat );
    FLASH_Obj_flatten( B, B_flat );
    FLASH_Obj_flatten( C, C_flat );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Symm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_ONE, A_flat, B_flat, FLA_ZERO, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Symm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_ONE, A, B, FLA_ZERO, C );
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
        REF_Symm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_ONE, A_flat, B_flat, FLA_ZERO, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Symm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_ONE, A, B, FLA_ZERO, C );
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
        REF_Symm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_ONE, A_flat, B_flat, FLA_ZERO, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Symm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 3
    case 3:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Symm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_ONE, A_flat, B_flat, FLA_ZERO, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Symm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_ONE, A, B, FLA_ZERO, C );
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


  if ( type == FLA_ALG_REFERENCE )
  {
    FLASH_Obj_hierarchify( C_flat, C_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLASH_Max_elemwise_diff( C, C_ref );
  }

  *gflops = 2.0 * 
            FLASH_Obj_scalar_length( C ) * 
            FLASH_Obj_scalar_width( C ) * 
            FLASH_Obj_scalar_width( A ) / 
            dtime_old / 
            1.0e9;

  if ( FLA_Obj_is_complex( C ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( C_old, C );

  FLASH_Obj_free( &C_old );
  FLASH_Obj_free( &A_flat );
  FLASH_Obj_free( &B_flat );
  FLASH_Obj_free( &C_flat );
}

