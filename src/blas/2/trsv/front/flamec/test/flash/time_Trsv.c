/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Trsv( FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
void time_Trsv(
               integer param_combo, integer type, integer nrepeats, integer m,
               FLA_Obj A, FLA_Obj x, FLA_Obj x_ref,
               double *dtime, double *diff, double *gflops );


void time_Trsv( 
               integer param_combo, integer type, integer nrepeats, integer m,
               FLA_Obj A, FLA_Obj x, FLA_Obj x_ref,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    x_old, A_flat, x_flat;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, x, &x_old );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, A, &A_flat );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, x, &x_flat );

  FLASH_Copy( x, x_old );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( x_old, x );
    FLASH_Obj_flatten( A, A_flat );
    FLASH_Obj_flatten( x, x_flat );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trsv( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, A_flat, x_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trsv( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
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
        REF_Trsv( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A_flat, x_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trsv( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
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
        REF_Trsv( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A_flat, x_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trsv( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
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
        REF_Trsv( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, A_flat, x_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trsv( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 4
    case 4:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A_flat, x_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 5
    case 5:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trsv( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A_flat, x_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trsv( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
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
    FLASH_Obj_hierarchify( x_flat, x_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLASH_Max_elemwise_diff( x, x_ref );
  }

  *gflops = 1.0 * m * m /
            dtime_old / 
            1.0e9;

  if ( param_combo == 0 ||
       param_combo == 3 )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( x_old, x );

  FLASH_Obj_free( &x_old );
  FLASH_Obj_free( &A_flat );
  FLASH_Obj_free( &x_flat );
}

