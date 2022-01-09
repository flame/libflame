/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Trmm( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B );
void time_Trmm(
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Trmm( 
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    C_old, A_flat, C_flat;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, A, &A_flat );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, C, &C_flat );

  FLASH_Copy( C, C_old );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( C_old, C );
    FLASH_Obj_flatten( A, A_flat );
    FLASH_Obj_flatten( C, C_flat );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
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
        REF_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
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
        REF_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
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
        REF_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
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
        REF_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
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
        REF_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 6
    case 6:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 7
    case 7:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 8
    case 8:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 9
    case 9:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 10
    case 10:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 11
    case 11:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A_flat, C_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
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
    FLASH_Obj_hierarchify( C_flat, C_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLASH_Max_elemwise_diff( C, C_ref );
  }

  *gflops = 1.0 * 
            FLASH_Obj_scalar_length( C ) * 
            FLASH_Obj_scalar_width( C ) * 
            FLASH_Obj_scalar_width( A ) / 
            dtime_old / 
            1.0e9;

  if ( param_combo == 0 ||
       param_combo == 3 ||
       param_combo == 6 ||
       param_combo == 9 )
  *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( C_old, C );

  FLASH_Obj_free( &C_old );
  FLASH_Obj_free( &A_flat );
  FLASH_Obj_free( &C_flat );
}

