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
               int param_combo, int type, int nrepeats, int m,
               FLA_Obj A, FLA_Obj x, FLA_Obj x_ref,
               double *dtime, double *diff, double *gflops );


void time_Trsv( 
               int param_combo, int type, int nrepeats, int m,
               FLA_Obj A, FLA_Obj x, FLA_Obj x_ref,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    x_old;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, x, &x_old );

  FLA_Copy_external( x, x_old );


  for ( irep = 0 ; irep < nrepeats; irep++ ){

    FLA_Copy_external( x_old, x );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trsv( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
        break;
      case FLA_ALG_FRONT:
        FLA_Trsv( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
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
        REF_Trsv( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
        break;
      case FLA_ALG_FRONT:
        FLA_Trsv( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
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
        REF_Trsv( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
        break;
      case FLA_ALG_FRONT:
        FLA_Trsv( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
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
        REF_Trsv( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
        break;
      case FLA_ALG_FRONT:
        FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
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
        REF_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
        break;
      case FLA_ALG_FRONT:
        FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
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
        REF_Trsv( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
        break;
      case FLA_ALG_FRONT:
        FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
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
    FLA_Copy_external( x, x_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLA_Max_elemwise_diff( x, x_ref );
  }

  *gflops = 1.0 * 
            FLA_Obj_length( A ) * 
            FLA_Obj_width( A ) / 
            dtime_old / 
            1.0e9;

  if ( param_combo == 0 ||
       param_combo == 3 )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( x_old, x );

  FLA_Obj_free( &x_old );
}

