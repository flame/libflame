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
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj x, FLA_Obj y, FLA_Obj y_ref,
               double *dtime, double *diff, double *gflops );


void time_Gemv( 
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj x, FLA_Obj y, FLA_Obj y_ref,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    y_old;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, y, &y_old );

  FLA_Copy_external( y, y_old );


  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLA_Copy_external( y_old, y );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A, x, FLA_ZERO, y );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A, x, FLA_ZERO, y );
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
        REF_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A, x, FLA_ZERO, y );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A, x, FLA_ZERO, y );
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
        REF_Gemv( FLA_TRANSPOSE, FLA_ONE, A, x, FLA_ZERO, y );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemv( FLA_TRANSPOSE, FLA_ONE, A, x, FLA_ZERO, y );
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
    FLA_Copy_external( y, y_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLA_Max_elemwise_diff( y, y_ref );
  }

  *gflops = 2.0 * m * n / 
            dtime_old / 
            1.0e9;

  if ( param_combo == 0 )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( y_old, y );

  FLA_Obj_free( &y_old );
}

