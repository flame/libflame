/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Sylv( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );

void time_Sylv(
                integer param_combo, integer type, integer nrepeats, integer m, integer n,
                FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref, FLA_Obj scale,
                double *dtime, double *diff, double *gflops );


void time_Sylv(
                integer param_combo, integer type, integer nrepeats, integer m, integer n,
                FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref, FLA_Obj scale,
                double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    C_old;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );

  FLA_Copy_external( C, C_old );


  for ( irep = 0 ; irep < nrepeats; irep++ ){
    FLA_Copy_external( C_old, C );

    *dtime = FLA_Clock();

    switch( param_combo ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Sylv( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
        break;
      case FLA_ALG_FRONT:
        FLA_Sylv( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 1:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Sylv( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, isgn, A, B, C, scale );
        break;
      case FLA_ALG_FRONT:
        FLA_Sylv( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, isgn, A, B, C, scale );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Sylv( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
        break;
      case FLA_ALG_FRONT:
        FLA_Sylv( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 3:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Sylv( FLA_TRANSPOSE, FLA_TRANSPOSE, isgn, A, B, C, scale );
        break;
      case FLA_ALG_FRONT:
        FLA_Sylv( FLA_TRANSPOSE, FLA_TRANSPOSE, isgn, A, B, C, scale );
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
    FLA_Copy_external( C, C_ref );
    *diff = 0.0;
  }
  else{
    *diff = FLA_Max_elemwise_diff( C, C_ref );
  }

  *gflops = ( m * m * n + n * n * m ) / 
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( C ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

