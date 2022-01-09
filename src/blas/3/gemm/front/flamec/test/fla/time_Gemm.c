/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Gemm( FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C );
void time_Gemm(
               integer param_combo, integer type, integer nrepeatm, integer m, integer k, integer n,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Gemm( 
               integer param_combo, integer type, integer nrepeats, integer m, integer k, integer n,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    C_old;

  if ( param_combo != 4 )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );

  FLA_Copy_external( C, C_old );


  for ( irep = 0 ; irep < nrepeats; irep++ ){
    FLA_Copy_external( C_old, C );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Gemm( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemm( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
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
        REF_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
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
        REF_Gemm( FLA_CONJ_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemm( FLA_CONJ_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
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
        REF_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
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
        REF_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        //FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        //FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ONE, C );
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
        REF_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
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
        REF_Gemm( FLA_TRANSPOSE, FLA_CONJ_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemm( FLA_TRANSPOSE, FLA_CONJ_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
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
        REF_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
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
        REF_Gemm( FLA_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Gemm( FLA_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
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


/*
  if ( type == FLA_ALG_REFERENCE )
  {
    FLA_Copy_external( C, C_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLA_Max_elemwise_diff( C, C_ref );
  }
*/

  *gflops = 2.0 * m * k * n / 
            dtime_old / 
            1.0e9;

  if ( param_combo == 0 ||
       param_combo == 1 ||
       param_combo == 2 ||
       param_combo == 3 ||
       param_combo == 6 )
  *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

