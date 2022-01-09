/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE   0
#define FLA_ALG_OPENMP_BVAR 1
#define FLA_ALG_OPENMP_CVAR 2




void time_Gemm_nn(
               integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
               double *dtime, double *diff, double *gflops );


void time_Gemm_nn( 
               integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep,
    info, lwork;

  double
    dtime_old,
    d_minus_one = -1.0, d_one = 1.0;

  FLA_Obj
    Cold;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &Cold );

  FLA_Copy_external( C, Cold );

  for ( irep = 0 ; irep < nrepeats; irep++ ){
    FLA_Copy_external( Cold, C );

    *dtime = FLA_Clock();

    switch( variant ){
    case 0:
      // Time reference implementation
      REF_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                FLA_ONE, A, B, FLA_ONE, C );
      break;

    case 1:{
      // Time variant 1
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var1( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 2:{
      // Time variant 2
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var2( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 3:{
      // Time variant 3
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var3( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 4:{
      // Time variant 4
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var4( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 5:{
      // Time variant 5
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var5( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 6:{
      // Time variant 6
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var6( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 13:{
      // Time variant 1->3
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var13( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 15:{
      // Time variant 1->5
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var15( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 31:{
      // Time variant 3->1 
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var31( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 35:{
      // Time variant 3->5 
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var35( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 51:{
      // Time variant 5->1 
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var51( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 53:{
      // Time variant 5->3 
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var53( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    }

    if ( irep == 0 )
      dtime_old = FLA_Clock() - *dtime;
    else{
      *dtime = FLA_Clock() - *dtime;
      dtime_old = min( *dtime, dtime_old );
    }
  }


  if ( variant == 0 ){
    FLA_Copy_external( C, Cref );
    *diff = 0.0;
  }
  else{
    *diff = FLA_Max_elemwise_diff( C, Cref );
    //FLA_Obj_show( "C:", C, "%f", "\n");
  }

  *gflops = 2.0 * 
            FLA_Obj_length( C ) * 
            FLA_Obj_width( C ) * 
            FLA_Obj_width( A ) / 
            dtime_old / 
            1e9;

  *dtime = dtime_old;

  FLA_Copy_external( Cold, C );

  FLA_Obj_free( &Cold );
}

