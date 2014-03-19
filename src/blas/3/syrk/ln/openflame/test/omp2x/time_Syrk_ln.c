/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_OPENMP_1TASK     5
#define FLA_ALG_OPENMP_2TASKS    6
#define FLA_ALG_OPENMP_2LOOPS    7
#define FLA_ALG_OPENMP_2LOOPSPLUS 8




void time_Syrk_ln(
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Syrk_ln( 
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old;

  FLA_Obj
    C_old;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );

  FLA_Copy_external( C, C_old );

  for ( irep = 0 ; irep < nrepeats; irep++ ){

    FLA_Copy_external( C_old, C );

    *dtime = FLA_Clock();

    switch( variant ){
    case 0:
      // Time reference implementation
      REF_Syrk_ln( FLA_ONE, A, FLA_ONE, C );
      break;

    case 1:{
      // Time variant 1
      switch( type ){
      case FLA_ALG_OPENMP_1TASK:
        FLA_Syrk_ln_omp1t_var1( A, C );
        break;
      case FLA_ALG_OPENMP_2TASKS:
        FLA_Syrk_ln_omp2t_var1( A, C );
        break;
      case FLA_ALG_OPENMP_2LOOPS:
        FLA_Syrk_ln_omp2l_var1( A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 2:{
      // Time variant 2
      switch( type ){
      case FLA_ALG_OPENMP_1TASK:
        FLA_Syrk_ln_omp1t_var2( A, C );
        break;
      case FLA_ALG_OPENMP_2TASKS:
        FLA_Syrk_ln_omp2t_var2( A, C );
        break;
      case FLA_ALG_OPENMP_2LOOPS:
        FLA_Syrk_ln_omp2l_var2( A, C );
        break;
      case FLA_ALG_OPENMP_2LOOPSPLUS:
        FLA_Syrk_ln_omp2x_var2( A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    } 
    case 3:{
      // Time variant 3 
      switch( type ){
      case FLA_ALG_OPENMP_1TASK:
        FLA_Syrk_ln_omp1t_var3( A, C );
        break;
      case FLA_ALG_OPENMP_2TASKS:
        FLA_Syrk_ln_omp2t_var3( A, C );
        break;
      case FLA_ALG_OPENMP_2LOOPS:
        FLA_Syrk_ln_omp2l_var3( A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 4:{
      // Time variant 4
      switch( type ){
      case FLA_ALG_OPENMP_1TASK:
        FLA_Syrk_ln_omp1t_var4( A, C );
        break;
      case FLA_ALG_OPENMP_2TASKS:
        FLA_Syrk_ln_omp2t_var4( A, C );
        break;
      case FLA_ALG_OPENMP_2LOOPS:
        FLA_Syrk_ln_omp2l_var4( A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 5:{
      // Time variant 5
      switch( type ){
      case FLA_ALG_OPENMP_1TASK:
        FLA_Syrk_ln_omp1t_var5( A, C );
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
    FLA_Copy_external( C, C_ref );
    *diff = 0.0;
  }
  else{
    *diff = FLA_Max_elemwise_diff( C, C_ref );
    //FLA_Obj_show( "C:", C, "%f", "\n");
  }

  *gflops = 1.0 * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) * 
            FLA_Obj_width( A ) / 
            dtime_old / 
            1e9;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

