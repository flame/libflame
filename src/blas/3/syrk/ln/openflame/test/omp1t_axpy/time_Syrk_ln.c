/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_OPENMP_1TASK_FLA_PIP     1
#define FLA_ALG_OPENMP_1TASK_FLA_CIR     2
#define FLA_ALG_OPENMP_1TASK_REF_PIP     3
#define FLA_ALG_OPENMP_1TASK_REF_CIR     4



void time_Syrk_ln(
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
               double *dtime, double *diff, double *gflops );


void time_Syrk_ln( 
	       int variant, int type, int nrepeats, int n, int nb_alg,
	       FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
	       double *dtime, double *diff, double *gflops )
{
  int
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
      REF_Syrk_ln( A, C );
      break;

    case 1:{
      // Time variant 1
      switch( type ){
      default:
	printf("trouble\n");
      }

      break;
    }
    case 2:{
      // Time variant 2
      switch( type ){
      default:
	printf("trouble\n");
      }

      break;
    } 
    case 3:{
      // Time variant 3 
      switch( type ){
      default:
	printf("trouble\n");
      }

      break;
    }
    case 4:{
      // Time variant 5
      switch( type ){
      default:
	printf("trouble\n");
      }

      break;
    }
    case 5:{
      // Time variant 4
      switch( type ){
      case FLA_ALG_OPENMP_1TASK_FLA_PIP:
	FLA_Syrk_ln_omp1t_var5_fp( A, C, nb_alg );
	break;
      case FLA_ALG_OPENMP_1TASK_FLA_CIR:
	FLA_Syrk_ln_omp1t_var5_fc( A, C, nb_alg );
	break;
      case FLA_ALG_OPENMP_1TASK_REF_PIP:
	FLA_Syrk_ln_omp1t_var5_rp( A, C, nb_alg );
	break;
      case FLA_ALG_OPENMP_1TASK_REF_CIR:
	FLA_Syrk_ln_omp1t_var5_rc( A, C, nb_alg );
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
  }

  *gflops = 1.0 * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) * 
            FLA_Obj_width( A ) / 
            dtime_old / 
            1e9;

  *dtime = dtime_old;

  FLA_Copy_external( Cold, C );

  FLA_Obj_free( &Cold );
}

