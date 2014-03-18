
#include "FLAME.h"


#define FLA_ALG_REFERENCE 0
#define FLA_ALG_UNBLOCKED 1
#define FLA_ALG_BLOCKED   2
#define FLA_ALG_RECURSIVE 3
#define FLA_ALG_OPTIMIZED 4


void time_Gemm_pp_nn(
		     int variant, int type, int nrepeats, int n, int nb_alg,
		     FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
		     double *dtime, double *diff, double *mflops );


void time_Gemm_pp_nn( 
		     int variant, int type, int nrepeats, int n, int nb_alg,
		     FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
		     double *dtime, double *diff, double *mflops )
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
      REF_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
		ONE, A, B, FLA_ONE, C );
      break;

    case 1:{
      // Time variant 1
      switch( type ){
      case FLA_ALG_UNBLOCKED:
	FLA_Gemm_pp_nn_var1( FLA_ONE, A, B, C, nb_alg );
	break;
      case FLA_ALG_BLOCKED:
        REF_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
		  ONE, A, B, FLA_ONE, C );
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

  *mflops = 2.0 * 
            FLA_Obj_length( C ) * 
            FLA_Obj_width( C ) * 
            FLA_Obj_width( A ) / 
            dtime_old / 
            1000000;

  *dtime = dtime_old;

  FLA_Copy_external( Cold, C );

  FLA_Obj_free( &Cold );
}

