

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Axpyt( FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B );
void time_Axpyt(
               int param_combo, int type, int nrepeatm, int m, int n,
               FLA_Obj A, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Axpyt( 
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops )
{
  int
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

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Axpyt( FLA_NO_TRANSPOSE, FLA_ONE, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Axpyt( FLA_NO_TRANSPOSE, FLA_ONE, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 1:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Axpyt( FLA_TRANSPOSE, FLA_ONE, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Axpyt( FLA_TRANSPOSE, FLA_ONE, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Axpyt( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 3:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Axpyt( FLA_CONJ_TRANSPOSE, FLA_ONE, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_ONE, A, C );
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
    FLA_Copy_external( C, C_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLA_Max_elemwise_diff( C, C_ref );
  }

  *gflops = 2.0 * m * n / 
            dtime_old / 
            1.0e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

