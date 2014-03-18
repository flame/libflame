
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Hess( FLA_Obj A, FLA_Obj t, int ilo, int ihi );
void time_Hess(
              int variant, int type, int nrepeats, int m,
              int nfc, int nlc,
              FLA_Obj A, FLA_Obj A_ref, FLA_Obj t,
              double *dtime, double *diff, double *gflops );


void time_Hess(
              int variant, int type, int nrepeats, int m,
              int nfc, int nlc,
              FLA_Obj A, FLA_Obj A_ref, FLA_Obj t,
              double *dtime, double *diff, double *gflops )
{
  int
    irep, rn, ilo, ihi;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_old;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_old );

  FLA_Copy_external( A, A_old );


  ilo = nfc + 1;
  ihi = m - nlc;

  for ( irep = 0 ; irep < nrepeats; irep++ ){
    FLA_Copy_external( A_old, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Hess( A, t, ilo, ihi );
        break;
      case FLA_ALG_FRONT:
        FLA_Hess( A, t, ilo-1, ihi-1 );
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
    FLA_Copy_external( A, A_ref );
    *diff = 0.0;
  }
  else{
    *diff = FLA_Max_elemwise_diff( A, A_ref );
  }

  rn = m - nfc - nlc;
  *gflops = 10.0 / 3.0 * rn * rn * rn /
            dtime_old / 1e9;

  *dtime = dtime_old;

  FLA_Copy_external( A_old, A );

  FLA_Obj_free( &A_old );
}

