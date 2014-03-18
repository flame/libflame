

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


void time_Apply_Q_UT_inc(
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj TW, FLA_Obj W1, FLA_Obj B, FLA_Obj B_ref,
               FLA_Obj A_flat, FLA_Obj T_flat, FLA_Obj W_flat, FLA_Obj B_flat,
               double *dtime, double *diff, double *gflops );


void time_Apply_Q_UT_inc(
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj TW, FLA_Obj W1, FLA_Obj B, FLA_Obj B_ref,
               FLA_Obj A_flat, FLA_Obj T_flat, FLA_Obj W_flat, FLA_Obj B_flat,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    B_save,
    B_flat_save;

  if ( type == FLA_ALG_REFERENCE )
  {
    FLA_QR_UT( A_flat, T_flat );
  }
  else
  {
    FLASH_QR_UT_inc( A, TW );
  }

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B_flat, &B_flat_save );
  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, B, &B_save );

  FLA_Copy( B_flat, B_flat_save );
  FLASH_Copy( B, B_save );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLA_Copy( B_flat_save, B_flat );
    FLASH_Copy( B_save, B );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                        A_flat, T_flat, W_flat, B_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Apply_Q_UT_inc( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                              A, TW, W1, B );
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
    FLA_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, A_flat, B_flat );
    FLASH_Obj_hierarchify( B_flat, B_ref );

    *diff = 0.0;
  }
  else
  {
    FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                FLA_ONE, A, B );

    *diff = FLASH_Max_elemwise_diff( B, B_ref );
  }

  *gflops = 2.0 * m * m * n /
            dtime_old / 
            1.0e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy( B_flat_save, B_flat );
  FLASH_Copy( B_save, B );

  FLA_Obj_free( &B_flat_save );
  FLASH_Obj_free( &B_save );
}

