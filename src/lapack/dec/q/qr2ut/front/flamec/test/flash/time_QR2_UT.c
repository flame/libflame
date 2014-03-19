/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


void time_QR2_UT(
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj A_flat_ref, FLA_Obj B, FLA_Obj B_flat, FLA_Obj D, FLA_Obj D_flat, FLA_Obj A_flat, FLA_Obj t, FLA_Obj T, FLA_Obj T_flat,
               double *dtime, double *diff, double *gflops );


void time_QR2_UT(
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj A_flat_ref, FLA_Obj B, FLA_Obj B_flat, FLA_Obj D, FLA_Obj D_flat, FLA_Obj A_flat, FLA_Obj t, FLA_Obj T, FLA_Obj T_flat,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    BD_flat, A_save, A_flat_save;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A_flat, &A_flat_save );

  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, A, &BD_flat );

  FLASH_Copy( A, A_save );
  FLASH_Copy( A_flat, A_flat_save );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( A_save, A );
    FLA_Copy( A_flat_save, A_flat );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
//FLA_Obj_show( "A_flat:", A_flat, "%12.4e", "" );
        //FLA_QR_UT( A_flat, t, T_flat );
        FLA_QR2_UT( B_flat,
                      D_flat, T_flat );
//FLA_Obj_show( "A_flat_after:", A_flat, "%12.4e", "" );
//FLA_Obj_show( "T_flat_after:", T_flat, "%12.4e", "" );
        break;
      case FLA_ALG_FRONT:
//FLASH_Obj_show( "A_orig:", A, "%12.4e", "" );
        FLASH_QR2_UT( B,
                        D, T );
//FLASH_Obj_show( "A_after:", A, "%12.4e", "" );
//FLASH_Obj_show( "T_after:", T, "%12.4e", "" );
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
    FLA_Copy( A_flat, A_flat_ref );
    *diff = 0.0;
  }
  else
  {
    FLASH_Obj_flatten( A, BD_flat );

    *diff = FLA_Max_elemwise_diff( BD_flat, A_flat_ref );
  }

  *gflops = 2.0 * n * n * (m - n/3.0) /
            dtime_old / 
            1.0e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( A_save, A );
  FLA_Copy( A_flat_save, A_flat );

  FLASH_Obj_free( &A_save );
  FLA_Obj_free( &A_flat_save );
  FLA_Obj_free( &BD_flat );
}

