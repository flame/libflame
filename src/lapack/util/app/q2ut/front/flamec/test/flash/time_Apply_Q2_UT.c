/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


void time_Apply_Q2_UT(
               integer param_combo, integer type, integer nrepeats, integer m,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj D, FLA_Obj E,
               FLA_Obj T, FLA_Obj TB, FLA_Obj TD, FLA_Obj TE, FLA_Obj W,
               double *dtime, double *diff, double *gflops );


void time_Apply_Q2_UT(
               integer param_combo, integer type, integer nrepeats, integer m,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj D, FLA_Obj E,
               FLA_Obj T, FLA_Obj TB, FLA_Obj TD, FLA_Obj TE, FLA_Obj W,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj A_save;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );

  FLASH_Copy( A, A_save );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( A_save, A );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
//FLA_Obj_show( "A_flat:", A_flat, "%12.4e", "" );
        //FLA_QR_UT( A_flat, t, T_flat );
        //FLA_QR2_UT( B_flat,
        //              D_flat, T_flat );
        FLASH_QR2_UT( B,
                      D, TD );
        FLASH_Apply_Q2_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                             D, TD, W, C, E );
//FLA_Obj_show( "A_flat_after:", A_flat, "%12.4e", "" );
//FLA_Obj_show( "T_flat_after:", T_flat, "%12.4e", "" );
        break;
      case FLA_ALG_FRONT:
//FLASH_Obj_show( "A_orig:", A, "%12.4e", "" );
        FLASH_QR2_UT( B,
                      D, TD );
        FLASH_Apply_Q2_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                             D, TD, W, C, E );
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
    //FLA_Copy( A_flat, A_flat_ref );
    *diff = 0.0;
  }
  else
  {
    //FLASH_Obj_flatten( A, BD_flat );

    //*diff = FLA_Max_elemwise_diff( BD_flat, A_flat_ref );
    *diff = 0.0;
  }

  *gflops = 2.0 * m * m * (m - m/3.0) /
            dtime_old / 
            1.0e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( A_save, A );

  FLASH_Obj_free( &A_save );
}

