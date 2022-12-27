/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#include "FLAME.h"

#define FLA_ALG_FRONT     0
#define FLA_ALG_FRONT_ALT 1


void time_QR_UT(
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj TW, FLA_Obj b, FLA_Obj x,
               double *dtime, double *diff, double *gflops );


void time_QR_UT(
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj TW, FLA_Obj b, FLA_Obj x,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj A_save;

  FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_save );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( A_save, A );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_FRONT:
        FLASH_QR_UT( A, TW );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    }
	
    *dtime = FLA_Clock() - *dtime;
    dtime_old = fla_min( *dtime, dtime_old );
  }

  {
    FLA_Obj A_save_flat, x_flat, b_flat, y_flat;
    FLA_Obj norm;

    FLASH_Obj_create_flat_copy_of_hier( A_save, &A_save_flat );
    FLASH_Obj_create_flat_copy_of_hier( b, &b_flat );
    FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, x, &x_flat );
    FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, x, &y_flat );
    FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );
/*
{
FLA_Obj T_flat;
FLASH_Obj_show( "AH", A, "%10.3e", "" );
FLA_Obj_create( FLA_DOUBLE, 2, 8, 0, 0, &T_flat );
FLA_QR_UT( A_save_flat, T_flat );
FLA_Obj_show( "A_flat", A_save_flat, "%10.3e", "" );
}
*/
    FLASH_QR_UT_solve( A, TW, b, x );

    FLASH_Obj_flatten( x, x_flat );

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A_save_flat, x_flat, FLA_MINUS_ONE, b_flat );
    FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A_save_flat, b_flat, FLA_ZERO, y_flat );
    FLA_Nrm2_external( y_flat, norm );
    FLA_Obj_extract_real_scalar( norm, diff );

    FLA_Obj_free( &A_save_flat );
    FLA_Obj_free( &b_flat );
    FLA_Obj_free( &x_flat );
    FLA_Obj_free( &y_flat );
    FLA_Obj_free( &norm );
  }

  *gflops = (         2.0   * m * n * n -
              ( 2.0 / 3.0 ) * n * n * n ) /
            dtime_old / 
            1.0e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Obj_free( &A_save );
}

