/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Apply_Q( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Obj A, FLA_Obj t, FLA_Obj B );
void time_Apply_Q(
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj B, FLA_Obj B_ref, FLA_Obj t, FLA_Obj T, FLA_Obj W,
               double *dtime, double *diff, double *gflops );


void time_Apply_Q( 
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj B, FLA_Obj B_ref, FLA_Obj t, FLA_Obj T, FLA_Obj W,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    B_save, A_flat, B_flat;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, B, &B_save );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, A, &A_flat );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, B, &B_flat );

  FLASH_Copy( B, B_save );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( B_save, B );
    FLASH_Obj_flatten( A, A_flat );
    FLASH_Obj_flatten( B, B_flat );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Apply_Q( FLA_LEFT, FLA_TRANSPOSE, FLA_COLUMNWISE, A_flat, t, B_flat );
        break;
      case FLA_ALG_FRONT:
        FLASH_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, T, W, B );
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
    FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, FLA_ONE, A_flat, B_flat );

    FLASH_Obj_hierarchify( B_flat, B_ref );

    *diff = 0.0;
  }
  else
  {
    FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                FLA_ONE, A, B );

    *diff = FLASH_Max_elemwise_diff( B, B_ref );
  }

  *gflops = 1.0 * 
            FLASH_Obj_scalar_length( B ) * 
            FLASH_Obj_scalar_width( B ) * 
            FLASH_Obj_scalar_width( A ) / 
            dtime_old / 
            1.0e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( B_save, B );

  FLASH_Obj_free( &B_save );
  FLASH_Obj_free( &A_flat );
  FLASH_Obj_free( &B_flat );
}

