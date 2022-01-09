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

void time_Chol(
                integer param_combo, integer type, integer nrepeats, integer m, integer n,
                FLA_Obj A, FLA_Obj A_ref, FLA_Obj T, FLA_Obj t_ref, FLA_Obj B, FLA_Obj B_ref, FLA_Obj X, FLA_Obj X_ref, FLA_Obj W,
                double *dtime, double *diff, double *gflops );


void time_Chol(
                integer param_combo, integer type, integer nrepeats, integer m, integer n,
                FLA_Obj A, FLA_Obj A_ref, FLA_Obj T, FLA_Obj t_ref, FLA_Obj B, FLA_Obj B_ref, FLA_Obj X, FLA_Obj X_ref, FLA_Obj W,
                double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    B_save;
  FLA_Obj
    normx;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B, &B_save );

  if ( FLA_Obj_is_single_precision( A ) )
    FLA_Obj_create( FLA_FLOAT, 1, 1, 0, 0, &normx );
  else
    FLA_Obj_create( FLA_DOUBLE, 1, 1, 0, 0, &normx );

  FLA_Copy_external( B, B_save );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLA_Copy_external( B_save, B );
    FLA_Copy_external( B_save, B_ref );

    *dtime = FLA_Clock();

    switch( param_combo ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        FLA_Copy_external( B_ref, X_ref );
        //REF_Apply_Q( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_COLUMNWISE, A_ref, t_ref, X_ref );
        REF_Apply_Q( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_COLUMNWISE, A_ref, t_ref, X_ref );
        FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                           FLA_NONUNIT_DIAG, FLA_ONE, A_ref, X_ref );
        break;
      case FLA_ALG_FRONT:
        FLA_Copy_external( B, X );
        FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, T, W, X );
        FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                           FLA_NONUNIT_DIAG, FLA_ONE, A, X );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 1:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        FLA_Copy_external( B_ref, X_ref );
        FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                           FLA_NONUNIT_DIAG, FLA_ONE, A_ref, X_ref );
        //REF_Apply_Q( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_ROWWISE, A_ref, t_ref, X_ref );
        REF_Apply_Q( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_ROWWISE, A_ref, t_ref, X_ref );
        break;
      case FLA_ALG_FRONT:
        FLA_Copy_external( B, X );
        FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                           FLA_NONUNIT_DIAG, FLA_ONE, A, X );
        FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE, A, T, W, X );
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
    //FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A_ref, X_ref, FLA_ONE, B_ref );
    //FLA_Nrm2_external( B_ref, normx );
    //FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, normx, 1, 1, diff, 1, 1 );

    //FLA_Obj_show( "X_ref:", X_ref, "%12.4e", "" );

    *diff = 0.0;
  }
  else
  {
    //FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A, X, FLA_ONE, B );
    //FLA_Nrm2_external( B, normx );
    //FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, normx, 1, 1, diff, 1, 1 );

    //FLA_Obj_show( "X_fla:", X, "%12.4e", "" );

    *diff = FLA_Max_elemwise_diff( X, X_ref );
  }

  *gflops = 1.0 / 3.0  * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( B ) * 
            FLA_Obj_width( B ) / 
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( B_save, B );
  FLA_Copy_external( B_save, B_ref );

  FLA_Obj_free( &B_save );
  FLA_Obj_free( &normx );
}

