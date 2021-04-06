/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_UNBLOCKED 1
#define FLA_ALG_UNB_OPT   2


FLA_Error REF_Accum_T_UT_fc( FLA_Obj A, FLA_Obj t );
void time_Accum_T_UT_fc(
                 integer variant, integer type, integer nrepeats, integer m, integer n,
                 FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj W, FLA_Obj b, FLA_Obj b_ref,
                 double *dtime, double *diff, double *gflops );


void time_Accum_T_UT_fc(
                 integer variant, integer type, integer nrepeats, integer m, integer n,
                 FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj W, FLA_Obj b, FLA_Obj b_ref,
                 double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, b_save, norm;


  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b, &b_save );

  if ( FLA_Obj_is_single_precision( A ) )
    FLA_Obj_create( FLA_FLOAT, 1, 1, 0, 0, &norm );
  else
    FLA_Obj_create( FLA_DOUBLE, 1, 1, 0, 0, &norm );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( b, b_save );


  for ( irep = 0 ; irep < nrepeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        FLA_QR_UT( A, T );
        break;
      case FLA_ALG_UNBLOCKED:
        FLA_QR_UT( A, T );
        FLA_QR_UT_recover_tau( T, t );
        FLA_Set( FLA_ZERO, T );
        FLA_Accum_T_UT_fc_unb_var1( A, t, T );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_QR_UT( A, T );
        FLA_QR_UT_recover_tau( T, t );
        FLA_Set( FLA_ZERO, T );
        FLA_Accum_T_UT_fc_opt_var1( A, t, T );
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
    FLA_Copy_external( b, b_ref );
    FLA_Apply_Q_UT( FLA_LEFT, FLA_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, T, W, b );
    FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, FLA_ONE, A, b );
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A_save, b, FLA_ONE, b_ref );
    FLA_Nrm2_external( b_ref, norm );
    if ( FLA_Obj_is_single_precision( A ) )
      *diff = *(FLA_FLOAT_PTR(norm));
    else
      *diff = *(FLA_DOUBLE_PTR(norm));
  }
  else
  {
    FLA_Copy_external( b, b_ref );
    FLA_Apply_Q_UT( FLA_LEFT, FLA_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, T, W, b );
    FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, FLA_ONE, A, b );
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A_save, b, FLA_ONE, b_ref );
    FLA_Nrm2_external( b_ref, norm );
    if ( FLA_Obj_is_single_precision( A ) )
      *diff = *(FLA_FLOAT_PTR(norm));
    else
      *diff = *(FLA_DOUBLE_PTR(norm));
  }

  *gflops = 2.0 * n * n * (m - n/3.0) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( b_save, b );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &b_save );
  FLA_Obj_free( &norm );
}

