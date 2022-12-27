/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_QR_UT( FLA_Obj A, FLA_Obj t );
void time_QR_UT(
                 integer variant, integer type, integer nrepeats, integer m, integer n,
                 FLA_Obj A, FLA_Obj A_ref, FLA_Obj t, FLA_Obj T, FLA_Obj W, FLA_Obj b, FLA_Obj b_orig,
                 double *dtime, double *diff, double *gflops );


void time_QR_UT(
                 integer variant, integer type, integer nrepeats, integer m, integer n,
                 FLA_Obj A, FLA_Obj A_ref, FLA_Obj t, FLA_Obj T, FLA_Obj W, FLA_Obj b, FLA_Obj b_orig,
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
        REF_QR_UT( A, t );
        break;
      case FLA_ALG_FRONT:
        FLA_QR_UT( A, T );
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

  if ( type == FLA_ALG_REFERENCE )
  {
    FLA_Obj AT, AB;
    FLA_Obj bT, bB;
    FLA_Obj y;

    FLA_Obj_create( FLA_Obj_datatype( b ), n, 1, 0, 0, &y );

    FLA_Copy_external( b, b_orig );

    if ( FLA_Obj_is_real( A ) )
      FLA_Apply_Q_blk_external( FLA_LEFT, FLA_TRANSPOSE, FLA_COLUMNWISE, A, t, b );
    else
      FLA_Apply_Q_blk_external( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_COLUMNWISE, A, t, b );

    FLA_Part_2x1( A,    &AT,
                        &AB,    FLA_Obj_width( A ), FLA_TOP );
    FLA_Part_2x1( b,    &bT,
                        &bB,    FLA_Obj_width( A ), FLA_TOP );
    FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, FLA_ONE, AT, bT );
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A_save, bT, FLA_ONE, b_orig );
    FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A_save, b_orig, FLA_ZERO, y );
    FLA_Nrm2_external( y, norm );
    FLA_Obj_extract_real_scalar( norm, diff );

    FLA_Obj_free( &y );
  }
  else
  {
    FLA_Obj x, y;

    FLA_Obj_create( FLA_Obj_datatype( b ), n, 1, 0, 0, &y );
    FLA_Obj_create( FLA_Obj_datatype( b ), n, 1, 0, 0, &x );

    FLA_Copy_external( b, b_orig );

    FLA_QR_UT_solve( A, T, b, x );

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A_save, x, FLA_ONE, b_orig );
    FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A_save, b_orig, FLA_ZERO, y );
    FLA_Nrm2_external( y, norm );
    FLA_Obj_extract_real_scalar( norm, diff );

    FLA_Obj_free( &x );
    FLA_Obj_free( &y );
  }

  *gflops = (         2.0   * m * n * n -
              ( 2.0 / 3.0 ) * n * n * n ) /
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

