/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Hess_UT( FLA_Obj A, FLA_Obj t );
void time_Hess_UT(
                 integer variant, integer type, integer nrepeats, integer m,
                 FLA_Obj A, FLA_Obj A_ref, FLA_Obj t, FLA_Obj T, FLA_Obj W,
                 double *dtime, double *diff, double *gflops );


void time_Hess_UT(
                 integer variant, integer type, integer nrepeats, integer m,
                 FLA_Obj A, FLA_Obj A_ref, FLA_Obj t, FLA_Obj T, FLA_Obj W,
                 double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, norm;


  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );

  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

  FLA_Copy_external( A, A_save );

  for ( irep = 0 ; irep < nrepeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Hess_UT( A, t );
        break;
      case FLA_ALG_FRONT:
        FLA_Hess_UT( A, T );
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

  //if ( type == FLA_ALG_REFERENCE )
  //{
  //  ;
  //}
  //else
  {
    FLA_Obj    AT, AB;
    FLA_Obj Q, QT, QB;
    FLA_Obj E, ET, EB;
    FLA_Obj F;
    dim_t   m_A, m_T;

    m_A = FLA_Obj_length( A );
    m_T = FLA_Obj_length( T );

    FLA_Obj_create( FLA_Obj_datatype( A ), m_A, m_A, 0, 0, &Q );
    FLA_Set_to_identity( Q );

    FLA_Part_2x1( Q,    &QT,
                        &QB,    1, FLA_TOP );
    FLA_Part_2x1( A,    &AT,
                        &AB,    1, FLA_TOP );

    if ( type == FLA_ALG_REFERENCE )
    {
      if ( FLA_Obj_is_real( A ) )
        FLA_Apply_Q_blk_external( FLA_LEFT, FLA_TRANSPOSE, FLA_COLUMNWISE, AB, t, QB );
      else
        FLA_Apply_Q_blk_external( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_COLUMNWISE, AB, t, QB );
    }
    else
      FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, AB, T, W, QB );

    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &E );
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &F );

    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_ONE, A_save, Q, FLA_ZERO, E );
    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_ONE, Q, E, FLA_ZERO, F );

    FLA_Copy( A, E );
    FLA_Part_2x1( E,    &ET,
                        &EB,    1, FLA_TOP );
    FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, EB );

    *diff = FLA_Max_elemwise_diff( E, F );

    FLA_Obj_free( &Q );
    FLA_Obj_free( &E );
    FLA_Obj_free( &F );
  }

  *gflops = ( 10.0 / 3.0 * m * m * m ) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &norm );
}

