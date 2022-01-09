/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Bidiag_UT( FLA_Obj A, FLA_Obj tu, FLA_Obj tv );

void time_Bidiag_UT(
                 integer param_combo, integer type, integer nrepeats, integer m, integer n,
                 FLA_Obj A, FLA_Obj tu, FLA_Obj tv, FLA_Obj TU, FLA_Obj TV,
                 double *dtime, double *diff, double *gflops );


void time_Bidiag_UT(
                 integer param_combo, integer type, integer nrepeats, integer m, integer n,
                 FLA_Obj A, FLA_Obj tu, FLA_Obj tv, FLA_Obj TU, FLA_Obj TV,
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

    switch( param_combo ){

    case 0:
    {
      switch( type )
      {
      case FLA_ALG_REFERENCE:
        REF_Bidiag_UT( A, tu, tv );
        break;
      case FLA_ALG_FRONT:
        FLA_Bidiag_UT( A, TU, TV );
        break;
      }

      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );

  }

  {
    FLA_Obj AL, AR;
    FLA_Obj ATL, ATR,
            ABL, ABR;
    FLA_Obj QU;
    FLA_Obj QV, QVL, QVR;
    FLA_Obj E, EL, ER;
    FLA_Obj F;
    FLA_Obj WU, WV, eye;
    FLA_Obj tvT,
            tvB;
    dim_t   m_A, n_A, m_TU;


//FLA_Obj_show( "A_save", A_save, "%10.3e", "" );

    m_A = FLA_Obj_length( A );
    n_A = FLA_Obj_width( A );
    m_TU = FLA_Obj_length( TU );

    FLA_Obj_create( FLA_Obj_datatype( A ), m_A,  m_A, 0, 0, &QU );
    FLA_Obj_create( FLA_Obj_datatype( A ), n_A,  n_A, 0, 0, &QV );

    FLA_Obj_create( FLA_Obj_datatype( A ), m_TU,  m_A, 0, 0, &WU );
    FLA_Obj_create( FLA_Obj_datatype( A ), m_TU,  n_A, 0, 0, &WV );

    FLA_Set_to_identity( QU );
    FLA_Set_to_identity( QV );

    FLA_Part_1x2( QV,   &QVL, &QVR,   1, FLA_LEFT );
    FLA_Part_1x2( A,    &AL,  &AR,    1, FLA_LEFT );
    FLA_Part_2x2( A,    &ATL, &ATR,
                        &ABL, &ABR,   1, 1, FLA_BL );
    FLA_Part_2x1( tv,   &tvT,
                        &tvB,   1, FLA_BOTTOM );

    if ( type == FLA_ALG_REFERENCE )
    {
      if ( FLA_Obj_is_real( A ) )
        FLA_Apply_Q_blk_external( FLA_LEFT, FLA_TRANSPOSE, FLA_COLUMNWISE, A, tu, QU );
      else
        FLA_Apply_Q_blk_external( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_COLUMNWISE, A, tu, QU );
      //FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_ROWWISE, AR, tv, QVR );
      //
      // Need to apply backwards transformation, since vectors are stored columnwise.
      // QL? RQ?
      //
      //FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_ROWWISE, ATR, tvT, QVR );
      //FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_CONJ_TRANSPOSE, FLA_ROWWISE, ATR, tvT, QVR );
      //FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_CONJ_TRANSPOSE, FLA_ROWWISE, AR, tvT, QVR );
    }
    else
    {
      FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, TU, WU, QU );
      FLA_Apply_Q_UT( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE, AR, TV, WV, QVR );
    }

/*
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &eye );     
    FLA_Set_to_identity( eye );

    //FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
    //          FLA_ONE, QV, QV, FLA_MINUS_ONE, eye );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, QU, QU, FLA_MINUS_ONE, eye );

FLA_Obj_show( "eye", eye, "%10.3e", "" );
    FLA_Norm_frob( eye, norm );
    FLA_Obj_extract_real_scalar( norm, diff );
    FLA_Obj_free( &eye );
*/

    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &E );
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &F );

    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, A_save, QV, FLA_ZERO, E );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, QU, E, FLA_ZERO, F );

//FLA_Obj_show( "A_save", A_save, "%10.3e", "" );

    FLA_Copy( A, E );
    FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, E );
    FLA_Part_1x2( E,    &EL, &ER,      1, FLA_LEFT );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, ER );

//FLA_Obj_show( "B", E, "%10.3e", "" );
//FLA_Obj_show( "Q'AV", F, "%10.3e", "" );
//FLA_Obj_show( "B", E, "%10.3e + %10.3e ", "" );
//FLA_Obj_show( "Q'AV", F, "%10.3e + %10.3e ", "" );

    *diff = FLA_Max_elemwise_diff( E, F );
    FLA_Obj_free( &E );
    FLA_Obj_free( &F );

    FLA_Obj_free( &QU );
    FLA_Obj_free( &QV );
    FLA_Obj_free( &WU );
    FLA_Obj_free( &WV );
  }

  *gflops = 4.0 * n * n * ( m - n / 3.0 ) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &norm );
}

