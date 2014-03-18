#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Tridiag_UT( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t );

void time_Tridiag_UT(
                 int param_combo, int type, int nrepeats, int m,
                 FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj W,
                 double *dtime, double *diff, double *gflops );


void time_Tridiag_UT(
                 int param_combo, int type, int nrepeats, int m,
                 FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj W,
                 double *dtime, double *diff, double *gflops )
{
  int
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
        REF_Tridiag_UT( FLA_LOWER_TRIANGULAR, A, t );
        break;
      case FLA_ALG_FRONT:
        FLA_Tridiag_UT( FLA_LOWER_TRIANGULAR, A, T );
        break;
      }

      break;
    }

    case 1:
    {
      switch( type )
      {
      case FLA_ALG_REFERENCE:
        REF_Tridiag_UT( FLA_UPPER_TRIANGULAR, A, t );
        break;
      case FLA_ALG_FRONT:
        FLA_Tridiag_UT( FLA_UPPER_TRIANGULAR, A, T );
        break;
      }

      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );

  }

  {
    FLA_Obj    ATL, ATR,
               ABL, ABR;
    FLA_Obj Q, QTL, QTR,
               QBL, QBR;
    FLA_Obj    AT, AB;
    FLA_Obj    QT, QB;
    FLA_Obj E, ET, EB,
               EL, ER;
    FLA_Obj F;
    FLA_Obj eye;
    dim_t   m_A, m_T;

    m_A = FLA_Obj_length( A );
    m_T = FLA_Obj_length( T );

    FLA_Obj_create( FLA_Obj_datatype( A ), m_A, m_A, 0, 0, &Q );
    FLA_Set_to_identity( Q );

    FLA_Part_2x2( A,    &ATL, &ATR,
                        &ABL, &ABR,     1, 1, FLA_TR );
    FLA_Part_2x2( Q,    &QTL, &QTR,
                        &QBL, &QBR,     1, 1, FLA_TL );
    FLA_Part_2x1( A,    &AT,
                        &AB,    1, FLA_TOP );
    FLA_Part_2x1( Q,    &QT,
                        &QB,    1, FLA_TOP );

    // NOTE: will need to partition A, Q from left to right for upper triangular
    // case!

    if ( type == FLA_ALG_REFERENCE )
    {
//FLA_Obj_show( "A", A, "%10.3e ", "" );
//FLA_Obj_show( "t", t, "%10.3e ", "" );
      if ( FLA_Obj_is_real( A ) )
        FLA_Apply_Q_blk_external( FLA_LEFT, FLA_TRANSPOSE, FLA_COLUMNWISE, AB, t, QB );
      else
        FLA_Apply_Q_blk_external( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_COLUMNWISE, AB, t, QB );
//FLA_Obj_show( "Q", Q, "%10.3e", "" );
//FLA_Obj_show( "QBR", QBR, "%10.3e", "" );
    }
    else
    {
//FLA_Obj_show( "A", A, "%10.3e ", "" );
//FLA_Obj_show( "T", T, "%10.3e ", "" );
      FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, AB, T, W, QB );
//FLA_Obj_show( "Q", Q, "%10.3e", "" );
//FLA_Obj_show( "QB", QB, "%10.3e", "" );
    }

/*
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &eye );
    FLA_Set_to_identity( eye );
    FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_ONE, Q, Q, FLA_MINUS_ONE, eye );
    FLA_Norm_frob( eye, norm );
    FLA_Obj_extract_real_scalar( norm, diff );
    FLA_Obj_free( &eye );
*/
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &E );
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &F );

    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_ONE, A_save, Q, FLA_ZERO, E );
    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_ONE, Q, E, FLA_ZERO, F );

    FLA_Copy( A, E );

    if ( param_combo == 0 )
    {
      FLA_Part_2x1( E,    &ET,
                          &EB,    1, FLA_TOP );
      FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, EB );
      FLA_Hermitianize( FLA_LOWER_TRIANGULAR, E );
    }
    else if ( param_combo == 1 )
    {
      FLA_Part_1x2( E,    &EL, &ER,    1, FLA_LEFT );
      FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, ER );
      FLA_Hermitianize( FLA_UPPER_TRIANGULAR, E );
    }

//FLA_Obj_show( "R", E, "%10.3e + %10.3e ", "" );
//FLA_Obj_show( "Q' A Q", F, "%10.3e + %10.3e ", "" );
//FLA_Obj_show( "R", E, "%10.3e", "" );
//FLA_Obj_show( "Q' A Q", F, "%10.3e", "" );

    *diff = FLA_Max_elemwise_diff( E, F );
    FLA_Obj_free( &E );
    FLA_Obj_free( &F );
    FLA_Obj_free( &Q );
  }

  *gflops = ( 4.0 / 3.0 * m * m * m ) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &norm );
}

