/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define FLA_ALG_REFERENCE     0
#define FLA_ALG_UNBLOCKED     1
#define FLA_ALG_UNB_OPT       2
#define FLA_ALG_UNB_OPT_FUSED 3
#define FLA_ALG_BLOCKED       4
#define FLA_ALG_BLOCKED_FUSED 5

FLA_Error REF_Bidiag_UT_u( FLA_Obj A, FLA_Obj tu, FLA_Obj tv );

void time_Bidiag_UT_u(
               integer variant, integer type, integer n_repeats, integer m, integer n, integer nb_alg,
               FLA_Obj A, FLA_Obj TU, FLA_Obj TV, FLA_Obj TTU, FLA_Obj TTV, FLA_Obj tu, FLA_Obj tv,
               double *dtime, double *diff, double *gflops );


void time_Bidiag_UT_u(
               integer variant, integer type, integer n_repeats, integer m, integer n, integer nb_alg,
               FLA_Obj A, FLA_Obj TU, FLA_Obj TV, FLA_Obj TTU, FLA_Obj TTV, FLA_Obj tu, FLA_Obj tv,
               double *dtime, double *diff, double *gflops )
{
  integer irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, norm;

  if (
       ( variant == 1 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       ( variant == 5 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       ( variant == 1 && type == FLA_ALG_BLOCKED_FUSED ) ||
       ( variant == 5 && type == FLA_ALG_BLOCKED_FUSED ) ||
       FALSE
     )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }

  if (
       //( variant == 0 ) ||

       //( variant == 1 && type == FLA_ALG_UNBLOCKED ) ||
       //( variant == 2 && type == FLA_ALG_UNBLOCKED ) ||
       //( variant == 3 && type == FLA_ALG_UNBLOCKED ) ||
       //( variant == 4 && type == FLA_ALG_UNBLOCKED ) ||
       //( variant == 5 && type == FLA_ALG_UNBLOCKED ) ||

       //( variant == 1 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 2 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 3 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 4 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 5 && type == FLA_ALG_UNB_OPT ) ||

       //( variant == 2 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       //( variant == 3 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       //( variant == 4 && type == FLA_ALG_UNB_OPT_FUSED ) ||

       //( variant == 1 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 2 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 3 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 4 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 5 && type == FLA_ALG_BLOCKED ) ||

       //( variant == 2 && type == FLA_ALG_BLOCKED_FUSED ) ||
       //( variant == 3 && type == FLA_ALG_BLOCKED_FUSED ) ||
       //( variant == 4 && type == FLA_ALG_BLOCKED_FUSED ) ||
       FALSE
  )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );

  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

  FLA_Copy_external( A, A_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:
      REF_Bidiag_UT_u( A, tu, tv );
      break;

    // Time variant 1
    case 1:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Bidiag_UT_u_unb_var1( A, TU, TV );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Bidiag_UT_u_opt_var1( A, TU, TV );
        break;
      case FLA_ALG_UNB_OPT_FUSED:
        break;
      case FLA_ALG_BLOCKED:
        FLA_Bidiag_UT_u_blk_var1( A, TTU, TTV );
        break;
      case FLA_ALG_BLOCKED_FUSED:
        break;
      }
      break;
    }

    // Time variant 2
    case 2:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Bidiag_UT_u_unb_var2( A, TU, TV );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Bidiag_UT_u_opt_var2( A, TU, TV );
        break;
      case FLA_ALG_UNB_OPT_FUSED:
        FLA_Bidiag_UT_u_ofu_var2( A, TU, TV );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Bidiag_UT_u_blk_var2( A, TTU, TTV );
        break;
      case FLA_ALG_BLOCKED_FUSED:
        FLA_Bidiag_UT_u_blf_var2( A, TTU, TTV );
        break;
      }
      break;
    }

    // Time variant 3
    case 3:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Bidiag_UT_u_unb_var3( A, TU, TV );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Bidiag_UT_u_opt_var3( A, TU, TV );
        break;
      case FLA_ALG_UNB_OPT_FUSED:
        FLA_Bidiag_UT_u_ofu_var3( A, TU, TV );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Bidiag_UT_u_blk_var3( A, TTU, TTV );
        break;
      case FLA_ALG_BLOCKED_FUSED:
        FLA_Bidiag_UT_u_blf_var3( A, TTU, TTV );
        break;
      }
      break;
    }

    // Time variant 4
    case 4:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Bidiag_UT_u_unb_var4( A, TU, TV );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Bidiag_UT_u_opt_var4( A, TU, TV );
        break;
      case FLA_ALG_UNB_OPT_FUSED:
        FLA_Bidiag_UT_u_ofu_var4( A, TU, TV );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Bidiag_UT_u_blk_var4( A, TTU, TTV );
        break;
      case FLA_ALG_BLOCKED_FUSED:
        FLA_Bidiag_UT_u_blf_var4( A, TTU, TTV );
        break;
      }
      break;
    }

    // Time variant 5
    case 5:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Bidiag_UT_u_unb_var5( A, TU, TV );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Bidiag_UT_u_opt_var5( A, TU, TV );
        break;
      case FLA_ALG_UNB_OPT_FUSED:
        break;
      case FLA_ALG_BLOCKED:
        FLA_Bidiag_UT_u_blk_var5( A, TTU, TTV );
        break;
      case FLA_ALG_BLOCKED_FUSED:
        break;
      }
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = fla_min( *dtime, dtime_old );

  }

  if ( variant > 0 )
  //if ( FALSE )
  {

//FLA_Obj_show( "A_save", A_save, "%10.3e + %10.3e ", "" );

/*
if ( ( variant >= 3 && type == FLA_ALG_UNBLOCKED && m < 10 ) )
{
if ( FLA_Obj_is_complex( A ) )
{
FLA_Obj_show( "A", A, "%10.3e + %10.3e ", "" );
FLA_Obj_show( "TU", TU, "%10.3e + %10.3e ", "" );
FLA_Obj_show( "TV", TV, "%10.3e + %10.3e ", "" );
}
else
{
FLA_Obj_show( "A", A, "%10.3e", "" );
FLA_Obj_show( "TU", TU, "%10.3e", "" );
FLA_Obj_show( "TV", TV, "%10.3e", "" );
}
}
*/
    FLA_Obj AL, AR;
    FLA_Obj QU;
    FLA_Obj QV, QVL, QVR;
    FLA_Obj E, EL, ER;
    FLA_Obj F;
    FLA_Obj WU, WWU, WV, WWV;
    dim_t   m_A, n_A, m_TU, m_TTU;

//FLA_Obj_show( "A_save", A_save, "%10.3e", "" );

    m_A = FLA_Obj_length( A );
    n_A = FLA_Obj_width( A );
    m_TU = FLA_Obj_length( TU );
    m_TTU = FLA_Obj_length( TTU );

//FLA_Obj_show( "A bidiag", A, "%10.3e", "" );
//FLA_Copy( A_save, A );
//FLA_QR_UT( A, TU );
//FLA_Obj_show( "A QR", A, "%10.3e", "" );

    FLA_Obj_create( FLA_Obj_datatype( A ), m_A,  m_A, 0, 0, &QU );
    FLA_Obj_create( FLA_Obj_datatype( A ), n_A,  n_A, 0, 0, &QV );

    FLA_Obj_create( FLA_Obj_datatype( A ), m_TU,  m_A, 0, 0, &WU );
    FLA_Obj_create( FLA_Obj_datatype( A ), m_TTU, m_A, 0, 0, &WWU );

    FLA_Obj_create( FLA_Obj_datatype( A ), m_TU,  n_A, 0, 0, &WV );
    FLA_Obj_create( FLA_Obj_datatype( A ), m_TTU, n_A, 0, 0, &WWV );

    FLA_Set_to_identity( QU );
    FLA_Set_to_identity( QV );

    FLA_Part_1x2( QV,   &QVL, &QVR,   1, FLA_LEFT );
    FLA_Part_1x2( A,    &AL,  &AR,    1, FLA_LEFT );

    if ( type == FLA_ALG_BLOCKED || type == FLA_ALG_BLOCKED_FUSED )
    {
      FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, TTU, WWU, QU );
      FLA_Apply_Q_UT( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE, AR, TTV, WWV, QVR );
    }
    else
    {
      FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, TU, WU, QU );
      FLA_Apply_Q_UT( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE, AR, TV, WV, QVR );
    }
/*
    FLA_Obj eye;

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
    FLA_Obj_free( &WWU );
    FLA_Obj_free( &WV );
    FLA_Obj_free( &WWV );
  }
  else
  {
    *diff = 0.0;
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

