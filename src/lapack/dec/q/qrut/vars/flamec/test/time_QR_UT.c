#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_BLOCKED   1
#define FLA_ALG_UNBLOCKED 2
#define FLA_ALG_UNB_OPT1  3

extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_trmm_t*  fla_trmm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
extern fla_copyt_t* fla_copyt_cntl_blas;
extern fla_axpyt_t* fla_axpyt_cntl_blas;

extern fla_apqut_t* fla_apqut_cntl_leaf;

FLA_Error REF_QR_UT( FLA_Obj A, FLA_Obj t );
void time_QR(
               int variant, int type, int n_repeats, int m, int n, int nb_alg,
               FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj TT, FLA_Obj w, FLA_Obj W, FLA_Obj WW, FLA_Obj b, FLA_Obj x, FLA_Obj y,
               double *dtime, double *diff, double *gflops );


void time_QR(
               int variant, int type, int n_repeats, int m, int n, int nb_alg,
               FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj TT, FLA_Obj w, FLA_Obj W, FLA_Obj WW, FLA_Obj b, FLA_Obj x, FLA_Obj y,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;
  int nb_alg_sm;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, b_save, norm;

  fla_blocksize_t* bp_la;
  fla_blocksize_t* bp_sm;
  fla_apqut_t*     cntl_apqut;
  fla_qrut_t*      cntl_qrut_unb;
  fla_qrut_t*      cntl_qrut_var1;
  fla_qrut_t*      cntl_qrut_var2;

  nb_alg_sm         = max( nb_alg/4, 1 );

  bp_la             = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  bp_sm             = FLA_Blocksize_create( nb_alg_sm, nb_alg_sm, nb_alg_sm, nb_alg_sm );


  cntl_apqut       = FLA_Cntl_apqut_obj_create( FLA_FLAT,
                                                FLA_BLOCKED_VARIANT1,
                                                bp_la,
                                                NULL,
                                                fla_trmm_cntl_blas,
                                                fla_trmm_cntl_blas,
                                                fla_gemm_cntl_blas,
                                                fla_gemm_cntl_blas,
                                                fla_trsm_cntl_blas,
                                                fla_copyt_cntl_blas,
                                                fla_axpyt_cntl_blas );

  cntl_qrut_unb    = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                               FLA_UNB_OPT_VARIANT1,
                                               NULL,
                                               NULL,
                                               NULL );
  cntl_qrut_var1   = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                               FLA_BLOCKED_VARIANT1,
                                               bp_la,
                                               cntl_qrut_unb,
                                               cntl_apqut );
  cntl_qrut_var2   = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                               FLA_BLOCKED_VARIANT2,
                                               bp_la,
                                               cntl_qrut_unb,
                                               cntl_apqut );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b, &b_save );

  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( b, b_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant )
    {

    case 0:
      REF_QR_UT( A, t );
      break;

    case 1:
    {
      // Time variant 1 
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_QR_UT_unb_var1( A, T );
        break;
      case FLA_ALG_UNB_OPT1:
        FLA_QR_UT_opt_var1( A, T );
        break;
      case FLA_ALG_BLOCKED:
        FLA_QR_UT_internal( A, T, cntl_qrut_var1 );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:
    {
      // Time variant 2 
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_QR_UT_unb_var2( A, t );
        break;
      case FLA_ALG_UNB_OPT1:
        FLA_QR_UT_opt_var2( A, t );
        break;
      case FLA_ALG_BLOCKED:
        FLA_QR_UT_internal( A, TT, cntl_qrut_var2 );
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


  FLA_Cntl_obj_free( cntl_apqut );
  FLA_Cntl_obj_free( cntl_qrut_unb );
  FLA_Cntl_obj_free( cntl_qrut_var1 );
  FLA_Cntl_obj_free( cntl_qrut_var2 );
  FLA_Blocksize_free( bp_sm );
  FLA_Blocksize_free( bp_la );


/*
  if ( variant == 0 )
  {
    FLA_Copy_external( b, b_ref );
    if ( FLA_Obj_is_real( A ) )
      FLA_Apply_Q_blk_external( FLA_LEFT, FLA_TRANSPOSE, FLA_COLUMNWISE, A, t, b );
    else
      FLA_Apply_Q_blk_external( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_COLUMNWISE, A, t, b );
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
*/
  {
//FLA_Obj_show( "A", A, "%10.3e + %10.3e ", "" );
//FLA_Obj_show( "T", T, "%10.3e + %10.3e ", "" );

/*
    FLA_Copy_external( b, x );
    FLA_Apply_Q_UT( FLA_LEFT, FLA_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, T, W, x );
    FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, FLA_ONE, A, x );
*/
    if ( variant == 1 )
      FLA_QR_UT_solve( A, T, b, x );
    else if ( variant == 2 && type == FLA_ALG_BLOCKED )
      FLA_QR_UT_solve( A, TT, b, x );
    else if ( variant == 2 && type != FLA_ALG_BLOCKED )
      FLA_QR_UT_solve( A, t, b, x );

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A_save, x, FLA_MINUS_ONE, b );
    FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A_save, b, FLA_ZERO, y );
    FLA_Nrm2_external( y, norm );
    FLA_Obj_extract_real_scalar( norm, diff );

/*
if ( variant == 2 && type == FLA_ALG_UNBLOCKED )
{
FLA_Obj eye, w2, b2, Y, YT, YB, AT, AB;

FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b_save, &b2 );
FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b_save, &Y );
FLA_Obj_create( FLA_Obj_datatype( A ), FLA_Obj_length( A ), FLA_Obj_length( A ), 0, 0, &eye );
FLA_Set_to_identity( eye );
FLA_Apply_Q_UT_create_workspace( t, eye, &w2 );
FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, t, w2, eye );
FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, eye, b_save, FLA_ZERO, b2 );
FLA_Copy( b2, Y );
  FLA_Part_2x1( Y,   &YT,
                     &YB,    FLA_Obj_width( A ), FLA_TOP );
  FLA_Part_2x1( A,   &AT,
                     &AB,    FLA_Obj_width( A ), FLA_TOP );
FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                   FLA_NONUNIT_DIAG, FLA_ONE, AT, YT );

FLA_Obj_show( "A_orig", A_save, "%10.3e ", "" );
FLA_Obj_show( "A", A, "%10.3e ", "" );
FLA_Obj_show( "x", x, "%10.3e ", "" );
FLA_Obj_show( "b", b_save, "%10.3e ", "" );
FLA_Obj_show( "Q'", eye, "%10.3e ", "" );
FLA_Obj_show( "Q'*b", b2, "%10.3e ", "" );
FLA_Obj_show( "Y", Y, "%10.3e ", "" );

FLA_Obj_free( &Y );
FLA_Obj_free( &w2 );
FLA_Obj_free( &b2 );
FLA_Obj_free( &eye );
}
*/
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

