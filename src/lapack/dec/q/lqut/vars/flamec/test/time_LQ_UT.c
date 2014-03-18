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

FLA_Error REF_LQ_UT( FLA_Obj A, FLA_Obj t );
void time_LQ(
               int variant, int type, int n_repeats, int m, int n, int nb_alg,
               FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj TT, FLA_Obj w, FLA_Obj W, FLA_Obj WW, FLA_Obj b, FLA_Obj b_ref,
               double *dtime, double *diff, double *gflops );


void time_LQ(
               int variant, int type, int n_repeats, int m, int n, int nb_alg,
               FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj TT, FLA_Obj w, FLA_Obj W, FLA_Obj WW, FLA_Obj b, FLA_Obj b_ref,
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
  fla_lqut_t*      cntl_lqut_unb;
  fla_lqut_t*      cntl_lqut_var1;
  fla_lqut_t*      cntl_lqut_var2;

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
  cntl_lqut_unb    = FLA_Cntl_lqut_obj_create( FLA_FLAT,
                                               FLA_UNB_OPT_VARIANT1,
                                               NULL,
                                               NULL,
                                               NULL );
  cntl_lqut_var1   = FLA_Cntl_lqut_obj_create( FLA_FLAT,
                                               FLA_BLOCKED_VARIANT1,
                                               bp_la,
                                               cntl_lqut_unb,
                                               cntl_apqut );
  cntl_lqut_var2   = FLA_Cntl_lqut_obj_create( FLA_FLAT,
                                               FLA_BLOCKED_VARIANT2,
                                               bp_la,
                                               cntl_lqut_unb,
                                               cntl_apqut );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b, &b_save );

  if ( FLA_Obj_is_single_precision( A ) )
    FLA_Obj_create( FLA_FLOAT, 1, 1, 0, 0, &norm );
  else
    FLA_Obj_create( FLA_DOUBLE, 1, 1, 0, 0, &norm );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( b, b_save );


  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant )
    {

    case 0:
      REF_LQ_UT( A, t );
      break;

    case 1:
    {
      // Time variant 1 
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_LQ_UT_unb_var1( A, T );
        break;
      case FLA_ALG_UNB_OPT1:
        FLA_LQ_UT_opt_var1( A, T );
        break;
      case FLA_ALG_BLOCKED:
        FLA_LQ_UT_internal( A, T, cntl_lqut_var1 );
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
        FLA_LQ_UT_unb_var2( A, t );
        break;
      case FLA_ALG_UNB_OPT1:
        FLA_LQ_UT_opt_var2( A, t );
        break;
      case FLA_ALG_BLOCKED:
        FLA_LQ_UT_internal( A, TT, cntl_lqut_var2 );
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
  FLA_Cntl_obj_free( cntl_lqut_unb );
  FLA_Cntl_obj_free( cntl_lqut_var1 );
  FLA_Cntl_obj_free( cntl_lqut_var2 );
  FLA_Blocksize_free( bp_sm );
  FLA_Blocksize_free( bp_la );

  if ( variant == 0 )
  {
    FLA_Copy_external( b, b_ref );
    FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, FLA_ONE, A, b );

    if ( FLA_Obj_is_real( A ) )
      FLA_Apply_Q_blk_external( FLA_LEFT, FLA_TRANSPOSE, FLA_ROWWISE, A, t, b );
    else
      FLA_Apply_Q_blk_external( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_ROWWISE, A, t, b );

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

    if ( variant == 1 )
      FLA_LQ_UT_solve( A, T, b_ref, b );
    else if ( variant == 2 && type == FLA_ALG_BLOCKED )
      FLA_LQ_UT_solve( A, TT, b_ref, b );
    else if ( variant == 2 && type != FLA_ALG_BLOCKED )
      FLA_LQ_UT_solve( A, t, b_ref, b );

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A_save, b, FLA_ONE, b_ref );
    FLA_Nrm2_external( b_ref, norm );
    FLA_Obj_extract_real_scalar( norm, diff );
  }

/*
if ( type == FLA_ALG_UNBLOCKED )
{
FLA_Obj AT, T2;

FLA_Obj_show( "A_orig", A_save, "%10.3e + %10.3e ", "" );

FLA_Obj_show( "LQ(A)", A, "%10.3e + %10.3e ", "" );
FLA_Obj_show( "T    ", T, "%10.3e + %10.3e ", "" );

FLA_Obj_create_copy_of( FLA_CONJ_TRANSPOSE, A_save, &AT );
FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, T, &T2 );

FLA_QR_UT_unb_var1( AT, T2 );
FLA_Transpose( AT );

FLA_Obj_show( "(QR(A^T))^T", AT, "%10.3e + %10.3e ", "" );
FLA_Obj_show( "T    ", T2, "%10.3e + %10.3e ", "" );

FLA_Obj_free( &AT );
FLA_Obj_free( &T2 );
}
*/

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

