
#include "FLAME.h"

#define FLA_ALG_REFERENCE     0
#define FLA_ALG_UNBLOCKED     1
#define FLA_ALG_UNB_OPT       2
#define FLA_ALG_BLOCKED       3

FLA_Error REF_Svd_uv( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
                      double* dtime_bred, double* dtime_bsvd, double* dtime_appq,
                      double* dtime_qrfa, double* dtime_gemm );
FLA_Error REF_Svdd_uv( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
                       double* dtime_bred, double* dtime_bsvd, double* dtime_appq,
                       double* dtime_qrfa, double* dtime_gemm );
FLA_Error REF_Svd_uv_components( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
                                 double* dtime_bred, double* dtime_bsvd, double* dtime_appq,
                                 double* dtime_qrfa, double* dtime_gemm );
FLA_Error REF_Svdd_uv_components( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
                                  double* dtime_bred, double* dtime_bsvd, double* dtime_appq,
                                  double* dtime_qrfa, double* dtime_gemm );
FLA_Error FLA_Svd_uv_var1_components( dim_t n_iter_max, dim_t k_accum, dim_t b_alg,
                                      FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
                                      double* dtime_bred, double* dtime_bsvd, double* dtime_appq,
                                      double* dtime_qrfa, double* dtime_gemm );
FLA_Error FLA_Svd_uv_var2_components( dim_t n_iter_max, dim_t k_accum, dim_t b_alg,
                                      FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
                                      double* dtime_bred, double* dtime_bsvd, double* dtime_appq,
                                      double* dtime_qrfa, double* dtime_gemm );

void time_Svd_uv_components(
               int variant, int type, int n_repeats, int m, int n, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
               double* dtime, double* diff1, double* diff2, double* gflops,
               double* dtime_bred, double* gflops_bred,
               double* dtime_bsvd, double* gflops_bsvd,
               double* dtime_appq, double* gflops_appq,
               double* dtime_qrfa, double* gflops_qrfa,
               double* dtime_gemm, double* gflops_gemm, int* k_perf );


void time_Svd_uv_components(
               int variant, int type, int n_repeats, int m, int n, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
               double* dtime, double* diff1, double* diff2, double* gflops,
               double* dtime_bred, double* gflops_bred,
               double* dtime_bsvd, double* gflops_bsvd,
               double* dtime_appq, double* gflops_appq,
               double* dtime_qrfa, double* gflops_qrfa,
               double* dtime_gemm, double* gflops_gemm, int* k_perf )
{
  int     i;
  double  k;
  double  dtime_save      = 1.0e9;
  double  dtime_bred_save = 1.0e9;
  double  dtime_bsvd_save = 1.0e9;
  double  dtime_appq_save = 1.0e9;
  double  dtime_qrfa_save = 1.0e9;
  double  dtime_gemm_save = 1.0e9;
  double  flops_bred;
  double  flops_bsvd;
  double  flops_appq;
  double  flops_qrfa;
  double  flops_gemm;
  double  mult_bred;
  double  mult_bsvd;
  double  mult_appq;
  double  mult_qrfa;
  double  mult_gemm;

  FLA_Obj A_save;

  if (
       ( variant == -2 ) ||
       ( variant == -3 ) ||
       //( variant == 0 ) ||
       //( variant == -1 ) ||
       //( variant == 1 ) ||
       //( variant == 2 ) ||
       FALSE
     )
  {
    *gflops      = 0.0;
    *dtime       = 0.0;
    *diff1       = 0.0;
    *diff2       = 0.0;
    *dtime_bred  = 0.0;
    *dtime_bsvd  = 0.0;
    *dtime_appq  = 0.0;
    *dtime_qrfa  = 0.0;
    *dtime_gemm  = 0.0;
    *gflops_bred = 0.0;
    *gflops_bsvd = 0.0;
    *gflops_appq = 0.0;
    *gflops_qrfa = 0.0;
    *gflops_gemm = 0.0;
    *k_perf      = 0;
    return;
  }

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );

  FLA_Copy_external( A, A_save );

  for ( i = 0 ; i < n_repeats; i++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case -2:
    {
      *k_perf = 0;
      REF_Svd_uv( A, s, U, V,
                  dtime_bred, dtime_bsvd, dtime_appq,
                  dtime_qrfa, dtime_gemm );
      break;
    }

    case -3:
    {
      *k_perf = 0;
      REF_Svdd_uv( A, s, U, V,
                   dtime_bred, dtime_bsvd, dtime_appq,
                   dtime_qrfa, dtime_gemm );
      break;
    }

    case 0:
    {
      *k_perf = 0;
      REF_Svd_uv_components( A, s, U, V,
                             dtime_bred, dtime_bsvd, dtime_appq,
                             dtime_qrfa, dtime_gemm );
      break;
    }

    case -1:
    {
      *k_perf = 0;
      REF_Svdd_uv_components( A, s, U, V,
                              dtime_bred, dtime_bsvd, dtime_appq,
                              dtime_qrfa, dtime_gemm );
      break;
    }

    // Time variant 1
    case 1:
    {
      *k_perf = FLA_Svd_uv_var1_components( n_iter_max, k_accum, b_alg,
                                            A, s, U, V,
                                            dtime_bred, dtime_bsvd, dtime_appq,
                                            dtime_qrfa, dtime_gemm );
      break;
    }

    // Time variant 2
    case 2:
    {
      *k_perf = FLA_Svd_uv_var2_components( n_iter_max, k_accum, b_alg,
                                            A, s, U, V,
                                            dtime_bred, dtime_bsvd, dtime_appq,
                                            dtime_qrfa, dtime_gemm );
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    if ( *dtime < dtime_save )
    {
      dtime_save      = *dtime;
      dtime_bred_save = *dtime_bred;
      dtime_bsvd_save = *dtime_bsvd;
      dtime_appq_save = *dtime_appq;
      dtime_qrfa_save = *dtime_qrfa;
      dtime_gemm_save = *dtime_gemm;
    }
  }

  *dtime      = dtime_save;
  *dtime_bred = dtime_bred_save;
  *dtime_bsvd = dtime_bsvd_save;
  *dtime_appq = dtime_appq_save;
  *dtime_qrfa = dtime_qrfa_save;
  *dtime_gemm = dtime_gemm_save;

  {
    FLA_Obj USVh, eyeU, eyeV, normU, normV;
    double  diffU, diffV;

    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &USVh );
    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, U, &eyeU );
    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, V, &eyeV );
    FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &normU );
    FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &normV );

    FLA_Set_to_identity( eyeU );
    FLA_Set_to_identity( eyeV );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, U, U, FLA_MINUS_ONE, eyeU );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, V, V, FLA_MINUS_ONE, eyeV );
    FLA_Norm_frob( eyeU, normU );
    FLA_Norm_frob( eyeV, normV );
    FLA_Obj_extract_real_scalar( normU, &diffU );
    FLA_Obj_extract_real_scalar( normV, &diffV );
    *diff2 = max( diffU, diffV );

    if ( variant < 1 ) // LAPACK computes and stores V'
    {
      FLA_Obj UL, UR;
      FLA_Obj VT,
              VB;

      FLA_Part_1x2( U,   &UL, &UR,   min( m, n ), FLA_LEFT );
      FLA_Part_2x1( V,   &VT,
                         &VB,        min( m, n ), FLA_TOP );

      FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, s, UL );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                FLA_ONE, UL, VT, FLA_ZERO, USVh );
    }
    else // libflame computes and stores V
    {
      FLA_Obj UL, UR;
      FLA_Obj VL, VR;

      FLA_Part_1x2( U,   &UL, &UR,   min( m, n ), FLA_LEFT );
      FLA_Part_1x2( V,   &VL, &VR,   min( m, n ), FLA_LEFT );

      FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, s, UL );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                FLA_ONE, UL, VL, FLA_ZERO, USVh );
    }

    FLA_Axpy( FLA_MINUS_ONE, A_save, USVh );
    FLA_Norm_frob( USVh, normU );
    FLA_Obj_extract_real_scalar( normU, diff1 );

    FLA_Obj_free( &USVh );
    FLA_Obj_free( &eyeU );
    FLA_Obj_free( &eyeV );
    FLA_Obj_free( &normU );
    FLA_Obj_free( &normV );

  }

  k = 2.00;

  if ( m < 1.6 * n )
  {
    flops_bred = ( 4.0           * m * n * n -
                   4.0 / 3.0     * n * n * n );
    flops_bsvd = ( 13.0      * k * n * n     +
                    3.0      * k * m * n * n +
                    3.0      * k * n * n * n );
    flops_appq = ( 4.0           * m * m * n -
                   4.0           * m * n * n +
                   4.0 / 3.0     * n * n * n ) +
                 ( 4.0 / 3.0     * n * n * n );
    flops_qrfa = 0.0;
    flops_gemm = 0.0;
  }
  else
  {
    flops_bred = ( 4.0           * n * n * n -
                   4.0 / 3.0     * n * n * n );
    flops_bsvd = ( 13.0      * k * n * n     +
                    3.0      * k * n * n * n +
                    3.0      * k * n * n * n );
    flops_appq = ( 4.0           * m * m * n -
                   4.0           * m * n * n +
                   4.0 / 3.0     * n * n * n ) +
                 ( 4.0 / 3.0     * n * n * n ) +
                 ( 4.0 / 3.0     * n * n * n );
    flops_qrfa = ( 2.0           * m * n * n -
                   2.0 / 3.0     * n * n * n );
    flops_gemm = ( 2.0           * m * n * n );
  }


  if ( FLA_Obj_is_complex( A ) )
  {
    mult_bred = 4.0;
    mult_bsvd = 2.0;
    mult_appq = 4.0;
    mult_qrfa = 4.0;
    mult_gemm = 4.0;
  }
  else
  {
    mult_bred = 1.0;
    mult_bsvd = 1.0;
    mult_appq = 1.0;
    mult_qrfa = 1.0;
    mult_gemm = 1.0;
  }

  *gflops = ( mult_bred * flops_bred + 
              mult_bsvd * flops_bsvd + 
              mult_appq * flops_appq + 
              mult_qrfa * flops_qrfa + 
              mult_gemm * flops_gemm ) / *dtime / 1e9;

  *gflops_bred = ( mult_bred * flops_bred ) / *dtime_bred / 1e9;
  *gflops_bsvd = ( mult_bsvd * flops_bsvd ) / *dtime_bsvd / 1e9;
  *gflops_appq = ( mult_appq * flops_appq ) / *dtime_appq / 1e9;

  if ( *dtime_qrfa > 0.0 && *dtime_gemm > 0.0 )
  {
    *gflops_qrfa = ( mult_qrfa * flops_qrfa ) / *dtime_qrfa / 1e9;
    *gflops_gemm = ( mult_gemm * flops_gemm ) / *dtime_gemm / 1e9;
  }
  else
  {
    *gflops_qrfa = 0.0;
    *gflops_gemm = 0.0;
  }

  FLA_Copy_external( A_save, A );

  FLA_Obj_free( &A_save );
}

