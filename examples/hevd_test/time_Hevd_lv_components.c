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
#define FLA_ALG_BLOCKED       3

FLA_Error REF_Hevd_lv( FLA_Obj A, FLA_Obj l,
                       double* dtime_tred, double* dtime_tevd, double* dtime_appq );
FLA_Error REF_Hevdd_lv( FLA_Obj A, FLA_Obj l,
                        double* dtime_tred, double* dtime_tevd, double* dtime_appq );
FLA_Error REF_Hevdr_lv( FLA_Obj A, FLA_Obj l, FLA_Obj Z,
                        double* dtime_tred, double* dtime_tevd, double* dtime_appq );
FLA_Error REF_Hevd_lv_components( FLA_Obj A, FLA_Obj l,
                                  double* dtime_tred, double* dtime_tevd, double* dtime_appq );
FLA_Error REF_Hevdd_lv_components( FLA_Obj A, FLA_Obj l,
                                   double* dtime_tred, double* dtime_tevd, double* dtime_appq );
FLA_Error REF_Hevdr_lv_components( FLA_Obj A, FLA_Obj l, FLA_Obj Z,
                                   double* dtime_tred, double* dtime_tevd, double* dtime_appq );
FLA_Error FLA_Hevd_lv_var1_components( dim_t n_iter_max, FLA_Obj A, FLA_Obj l, dim_t k_accum, dim_t b_alg,
                                       double* dtime_tred, double* dtime_tevd, double* dtime_appq );
FLA_Error FLA_Hevd_lv_var2_components( dim_t n_iter_max, FLA_Obj A, FLA_Obj l, dim_t k_accum, dim_t b_alg,
                                       double* dtime_tred, double* dtime_tevd, double* dtime_appq );
FLA_Error FLA_Hevd_lv_var3_components( dim_t n_iter_max, FLA_Obj A, FLA_Obj l, dim_t k_accum, dim_t b_alg,
                                       double* dtime_tred, double* dtime_tevd, double* dtime_appq );
FLA_Error FLA_Hevd_lv_var4_components( dim_t n_iter_max, FLA_Obj A, FLA_Obj l, dim_t k_accum, dim_t b_alg,
                                       double* dtime_tred, double* dtime_tevd, double* dtime_appq );

void time_Hevd_lv_components(
               int variant, int type, int n_repeats, int m, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj l,
               double* dtime, double* diff1, double* diff2, double* gflops,
               double* dtime_tred, double* gflops_tred,
               double* dtime_tevd, double* gflops_tevd,
               double* dtime_appq, double* gflops_appq, int* k_perf );


void time_Hevd_lv_components(
               int variant, int type, int n_repeats, int m, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj l,
               double* dtime, double* diff1, double* diff2, double* gflops,
               double* dtime_tred, double* gflops_tred,
               double* dtime_tevd, double* gflops_tevd,
               double* dtime_appq, double* gflops_appq, int* k_perf )
{
  int     i;
  double  k;
  double  dtime_save      = 1.0e9;
  double  dtime_tred_save = 1.0e9;
  double  dtime_tevd_save = 1.0e9;
  double  dtime_appq_save = 1.0e9;
  double  flops_tred;
  double  flops_tevd;
  double  flops_appq;
  double  mult_tred;
  double  mult_tevd;
  double  mult_appq;

  FLA_Obj A_save, Z;

  if (
       ( variant == -3 ) ||
       ( variant == -4 ) ||
       ( variant == -5 ) ||
       //( variant == 0 ) ||
       //( variant == -1 ) ||
       //( variant == -2 ) ||
       //( variant == 1 ) ||
       //( variant == 2 ) ||
       //( variant == 3 ) ||
       //( variant == 4 ) ||
       FALSE
     )
  {
    *gflops      = 0.0;
    *dtime       = 0.0;
    *diff1       = 0.0;
    *diff2       = 0.0;
    *dtime_tred  = 0.0;
    *dtime_tevd  = 0.0;
    *dtime_appq  = 0.0;
    *gflops_tred = 0.0;
    *gflops_tevd = 0.0;
    *gflops_appq = 0.0;
    *k_perf      = 0;
    return;
  }

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Z );

  FLA_Copy_external( A, A_save );

  for ( i = 0 ; i < n_repeats; i++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case -3:
    {
      *k_perf = 0;
      REF_Hevd_lv( A, l,
                   dtime_tred, dtime_tevd, dtime_appq );
      break;
    }

    case -4:
    {
      *k_perf = 0;
      REF_Hevdd_lv( A, l,
                    dtime_tred, dtime_tevd, dtime_appq );
      break;
    }

    case -5:
    {
      *k_perf = 0;
      REF_Hevdr_lv( A, l, Z,
                    dtime_tred, dtime_tevd, dtime_appq );
      break;
    }

    case 0:
    {
      *k_perf = 0;
      REF_Hevd_lv_components( A, l,
                              dtime_tred, dtime_tevd, dtime_appq );
      break;
    }

    case -1:
    {
      *k_perf = 0;
      REF_Hevdd_lv_components( A, l,
                               dtime_tred, dtime_tevd, dtime_appq );
      break;
    }

    case -2:
    {
      *k_perf = 0;
      REF_Hevdr_lv_components( A, l, Z,
                               dtime_tred, dtime_tevd, dtime_appq );
      break;
    }

    // Time variant 1
    case 1:
    {
      *k_perf = FLA_Hevd_lv_var1_components( n_iter_max, A, l, k_accum, b_alg,
                                             dtime_tred, dtime_tevd, dtime_appq );
      break;
    }

    // Time variant 2
    case 2:
    {
      *k_perf = FLA_Hevd_lv_var2_components( n_iter_max, A, l, k_accum, b_alg,
                                             dtime_tred, dtime_tevd, dtime_appq );
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    if ( *dtime < dtime_save )
    {
      dtime_save      = *dtime;
      dtime_tred_save = *dtime_tred;
      dtime_tevd_save = *dtime_tevd;
      dtime_appq_save = *dtime_appq;
    }
  }

  *dtime      = dtime_save;
  *dtime_tred = dtime_tred_save;
  *dtime_tevd = dtime_tevd_save;
  *dtime_appq = dtime_appq_save;

//if ( variant == -3 || variant == 0 )
//printf( "\ndtime is %9.3e\n", *dtime );

  {
    FLA_Obj V, A_rev_evd, norm, eye;

    if ( variant == -2 || variant == -5 ) FLA_Copy( Z, A );

    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &V ); 
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_rev_evd ); 
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &eye ); 
    FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

    FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, l, A );

    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, A, V, FLA_ZERO, A_rev_evd );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A_rev_evd );

//FLA_Obj_show( "A_rev_evd", A_rev_evd, "%9.2e + %9.2e ", "" );
 
    FLA_Axpy( FLA_MINUS_ONE, A_save, A_rev_evd );
    FLA_Norm_frob( A_rev_evd, norm );
    FLA_Obj_extract_real_scalar( norm, diff1 );

    FLA_Set_to_identity( eye );
	FLA_Copy( V, A_rev_evd );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, V, A_rev_evd, FLA_MINUS_ONE, eye );
    FLA_Norm_frob( eye, norm );
    FLA_Obj_extract_real_scalar( norm, diff2 );

    FLA_Obj_free( &V );
    FLA_Obj_free( &A_rev_evd );
    FLA_Obj_free( &eye );
    FLA_Obj_free( &norm );
  }

  k = 2.00;

  flops_tred = ( ( 4.0 / 3.0 )   * m * m * m );
  flops_tevd = (   4.5           * k * m * m     +
                   3.0           * k * m * m * m );

  if ( variant == -1 || variant == -2 ||
       variant == -4 || variant == -5 )
    flops_appq = ( 2.0           * m * m * m );
  else
    flops_appq = ( 4.0 / 3.0     * m * m * m );

/*
  if ( FLA_Obj_is_complex( A ) )
  {
    *gflops      = ( 4.0 * flops_tred + 
                     2.0 * flops_tevd + 
                     4.0 * flops_appq ) / *dtime      / 1e9;

    *gflops_tred = ( 4.0 * flops_tred ) / *dtime_tred / 1e9;
    *gflops_tevd = ( 2.0 * flops_tevd ) / *dtime_tevd / 1e9;
    *gflops_appq = ( 4.0 * flops_appq ) / *dtime_appq / 1e9;
  }
  else
  {
    *gflops      = ( 1.0 * flops_tred + 
                     1.0 * flops_tevd + 
                     1.0 * flops_appq ) / *dtime      / 1e9;

    *gflops_tred = ( 1.0 * flops_tred ) / *dtime_tred / 1e9;
    *gflops_tevd = ( 1.0 * flops_tevd ) / *dtime_tevd / 1e9;
    *gflops_appq = ( 1.0 * flops_appq ) / *dtime_appq / 1e9;
  }
*/

  if ( FLA_Obj_is_complex( A ) )
  {
    mult_tred = 4.0;
    mult_tevd = 2.0;
    mult_appq = 4.0;
  }
  else
  {
    mult_tred = 1.0;
    mult_tevd = 1.0;
    mult_appq = 1.0;
  }

  *gflops = ( mult_tred * flops_tred + 
              mult_tevd * flops_tevd + 
              mult_appq * flops_appq ) / *dtime / 1e9;

  *gflops_tred = ( mult_tred * flops_tred ) / *dtime_tred / 1e9;
  *gflops_tevd = ( mult_tevd * flops_tevd ) / *dtime_tevd / 1e9;
  *gflops_appq = ( mult_appq * flops_appq ) / *dtime_appq / 1e9;

  FLA_Copy_external( A_save, A );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &Z );
}

