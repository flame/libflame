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

FLA_Error REF_Svd_uv( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V );
FLA_Error REF_Svdd_uv( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V );

void time_Svd_uv(
               int variant, int type, int n_repeats, int m, int n, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
               double *dtime, double *diff1, double* diff2, double *gflops, int* k_perf );


void time_Svd_uv(
               int variant, int type, int n_repeats, int m, int n, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
               double *dtime, double *diff1, double* diff2, double *gflops, int* k_perf )
{
  int irep;

  double
    k, dtime_old = 1.0e9;

  FLA_Obj
    A_save;

  if (
       //( variant == 0 ) ||
       //( variant == -1 ) ||
       //( variant == 1 ) ||
       //( variant == 2 ) ||
       FALSE
     )
  {
    *gflops = 0.0;
    *dtime  = 0.0;
    *diff1  = 0.0;
    *diff2  = 0.0;
    *k_perf = 0;
    return;
  }

/*
  if ( variant == 0 && max( m, n ) > 800 )
  {
    *gflops = 0.0;
    *dtime  = 0.0;
    *diff1  = 0.0;
    *diff2  = 0.0;
    *k_perf = 0;
    return;
  }
*/

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );

  FLA_Copy_external( A, A_save );


  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:
      REF_Svd_uv( A, s, U, V );
      *k_perf = 0;
      break;

    case -1:
      REF_Svdd_uv( A, s, U, V );
      *k_perf = 0;
      break;

    // Time variant 1
    case 1:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        *k_perf = FLA_Svd_uv_unb_var1( n_iter_max, A, s, U, V, k_accum, b_alg );
        break;
      }
      break;
    }

    // Time variant 2
    case 2:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        *k_perf = FLA_Svd_uv_unb_var2( n_iter_max, A, s, U, V, k_accum, b_alg );
        break;
      }
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );

  }
  {
    FLA_Obj USVh, eyeU, eyeV, normU, normV;
    double  diffU, diffV;

//FLA_Obj_show( "A_save", A_save, "%9.2e + %9.2e ", "" );
//if ( variant == 1 )
//FLA_Obj_show( "U", U, "%8.1e %8.1e ", "" );
//FLA_Obj_show( "V'", V, "%8.1e %8.1e ", "" );

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
    //*diff2 = diffU;
    //*diff2 = diffV;

    if ( variant < 1 ) // LAPACK computes and stores V'
    {
      FLA_Obj UL, UR;
      FLA_Obj VT,
              VB;

      FLA_Part_1x2( U,   &UL, &UR,   min( m, n ), FLA_LEFT );
      FLA_Part_2x1( V,   &VT,
                         &VB,        min( m, n ), FLA_TOP );
/*
FLA_Obj_show( "A_orig", A_save, "%8.1e + %8.1e ", "" );
FLA_Obj_show( "s", s, "%8.1e ", "" );
FLA_Obj_show( "U", U, "%8.1e + %8.1e ", "" );
FLA_Obj_show( "V'", V, "%8.1e + %8.1e ", "" );
*/

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

  if ( FLA_Obj_is_complex( A ) )
  {
    *gflops = ( 4.0 * ( 4.0           * m * n * n - 4.0 / 3.0 * n * n * n ) + 
                4.0 * ( 4.0           * m * m * n - 4.0       * m * n * n
                                                  + 4.0 / 3.0 * n * n * n ) + 
                4.0 * ( 4.0 / 3.0     * n * n * n ) + 
                      ( 13.0      * k * n * n     ) +
                2.0 * (       3.0 * k * m * n * n ) + 
                2.0 * (       3.0 * k * n * n * n ) ) / 
              dtime_old / 1e9;
  }
  else
  {
    *gflops = ( 1.0 * ( 4.0           * m * n * n - 4.0 / 3.0 * n * n * n ) + 
                1.0 * ( 4.0           * m * m * n - 4.0       * m * n * n
                                                  + 4.0 / 3.0 * n * n * n ) + 
                1.0 * ( 4.0 / 3.0     * n * n * n ) + 
                      ( 13.0      * k * n * n     ) +
                1.0 * (       3.0 * k * m * n * n ) + 
                1.0 * (       3.0 * k * n * n * n ) ) / 
              dtime_old / 1e9;
  }


  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );

  FLA_Obj_free( &A_save );
}

