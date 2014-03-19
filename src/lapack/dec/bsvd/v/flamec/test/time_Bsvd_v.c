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

FLA_Error REF_Bsvd_v( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, FLA_Obj U, FLA_Obj V );

void time_Bsvd_v(
               int variant, int type, int n_repeats, int m, int n, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A_orig, FLA_Obj d, FLA_Obj e, FLA_Obj s, FLA_Obj G, FLA_Obj H, FLA_Obj RG, FLA_Obj RH, FLA_Obj W, FLA_Obj U, FLA_Obj V,
               double *dtime, double *diff1, double* diff2, double *gflops );


void time_Bsvd_v(
               int variant, int type, int n_repeats, int m, int n, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A_orig, FLA_Obj d, FLA_Obj e, FLA_Obj s, FLA_Obj G, FLA_Obj H, FLA_Obj RG, FLA_Obj RH, FLA_Obj W, FLA_Obj U, FLA_Obj V,
               double *dtime, double *diff1, double* diff2, double *gflops )
{
  int irep;

  double
    k, dtime_old = 1.0e9;

  FLA_Obj
    d_save, e_save, U_save, V_save;

  if (
       //( variant == 0 ) ||
       //( variant == 1 ) ||
       //( variant == 2 ) ||
       FALSE
     )
  {
    *gflops = 0.0;
    *dtime  = 0.0;
    *diff1  = 0.0;
    *diff2  = 0.0;
    return;
  }

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, d, &d_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, e, &e_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, U, &U_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, V, &V_save );

  FLA_Copy_external( d, d_save );
  FLA_Copy_external( e, e_save );
  FLA_Copy_external( U, U_save );
  FLA_Copy_external( V, V_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( d_save, d );
    FLA_Copy_external( e_save, e );
    FLA_Copy_external( U_save, U );
    FLA_Copy_external( V_save, V );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:
      REF_Bsvd_v( FLA_UPPER_TRIANGULAR, d, e, U, V );
      break;

    // Time variant 1
    case 1:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Bsvd_v_opt_var1( n_iter_max, d, e, G, H, U, V, b_alg );
        break;
      }
      break;
    }

    // Time variant 2
    case 2:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Bsvd_v_opt_var2( n_iter_max, d, e, G, H, RG, RH, W, U, V, b_alg );
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

	FLA_Copy( d, s );

//FLA_Obj_show( "s before", s, "%9.2e ", "" );

    //FLA_Sort_svd( FLA_FORWARD, s, U, V );

//FLA_Obj_show( "s after", s, "%9.2e ", "" );

    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A_orig, &USVh ); 
    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, U,      &eyeU ); 
    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, V,      &eyeV ); 
    FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A_orig ), 1, 1, 0, 0, &normU );
    FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A_orig ), 1, 1, 0, 0, &normV );

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

//FLA_Obj_show( "V ref", V, "%8.1e + %8.1e ", "" );
/*
{
FLA_Obj Vh;
FLA_Obj_create_copy_of( FLA_CONJ_TRANSPOSE, V, &Vh );
FLA_Obj_show( "Vref'", Vh, "%8.1e %8.1e ", "" );
FLA_Obj_free( &Vh );
}
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
/*
FLA_Obj_show( "A_orig", A_orig, "%8.1e + %8.1e ", "" );
FLA_Obj_show( "s", s, "%8.1e ", "" );
FLA_Obj_show( "U", U, "%8.1e + %8.1e ", "" );
{
FLA_Obj Vh;
FLA_Obj_create_copy_of( FLA_CONJ_TRANSPOSE, V, &Vh );
FLA_Obj_show( "V'", Vh, "%8.1e + %8.1e ", "" );
FLA_Obj_free( &Vh );
}
*/

      FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, s, UL );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                FLA_ONE, UL, VL, FLA_ZERO, USVh );
    }

    FLA_Axpy( FLA_MINUS_ONE, A_orig, USVh );
    FLA_Norm_frob( USVh, normU );
    FLA_Obj_extract_real_scalar( normU, diff1 );

    FLA_Obj_free( &USVh );
    FLA_Obj_free( &eyeU );
    FLA_Obj_free( &eyeV );
    FLA_Obj_free( &normU );
    FLA_Obj_free( &normV );
  }

  k = 2.00;

  if ( FLA_Obj_is_complex( A_orig ) )
  { 
    *gflops = (       (      13.0 * k * n * n     ) +
                2.0 * (       3.0 * k * m * n * n ) + 
                2.0 * (       3.0 * k * n * n * n ) ) /
              dtime_old / 1e9;
  }
  else
  { 
    *gflops = (       (      13.0 * k * n * n     ) +
                1.0 * (       3.0 * k * m * n * n ) + 
                1.0 * (       3.0 * k * n * n * n ) ) /
              dtime_old / 1e9;
  }

  *dtime = dtime_old;

  FLA_Copy_external( d_save, d );
  FLA_Copy_external( e_save, e );
  FLA_Copy_external( U_save, U );
  FLA_Copy_external( V_save, V );

  FLA_Obj_free( &d_save );
  FLA_Obj_free( &e_save );
  FLA_Obj_free( &U_save );
  FLA_Obj_free( &V_save );
}

