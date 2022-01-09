/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_BLOCKED   1
#define FLA_ALG_UNBLOCKED 2
#define FLA_ALG_UNB_OPT   3


FLA_Error REF_Eig_gest_nu( FLA_Obj A, FLA_Obj B );
void time_Eig_gest_nu(
               integer variant, integer type, integer n_repeats, integer n, integer b_alg,
               FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj Y, FLA_Obj B,
               double *dtime, double *diff, double *gflops );

extern TLS_CLASS_SPEC fla_axpy_t*  fla_axpy_cntl_blas;
extern TLS_CLASS_SPEC fla_gemm_t*  fla_gemm_cntl_blas;
extern TLS_CLASS_SPEC fla_hemm_t*  fla_hemm_cntl_blas;
extern TLS_CLASS_SPEC fla_her2k_t* fla_her2k_cntl_blas;
extern TLS_CLASS_SPEC fla_trmm_t*  fla_trmm_cntl_blas;
extern TLS_CLASS_SPEC fla_trsm_t*  fla_trsm_cntl_blas;

void time_Eig_gest_nu(
               integer variant, integer type, integer n_repeats, integer n, integer b_alg,
               FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj Y, FLA_Obj B,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_save = 1.0e9;

  FLA_Obj
    A_save, B_save, norm;

  fla_blocksize_t* bp;
  fla_eig_gest_t*  cntl_eig_gest_var;
  fla_eig_gest_t*  cntl_eig_gest_unb;


  if ( ( type == FLA_ALG_UNBLOCKED || type == FLA_ALG_UNB_OPT ) &&
       n > 300 )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }

  if ( variant == 3 )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }


  bp                = FLA_Blocksize_create( b_alg, b_alg, b_alg, b_alg );
  cntl_eig_gest_unb = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
                                                    //FLA_UNBLOCKED_VARIANT1,
                                                    FLA_UNB_OPT_VARIANT1,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL );
  cntl_eig_gest_var = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
                                                    variant,
                                                    bp,
                                                    cntl_eig_gest_unb,
                                                    fla_axpy_cntl_blas,
                                                    fla_axpy_cntl_blas,
                                                    fla_gemm_cntl_blas,
                                                    fla_gemm_cntl_blas,
                                                    fla_gemm_cntl_blas,
                                                    fla_hemm_cntl_blas,
                                                    fla_her2k_cntl_blas,
                                                    fla_trmm_cntl_blas,
                                                    fla_trmm_cntl_blas,
                                                    fla_trsm_cntl_blas,
                                                    fla_trsm_cntl_blas );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B, &B_save );

  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( B, B_save );


  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );
    FLA_Copy_external( B_save, B );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:
      REF_Eig_gest_nu( A, B );
      break;

    case 1:
    {
      // Time variant 1
      switch( type )
      {
      case FLA_ALG_UNBLOCKED:
        FLA_Eig_gest_nu_unb_var1( A, Y, B );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Eig_gest_nu_opt_var1( A, Y, B );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Eig_gest_nu_blk_var1( A, Y, B, cntl_eig_gest_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:
    {
      // Time variant 2
      switch( type )
      {
      case FLA_ALG_UNBLOCKED:
        FLA_Eig_gest_nu_unb_var2( A, Y, B );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Eig_gest_nu_opt_var2( A, Y, B );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Eig_gest_nu_blk_var2( A, Y, B, cntl_eig_gest_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 3:
    {
      // Time variant 3
      switch( type )
      {
      case FLA_ALG_UNBLOCKED:
        //FLA_Eig_gest_nu_unb_var3( A, Y, B );
        break;
      case FLA_ALG_UNB_OPT:
        //FLA_Eig_gest_nu_opt_var3( A, Y, B );
        break;
      case FLA_ALG_BLOCKED:
        //FLA_Eig_gest_nu_blk_var3( A, Y, B, cntl_eig_gest_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 4:
    {
      // Time variant 4
      switch( type )
      {
      case FLA_ALG_UNBLOCKED:
        FLA_Eig_gest_nu_unb_var4( A, Y, B );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Eig_gest_nu_opt_var4( A, Y, B );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Eig_gest_nu_blk_var4( A, Y, B, cntl_eig_gest_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 5:
    {
      // Time variant 5
      switch( type )
      {
      case FLA_ALG_UNBLOCKED:
        FLA_Eig_gest_nu_unb_var5( A, Y, B );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Eig_gest_nu_opt_var5( A, Y, B );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Eig_gest_nu_blk_var5( A, Y, B, cntl_eig_gest_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_save = min( *dtime, dtime_save );
  }

  FLA_Cntl_obj_free( cntl_eig_gest_var );
  FLA_Cntl_obj_free( cntl_eig_gest_unb );
  FLA_Blocksize_free( bp );

  // Recover A.
  if ( inv == FLA_NO_INVERSE )
  {
    if ( uplo == FLA_LOWER_TRIANGULAR )
    {
      // A = L' * A_orig * L
      // A_orig = inv(L') * A * inv(L)
      FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A );
      FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
      FLA_Trsm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
    }
    else // if ( uplo == FLA_UPPER_TRIANGULAR )
    {
      // A = U * A_orig * U'
      // A_orig = inv(U) * A * inv(U')
      FLA_Hermitianize( FLA_UPPER_TRIANGULAR, A );
      FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
      FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
    }
  }
  else // if ( inv == FLA_INVERSE )
  {
    if ( uplo == FLA_LOWER_TRIANGULAR )
    {
      // A = inv(L) * A_orig * inv(L')
      // A_orig = L * A * L'
      FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A );
      FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
      FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
    }
    else // if ( uplo == FLA_UPPER_TRIANGULAR )
    {
      // A = inv(U') * A_orig * inv(U)
      // A_orig = U' * A * U
      FLA_Hermitianize( FLA_UPPER_TRIANGULAR, A );
      FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
      FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
    }
  }

  *diff = FLA_Max_elemwise_diff( A, A_save );

/*
if ( type == FLA_ALG_UNBLOCKED )
{
FLA_Obj_show( "A", A, "%10.3e", "" );
FLA_Obj_show( "A_orig", A_save, "%10.3e", "" );
}
*/

  *gflops = 1.0 * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) / 
            dtime_save / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_save;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( B_save, B );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &B_save );
  FLA_Obj_free( &norm );
}

