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


void time_Chol_l(
                  integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
                  FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj norm,
                  double *dtime, double *diff, double *gflops );


void time_Chol_l(
                  integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
                  FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj norm,
                  double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_save = 1.0e9;

  FLA_Obj
    A_save, b_save, b_orig_save;

  fla_blocksize_t* bp;
  fla_chol_t* cntl_chol_var;
  fla_chol_t* cntl_chol_unb;
  fla_herk_t* cntl_herk_blas;
  fla_trsm_t* cntl_trsm_blas;
  fla_gemm_t* cntl_gemm_blas;

/*
  if( type == FLA_ALG_UNBLOCKED && n > 400 )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }
*/

  bp               = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  cntl_chol_unb    = FLA_Cntl_chol_obj_create( FLA_FLAT, FLA_UNB_OPT_VARIANT2, NULL, NULL, NULL, NULL, NULL );
  cntl_herk_blas   = FLA_Cntl_herk_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL, NULL, NULL );
  cntl_trsm_blas   = FLA_Cntl_trsm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL, NULL, NULL );
  cntl_gemm_blas   = FLA_Cntl_gemm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL, NULL );
  cntl_chol_var    = FLA_Cntl_chol_obj_create( FLA_FLAT,
                                               variant,
                                               bp,
                                               cntl_chol_unb,
                                               cntl_herk_blas,
                                               cntl_trsm_blas,
                                               cntl_gemm_blas );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b, &b_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b_orig, &b_orig_save );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( b, b_save );
  FLA_Copy_external( b_orig, b_orig_save );


  for ( irep = 0 ; irep < nrepeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:

      REF_Chol_l( A );

      break;

    case 1:{

      // Time variant 1
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Chol_l_unb_var1( A );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Chol_l_opt_var1( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Chol_l_blk_var1( A, cntl_chol_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:{

      // Time variant 2
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Chol_l_unb_var2( A );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Chol_l_opt_var2( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Chol_l_blk_var2( A, cntl_chol_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    } 
    case 3:{

      // Time variant 3 
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Chol_l_unb_var3( A );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Chol_l_opt_var3( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Chol_l_blk_var3( A, cntl_chol_var );
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

  FLA_Cntl_obj_free( cntl_chol_var );
  FLA_Cntl_obj_free( cntl_chol_unb );
  FLA_Cntl_obj_free( cntl_herk_blas );
  FLA_Cntl_obj_free( cntl_trsm_blas );
  FLA_Cntl_obj_free( cntl_gemm_blas );
  FLA_Blocksize_free( bp );

  {
    FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, A, b );
    FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                       FLA_NONUNIT_DIAG, A, b );

    FLA_Hemv_external( FLA_LOWER_TRIANGULAR,
                       FLA_ONE, A_save, b, FLA_MINUS_ONE, b_orig );

    FLA_Nrm2_external( b_orig, norm );
    FLA_Obj_extract_real_scalar( norm, diff );
  }

  *gflops = 1.0 / 3.0 * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) / 
            dtime_save / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_save;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( b_save, b );
  FLA_Copy_external( b_orig_save, b_orig );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &b_save );
  FLA_Obj_free( &b_orig_save );
}

