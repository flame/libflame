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


void time_LU_piv(
                  integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
                  FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj p, FLA_Obj norm, 
                  double *dtime, double *diff, double *gflops );


void time_LU_piv(
                  integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
                  FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj p, FLA_Obj norm,
                  double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_save = 1.0e9;

  FLA_Obj
    A_save, b_save, b_orig_save;

  fla_blocksize_t*
    bp;
  fla_lu_t*
    cntl_lu_var;
  fla_lu_t*
    cntl_lu_lapack;
  fla_trsm_t*
    cntl_trsm_blas;
  fla_gemm_t*
    cntl_gemm_blas;
  fla_appiv_t*
    cntl_appiv_unb;


  bp             = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  cntl_appiv_unb = FLA_Cntl_appiv_obj_create( FLA_FLAT, FLA_UNB_OPT_VARIANT1, NULL, NULL );
  cntl_lu_lapack = FLA_Cntl_lu_obj_create( FLA_FLAT, FLA_UNB_OPT_VARIANT5, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL );
  cntl_trsm_blas = FLA_Cntl_trsm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL, NULL );
  cntl_gemm_blas = FLA_Cntl_gemm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL );
  cntl_lu_var    = FLA_Cntl_lu_obj_create( FLA_FLAT, variant, bp, cntl_lu_lapack, cntl_gemm_blas, cntl_gemm_blas, cntl_gemm_blas, cntl_trsm_blas, cntl_trsm_blas, cntl_appiv_unb, cntl_appiv_unb );


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
      REF_LU_piv( A, p );

      break;

    case 3:{

      // Time variant 3 
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_LU_piv_unb_var3( A, p );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_LU_piv_opt_var3( A, p );
        break;
      case FLA_ALG_BLOCKED:
        FLA_LU_piv_blk_var3( A, p, cntl_lu_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 4:{

      // Time variant 4 
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_LU_piv_unb_var4( A, p );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_LU_piv_opt_var4( A, p );
        break;
      case FLA_ALG_BLOCKED:
        FLA_LU_piv_blk_var4( A, p, cntl_lu_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 5:{

      // Time variant 5 
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_LU_piv_unb_var5( A, p );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_LU_piv_opt_var5( A, p );
        break;
      case FLA_ALG_BLOCKED:
        FLA_LU_piv_blk_var5( A, p, cntl_lu_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    }

    *dtime = FLA_Clock() - *dtime;
    dtime_save = fla_min( *dtime, dtime_save );
  }

  FLA_Cntl_obj_free( cntl_lu_var );
  FLA_Cntl_obj_free( cntl_lu_lapack );
  FLA_Cntl_obj_free( cntl_appiv_unb );
  FLA_Cntl_obj_free( cntl_trsm_blas );
  FLA_Cntl_obj_free( cntl_gemm_blas );
  FLA_Blocksize_free( bp );

  if ( type == FLA_ALG_REFERENCE )
  {
    FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, b );

    FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_UNIT_DIAG, A, b );
    FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, A, b );

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE,
                       A_save, b, FLA_MINUS_ONE, b_orig );

    FLA_Nrm2_external( b_orig, norm );
    FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, norm,
                               1, 1, diff, 1, 1 );
  }
  else
  {
/*
    FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, b );

    FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_UNIT_DIAG, A, b );
    FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, A, b );
*/
    FLA_LU_piv_solve( A, p, b_orig, b );
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE,
                       A_save, b, FLA_MINUS_ONE, b_orig );

    FLA_Nrm2_external( b_orig, norm );
    FLA_Obj_extract_real_scalar( norm, diff );
/*
    FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, norm,
                               1, 1, diff, 1, 1 );
*/
  }

  *gflops = 2.0 / 3.0 * 
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

