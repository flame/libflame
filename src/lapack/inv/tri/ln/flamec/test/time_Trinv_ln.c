
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_BLOCKED   1
#define FLA_ALG_UNBLOCKED 2
#define FLA_ALG_UNB_OPT   3


FLA_Error REF_Trinv_ln( FLA_Obj A );
void time_Trinv_ln(
                  int variant, int type, int nrepeats, int m, int nb_alg,
                  FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj norm,
                  double *dtime, double *diff, double *gflops );


void time_Trinv_ln(
                  int variant, int type, int nrepeats, int m, int nb_alg,
                  FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj norm,
                  double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, b_save, b_orig_save;

  fla_blocksize_t*
    bp;
  fla_trinv_t*
    cntl_trinv_var;
  fla_trinv_t*
    cntl_trinv_unb;
  fla_gemm_t*
    cntl_gemm_blas;
  fla_trmm_t*
    cntl_trmm_blas;
  fla_trsm_t*
    cntl_trsm_blas;


  bp                = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  cntl_trinv_unb    = FLA_Cntl_trinv_obj_create( FLA_FLAT, FLA_UNB_OPT_VARIANT3, NULL, NULL, NULL, NULL, NULL, NULL );
  cntl_trmm_blas    = FLA_Cntl_trmm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL, NULL );
  cntl_trsm_blas    = FLA_Cntl_trsm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL, NULL );
  cntl_gemm_blas    = FLA_Cntl_gemm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL );
  cntl_trinv_var    = FLA_Cntl_trinv_obj_create( FLA_FLAT, variant, bp, cntl_trinv_unb, cntl_trmm_blas, cntl_trsm_blas, cntl_trsm_blas, cntl_gemm_blas );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b, &b_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b_orig, &b_orig_save );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( b, b_save );
  FLA_Copy_external( b_orig, b_orig_save );


  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    // Time reference
    case 0:
      REF_Trinv_ln( A );
      break;

    // Time variant 1
    case 1:{
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Trinv_ln_unb_var1( A );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Trinv_ln_opt_var1( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Trinv_ln_blk_var1( A, cntl_trinv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time variant 2
    case 2:{
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Trinv_ln_unb_var2( A );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Trinv_ln_opt_var2( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Trinv_ln_blk_var2( A, cntl_trinv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    } 

    // Time variant 3 
    case 3:{
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Trinv_ln_unb_var3( A );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Trinv_ln_opt_var3( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Trinv_ln_blk_var3( A, cntl_trinv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time variant 4 
    case 4:{
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Trinv_ln_unb_var4( A );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Trinv_ln_opt_var4( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Trinv_ln_blk_var4( A, cntl_trinv_var );
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

  FLA_Cntl_obj_free( cntl_trinv_var );
  FLA_Cntl_obj_free( cntl_trinv_unb );
  FLA_Cntl_obj_free( cntl_gemm_blas );
  FLA_Cntl_obj_free( cntl_trmm_blas );
  FLA_Cntl_obj_free( cntl_trsm_blas );
  FLA_Blocksize_free( bp );

  {
    FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, A, b );

    FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, A_save, b );

    FLA_Axpy_external( FLA_MINUS_ONE, b_orig, b );

    FLA_Nrm2_external( b, norm );
    FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, norm,
                               1, 1, diff, 1, 1 );
  }

  *gflops = 1.0 / 3.0 * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) / 
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( b_save, b );
  FLA_Copy_external( b_orig_save, b_orig );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &b_save );
  FLA_Obj_free( &b_orig_save );
}

