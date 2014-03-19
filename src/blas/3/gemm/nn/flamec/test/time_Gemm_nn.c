/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_UNBLOCKED 1
#define FLA_ALG_BLOCKED   2


void time_Gemm_nn(
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Gemm_nn( 
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    C_old;

  fla_blocksize_t*
    bp;
  fla_gemm_t*
    cntl_gemm_blas;
  fla_gemm_t*
    cntl_gemm_var;

  bp             = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  cntl_gemm_blas = FLA_Cntl_gemm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL );
  cntl_gemm_var  = FLA_Cntl_gemm_obj_create( FLA_FLAT, variant, bp, cntl_gemm_blas );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );

  FLA_Copy_external( C, C_old );


  for ( irep = 0 ; irep < nrepeats; irep++ ){
    FLA_Copy_external( C_old, C );

    *dtime = FLA_Clock();

    switch( variant ){
    // Time reference implementation
    case 0:
      REF_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                FLA_ONE, A, B, FLA_ONE, C );
  //FLA_Ger_external( FLA_ONE, A, B, C );
  //FLA_Obj_show( "C_ref = [", C, "%f ", "]" );
      break;

    // Time variant 1
    case 1:{
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Gemm_nn_unb_var1( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Gemm_nn_blk_var1( FLA_ONE, A, B, FLA_ONE, C, cntl_gemm_var );
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
        FLA_Gemm_nn_unb_var2( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Gemm_nn_blk_var2( FLA_ONE, A, B, FLA_ONE, C, cntl_gemm_var );
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
        FLA_Gemm_nn_unb_var3( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Gemm_nn_blk_var3( FLA_ONE, A, B, FLA_ONE, C, cntl_gemm_var );
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
        FLA_Gemm_nn_unb_var4( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Gemm_nn_blk_var4( FLA_ONE, A, B, FLA_ONE, C, cntl_gemm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time variant 5
    case 5:{
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Gemm_nn_unb_var5( FLA_ONE, A, B, FLA_ONE, C );
 //FLA_Ger_external( FLA_ONE, A, B, C );
 //FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Gemm_nn_blk_var5( FLA_ONE, A, B, FLA_ONE, C, cntl_gemm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time variant 6
    case 6:{
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Gemm_nn_unb_var6( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Gemm_nn_blk_var6( FLA_ONE, A, B, FLA_ONE, C, cntl_gemm_var );
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

  FLA_Cntl_obj_free( cntl_gemm_var );
  FLA_Cntl_obj_free( cntl_gemm_blas );
  FLA_Blocksize_free( bp );

//FLA_Obj_show( "C_after = [", C, "%f ", "]" );

  if ( variant == 0 )
  {
    FLA_Copy_external( C, C_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLA_Max_elemwise_diff( C, C_ref );
  }

  *gflops = 2.0 * 
            FLA_Obj_length( C ) * 
            FLA_Obj_width( C ) * 
            FLA_Obj_width( A ) / 
            dtime_old / 
            1.0e9;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

