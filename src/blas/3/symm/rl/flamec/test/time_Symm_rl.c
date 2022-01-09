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


void time_Symm_rl(
               integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Symm_rl( 
               integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9; 

  FLA_Obj
    C_old;

  fla_blocksize_t*
    bp;
  fla_gemm_t*
    cntl_gemm_blas;
  fla_symm_t*
    cntl_symm_blas;
  fla_symm_t*
    cntl_symm_var;

  if ( type == FLA_ALG_UNBLOCKED && n > 300 )
  {
    *diff = 0.0;
    *gflops = 0.0;
    return;
  }

  bp             = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  cntl_gemm_blas = FLA_Cntl_gemm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL );
  cntl_symm_blas = FLA_Cntl_symm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL, NULL, NULL );
  cntl_symm_var  = FLA_Cntl_symm_obj_create( FLA_FLAT, variant, bp, cntl_symm_blas, cntl_gemm_blas, cntl_gemm_blas );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );

  FLA_Copy_external( C, C_old );


  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLA_Copy_external( C_old, C );

    *dtime = FLA_Clock();

    switch( variant ){
    case 0:
      // Time reference implementation
      REF_Symm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_ONE, A, B, FLA_ONE, C );
      break;

    case 1:{
      // Time variant 1
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_rl_unb_var1( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_rl_blk_var1( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
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
        FLA_Symm_rl_unb_var2( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_rl_blk_var2( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
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
        FLA_Symm_rl_unb_var3( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_rl_blk_var3( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
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
        FLA_Symm_rl_unb_var4( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_rl_blk_var4( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
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
        FLA_Symm_rl_unb_var5( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_rl_blk_var5( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 6:{
      // Time variant 6
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_rl_unb_var6( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_rl_blk_var6( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 7:{
      // Time variant 7
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_rl_unb_var7( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_rl_blk_var7( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 8:{
      // Time variant 8
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_rl_unb_var8( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_rl_blk_var8( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 9:{
      // Time variant 9
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_rl_unb_var9( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_rl_blk_var9( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 10:{
      // Time variant 10
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_rl_unb_var10( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_rl_blk_var10( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
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

  FLA_Cntl_obj_free( cntl_symm_var );
  FLA_Cntl_obj_free( cntl_symm_blas );
  FLA_Cntl_obj_free( cntl_gemm_blas );
  FLA_Blocksize_free( bp );

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
            1e9;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

