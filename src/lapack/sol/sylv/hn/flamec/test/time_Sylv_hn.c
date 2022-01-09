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


FLA_Error REF_Sylv_hn( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );

void time_Sylv_hn(
                   integer variant, integer type, integer n_repeats, integer m, integer n, integer nb_alg,
                   FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref, FLA_Obj scale,
                   double *dtime, double *diff, double *gflops );


void time_Sylv_hn(
                   integer variant, integer type, integer n_repeats, integer m, integer n, integer nb_alg,
                   FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref, FLA_Obj scale,
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
  fla_sylv_t*
    cntl_sylv_var;
  fla_sylv_t*
    cntl_sylv_unb;
  fla_gemm_t*
    cntl_gemm_blas;

/*
  if( type == FLA_ALG_UNBLOCKED && n > 400 )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }
*/

  bp               = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  cntl_sylv_unb    = FLA_Cntl_sylv_obj_create( FLA_FLAT, FLA_UNB_OPT_VARIANT1, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL );
  cntl_gemm_blas   = FLA_Cntl_gemm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL, NULL );
  cntl_sylv_var    = FLA_Cntl_sylv_obj_create( FLA_FLAT, variant, bp, cntl_sylv_unb, cntl_sylv_unb, cntl_sylv_unb, cntl_gemm_blas, cntl_gemm_blas, cntl_gemm_blas, cntl_gemm_blas, cntl_gemm_blas, cntl_gemm_blas, cntl_gemm_blas, cntl_gemm_blas );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );

  FLA_Copy_external( C, C_old );


  for ( irep = 0 ; irep < n_repeats; irep++ ){
    FLA_Copy_external( C_old, C );

    *dtime = FLA_Clock();

    switch( variant ){
    case 0:
      /* Time reference implementation */
      REF_Sylv_hn( isgn, A, B, C, scale );

      break;

    case 1:{

      /* Time variant 1 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var1( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var1( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:{

      /* Time variant 2 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var2( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var2( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 3:{

      /* Time variant 3 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var3( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var3( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 4:{

      /* Time variant 4 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var4( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var4( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 5:{

      /* Time variant 5 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var5( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var5( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 6:{

      /* Time variant 6 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var6( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var6( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 7:{

      /* Time variant 7 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var7( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var7( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 8:{

      /* Time variant 8 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var8( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var8( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 9:{

      /* Time variant 9 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var9( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var9( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 10:{

      /* Time variant 10 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var10( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var10( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 11:{

      /* Time variant 11 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var11( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var11( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 12:{

      /* Time variant 12 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var12( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var12( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 13:{

      /* Time variant 13 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var13( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var13( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 14:{

      /* Time variant 14 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var14( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var14( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 15:{

      /* Time variant 15 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var15( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var15( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 16:{

      /* Time variant 16 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var16( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var16( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 17:{

      /* Time variant 17 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var17( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var17( isgn, A, B, C, scale, cntl_sylv_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 18:{

      /* Time variant 18 */
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Sylv_hn_opt_var18( isgn, A, B, C, scale );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Sylv_hn_blk_var18( isgn, A, B, C, scale, cntl_sylv_var );
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

  FLA_Cntl_obj_free( cntl_sylv_var );
  FLA_Cntl_obj_free( cntl_sylv_unb );
  FLA_Cntl_obj_free( cntl_gemm_blas );
  FLA_Blocksize_free( bp );

  if ( variant == 0 ){
    FLA_Copy_external( C, C_ref );
    *diff = 0.0;
  }
  else{
    *diff = FLA_Max_elemwise_diff( C, C_ref );
  }

  *gflops = ( m * m * n + n * n * m ) / 
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( C ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

