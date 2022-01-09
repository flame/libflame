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


void time_Transpose(
                  integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
                  FLA_Obj A, FLA_Obj A_ref,
                  double *dtime, double *diff, double *gflops );


void time_Transpose(
                  integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
                  FLA_Obj A, FLA_Obj A_ref,
                  double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_old, A_tmp;

  fla_blocksize_t*
    bp;
  fla_transpose_t*
    cntl_trans_var_unb;
  fla_transpose_t*
    cntl_trans_var_blk;
  fla_swap_t*
    cntl_swap_var_blk;
  fla_swap_t*
    cntl_swap_blas;


  bp                 = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  cntl_swap_blas     = FLA_Cntl_swap_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL );
  cntl_swap_var_blk  = FLA_Cntl_swap_obj_create( FLA_FLAT, FLA_UNBLOCKED_VARIANT1, bp, cntl_swap_blas );
  cntl_trans_var_unb = FLA_Cntl_transpose_obj_create( FLA_FLAT, FLA_UNBLOCKED_VARIANT1, NULL, NULL, NULL );
  cntl_trans_var_blk = FLA_Cntl_transpose_obj_create( FLA_FLAT, variant, bp, cntl_trans_var_unb, cntl_swap_var_blk );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_old );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_tmp );

  FLA_Copy_external( A, A_old );


  for ( irep = 0 ; irep < nrepeats; irep++ ){
    FLA_Copy_external( A_old, A );

    *dtime = FLA_Clock();

    switch( variant ){
    case 0:

      //FLA_Copyt_external( FLA_TRANSPOSE, A, A_tmp );
      //FLA_Set( FLA_ZERO, A );
      //FLA_Copyt_external( FLA_NO_TRANSPOSE, A_tmp, A );
      FLA_Transpose( A );

      break;

    case 1:{

      /* Time variant 1 */
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Transpose_unb_var1( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Transpose_blk_var1( A, cntl_trans_var_blk );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:{

      /* Time variant 2 */
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Transpose_unb_var2( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Transpose_blk_var2( A, cntl_trans_var_blk );
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

  FLA_Cntl_obj_free( cntl_trans_var_blk );
  FLA_Cntl_obj_free( cntl_trans_var_unb );
  FLA_Cntl_obj_free( cntl_swap_var_blk );
  FLA_Cntl_obj_free( cntl_swap_blas );
  FLA_Blocksize_free( bp );

  if ( variant == 0 ){
    FLA_Copy_external( A, A_ref );
    *diff = 0.0;
  }
  else{
    *diff = FLA_Max_elemwise_diff( A, A_ref );
  }

  *gflops = 4 * n * n /
            dtime_old / 1e9;

  *dtime = dtime_old;

  FLA_Copy_external( A_old, A );

  FLA_Obj_free( &A_old );
  FLA_Obj_free( &A_tmp );
}

