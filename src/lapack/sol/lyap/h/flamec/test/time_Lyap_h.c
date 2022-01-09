/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "FLA_Lyap_h.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_BLOCKED   1
#define FLA_ALG_UNBLOCKED 2
#define FLA_ALG_UNB_OPT   3

FLA_Error REF_Lyap_h( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale );

void time_Lyap_h(
               integer variant, integer type, integer n_repeats, integer m, integer nb_alg,
               FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj C_ref, FLA_Obj scale,
               double *dtime, double *diff, double *gflops );

extern TLS_CLASS_SPEC fla_scal_t*  fla_scal_cntl_blas;
extern TLS_CLASS_SPEC fla_gemm_t*  fla_gemm_cntl_blas;
extern TLS_CLASS_SPEC fla_hemm_t*  fla_hemm_cntl_blas;
extern TLS_CLASS_SPEC fla_her2k_t* fla_her2k_cntl_blas;
extern TLS_CLASS_SPEC fla_sylv_t*  fla_sylv_cntl;
extern TLS_CLASS_SPEC fla_lyap_t*  fla_lyap_cntl_leaf;

void time_Lyap_h(
               integer variant, integer type, integer n_repeats, integer m, integer nb_alg,
               FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj C_ref, FLA_Obj scale,
               double *dtime, double *diff, double *gflops )
{
  integer irep;

  double dtime_old = 1.0e9;

  FLA_Obj C_save, norm;

  fla_blocksize_t* bp;
  fla_lyap_t*      cntl_lyap_unb;
  fla_lyap_t*      cntl_lyap_opt;
  fla_lyap_t*      cntl_lyap_blk;

  if ( type == FLA_ALG_UNB_OPT && variant > 4 )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }

  bp               = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  cntl_lyap_unb    = FLA_Cntl_lyap_obj_create( FLA_FLAT,
                                               FLA_UNB_VAR_OFFSET + variant,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL );
  cntl_lyap_opt    = FLA_Cntl_lyap_obj_create( FLA_FLAT,
                                               FLA_OPT_VAR_OFFSET + variant,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL );
  cntl_lyap_blk    = FLA_Cntl_lyap_obj_create( FLA_FLAT,
                                               FLA_BLK_VAR_OFFSET + variant,
                                               bp,
                                               fla_scal_cntl_blas,
                                               fla_lyap_cntl_leaf,
                                               fla_sylv_cntl,
                                               fla_gemm_cntl_blas,
                                               fla_gemm_cntl_blas,
                                               fla_hemm_cntl_blas,
                                               fla_her2k_cntl_blas );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_save );
  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( C ), 1, 1, 0, 0, &norm );

  FLA_Copy_external( C, C_save );


  for ( irep = 0 ; irep < n_repeats; irep++ )
  {
    FLA_Copy_external( C_save, C );

    *dtime = FLA_Clock();

    switch( variant )
    {

    case 0:
      REF_Lyap_h( isgn, A, C, scale );
      break;

    case 1:
    {
      switch( type )
      {
        case FLA_ALG_UNBLOCKED:
          FLA_Lyap_h_unb_var1( isgn, A, C );
          break;
        case FLA_ALG_UNB_OPT:
          FLA_Lyap_h_opt_var1( isgn, A, C );
          break;
        case FLA_ALG_BLOCKED:
          FLA_Lyap_h_blk_var1( isgn, A, C, scale, cntl_lyap_blk );
          break;
      }
      break;
    }

    case 2:
    {
      switch( type )
      {
        case FLA_ALG_UNBLOCKED:
          FLA_Lyap_h_unb_var2( isgn, A, C );
          break;
        case FLA_ALG_UNB_OPT:
          FLA_Lyap_h_opt_var2( isgn, A, C );
          break;
        case FLA_ALG_BLOCKED:
          FLA_Lyap_h_blk_var2( isgn, A, C, scale, cntl_lyap_blk );
          break;
      }
      break;
    }

    case 3:
    {
      switch( type )
      {
        case FLA_ALG_UNBLOCKED:
          FLA_Lyap_h_unb_var3( isgn, A, C );
          break;
        case FLA_ALG_UNB_OPT:
          FLA_Lyap_h_opt_var3( isgn, A, C );
          break;
        case FLA_ALG_BLOCKED:
          FLA_Lyap_h_blk_var3( isgn, A, C, scale, cntl_lyap_blk );
          break;
      }
      break;
    }

    case 4:
    {
      switch( type )
      {
        case FLA_ALG_UNBLOCKED:
          FLA_Lyap_h_unb_var4( isgn, A, C );
          break;
        case FLA_ALG_UNB_OPT:
          FLA_Lyap_h_opt_var4( isgn, A, C );
          break;
        case FLA_ALG_BLOCKED:
          FLA_Lyap_h_blk_var4( isgn, A, C, scale, cntl_lyap_blk );
          break;
      }
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );
  }


  FLA_Blocksize_free( bp );
  FLA_Cntl_obj_free( cntl_lyap_unb );
  FLA_Cntl_obj_free( cntl_lyap_opt );
  FLA_Cntl_obj_free( cntl_lyap_blk );


/*
  if ( variant == 0 )
  {
    FLA_Copy_external( C, C_ref );
    *diff = 0.0;
  }
  else
  {
    FLA_Hermitianize( FLA_UPPER_TRIANGULAR, C );
    *diff = FLA_Max_elemwise_diff( C, C_ref );
  }
*/
  {
    FLA_Obj X, W;

    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &X );
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &W );

    FLA_Copy( C, X );
    FLA_Hermitianize( FLA_UPPER_TRIANGULAR, X );

    FLA_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, X, FLA_ZERO, W );
    FLA_Gemm( FLA_NO_TRANSPOSE,   FLA_NO_TRANSPOSE, FLA_ONE, X, A, FLA_ONE,  W );
    FLA_Scal( isgn, W );
/*
    if ( variant == 3 && type == FLA_ALG_UNBLOCKED )
    {
      FLA_Obj_show( "W", W, "%10.3e + %10.3e ", "" );
      FLA_Obj_show( "C_save", C_save, "%10.3e + %10.3e ", "" );
    }
*/
    FLA_Axpy( FLA_MINUS_ONE, C_save, W );
    FLA_Norm1( W, norm );
    FLA_Obj_extract_real_scalar( norm, diff );

    FLA_Obj_free( &X );
    FLA_Obj_free( &W );
  }

  *gflops = ( 2.0 / 3.0 ) * ( m * m * m ) / 
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( C ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( C_save, C );

  FLA_Obj_free( &C_save );
  FLA_Obj_free( &norm );
}

