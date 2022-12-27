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

FLA_Error REF_Hevd_ln( FLA_Obj A, FLA_Obj l );
FLA_Error REF_Hevdd_ln( FLA_Obj A, FLA_Obj l );
void time_Hevd_ln(
               integer variant, integer type, integer n_repeats, integer m, integer b_alg,
               FLA_Obj A, FLA_Obj l,
               double *dtime, double *diff, double *gflops );


void time_Hevd_ln(
               integer variant, integer type, integer n_repeats, integer m, integer b_alg,
               FLA_Obj A, FLA_Obj l,
               double *dtime, double *diff, double *gflops )
{
  integer irep;

  double
    k, dtime_old = 1.0e9;

  FLA_Obj
    A_save;

  if (
       //( variant == 1 && type == FLA_ALG_UNBLOCKED ) ||
       //( variant == 1 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 1 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       //( variant == 1 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 2 && type == FLA_ALG_UNBLOCKED ) ||
       //( variant == 2 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 2 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       //( variant == 2 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 3 && type == FLA_ALG_UNBLOCKED ) ||
       //( variant == 3 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 3 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       //( variant == 3 && type == FLA_ALG_BLOCKED ||
       FALSE
     )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );

  FLA_Copy_external( A, A_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:
      REF_Hevd_ln( A, l );
      break;

    case -1:
      REF_Hevdd_ln( A, l );
      break;

    // Time variant 1
    case 1:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Hevd_ln_unb_var1( A, l );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Hevd_ln_unb_var1( A, l );
        break;
      }
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = fla_min( *dtime, dtime_old );

  }

  {
    FLA_Obj Av, lv, A_rev_evd, norm;

    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_rev_evd ); 
    FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A_save, &Av ); 
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, l, &lv ); 
    FLA_Hevd_lv_unb_var1( Av, lv );

//FLA_Obj_show( "l", l, "%9.2e", "" );
//FLA_Obj_show( "lv", lv, "%9.2e", "" );
//FLA_Obj_show( "V", Av, "%9.2e + %9.2e ", "" );

    FLA_Axpy( FLA_MINUS_ONE, l, lv );
    FLA_Norm_frob( lv, norm );
    FLA_Obj_extract_real_scalar( norm, diff );

    FLA_Obj_free( &Av );
    FLA_Obj_free( &lv );
    FLA_Obj_free( &A_rev_evd );
    FLA_Obj_free( &norm );
  }

  k = 2.05;
  *gflops = ( 4.0 / 3.0 * m * m * m + 
              9.0       * k * m * m ) / 
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );

  FLA_Obj_free( &A_save );
}

