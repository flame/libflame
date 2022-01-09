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
#define FLA_ALG_UNB_ASM       3
#define FLA_ALG_BLOCKED       4

void time_Apply_G_rf(
               integer variant, integer type, integer n_repeats, integer m, integer k, integer n, integer b_alg,
               FLA_Obj A, FLA_Obj A_ref, FLA_Obj G, FLA_Obj P,
               double *dtime, double *diff, double *gflops );

void time_Apply_G_rf(
               integer variant, integer type, integer n_repeats, integer m, integer k, integer n, integer b_alg,
               FLA_Obj A, FLA_Obj A_ref, FLA_Obj G, FLA_Obj P,
               double *dtime, double *diff, double *gflops )
{
  integer irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, G_save, norm;

  if ( FLA_Obj_is_real( A ) )
  {
    if (
       //( variant == 1 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 1 && type == FLA_ALG_UNB_ASM ) ||
       //( variant == 1 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 2 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 2 && type == FLA_ALG_UNB_ASM ) ||
       //( variant == 2 && type == FLA_ALG_BLOCKED ) ||

       //( variant == 3 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 3 && type == FLA_ALG_UNB_ASM ) ||
       //( variant == 3 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 6 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 6 && type == FLA_ALG_UNB_ASM ) ||
       //( variant == 6 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 9 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 9 && type == FLA_ALG_UNB_ASM ) ||
       //( variant == 9 && type == FLA_ALG_BLOCKED ) ||

       ( variant == 4 ) ||
       ( variant == 5 ) ||
       ( variant == 7 ) ||
       ( variant == 8 ) ||
       FALSE
    ) 
    {
      *gflops = 0.0;
      *diff   = 0.0;
      return;
    }
  }
  else if ( FLA_Obj_is_complex( A ) )
  {
    if (
       //( variant == 1 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 1 && type == FLA_ALG_UNB_ASM ) ||
       //( variant == 1 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 2 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 2 && type == FLA_ALG_UNB_ASM ) ||
       //( variant == 2 && type == FLA_ALG_BLOCKED ) ||

       //( variant == 3 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 3 && type == FLA_ALG_UNB_ASM ) ||
       //( variant == 3 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 6 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 6 && type == FLA_ALG_UNB_ASM ) ||
       //( variant == 6 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 9 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 9 && type == FLA_ALG_UNB_ASM ) ||
       //( variant == 9 && type == FLA_ALG_BLOCKED ) ||

       ( variant == 4 ) ||
       ( variant == 5 ) ||
       ( variant == 7 ) ||
       ( variant == 8 ) ||
       FALSE
    )
    {
      *gflops = 0.0;
      *diff   = 0.0;
      return;
    }
  }

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, G, &G_save );
  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

  //dim_t b_flash_m = b_alg;
  //dim_t b_flash_n = n;
  //FLASH_Obj_create_hier_copy_of_flat_ext( A, 1, &b_flash_m, &b_flash_n, &AH ); 

//printf ( "flash dims: %d x %d\n", FLA_Obj_length( AH ), FLA_Obj_width( AH ) );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( G, G_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );
    FLA_Copy_external( G_save, G );
    //FLASH_Obj_hierarchify( A_save, AH );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:
      break;

    // Time variant 1
    case 1:
    {
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Apply_G_rf_opt_var1( G, A );
        break;
      case FLA_ALG_UNB_ASM:
        FLA_Apply_G_rf_asm_var1( G, A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Apply_G_rf_blk_var1( G, A, b_alg );
        break;
      }
      break;
    }

    // Time variant 2
    case 2:
    {
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Apply_G_rf_opt_var2( G, A );
        break;
      case FLA_ALG_UNB_ASM:
        FLA_Apply_G_rf_asm_var2( G, A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Apply_G_rf_blk_var2( G, A, b_alg );
        break;
      }
      break;
    }

    // Time variant 3
    case 3:
    {
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Apply_G_rf_opt_var3( G, A );
        break;
      case FLA_ALG_UNB_ASM:
        FLA_Apply_G_rf_asm_var3( G, A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Apply_G_rf_blk_var3( G, A, b_alg );
        break;
      }
      break;
    }

    // Time variant 6
    case 6:
    {
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Apply_G_rf_opt_var6( G, A );
        break;
      case FLA_ALG_UNB_ASM:
        FLA_Apply_G_rf_asm_var6( G, A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Apply_G_rf_blk_var6( G, A, b_alg );
        break;
      }
      break;
    }

    // Time variant 9
    case 9:
    {
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Apply_G_rf_opt_var9( G, A );
        break;
      case FLA_ALG_UNB_ASM:
        FLA_Apply_G_rf_asm_var9( G, A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Apply_G_rf_blk_var9( G, A, b_alg );
        break;
      }
      break;
    }


    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );

  }

  if ( variant == 1 && type == FLA_ALG_UNB_OPT )
  {
    //FLA_Obj_show( "A_ref", A, "%9.2e + %9.2e ", "" );
    //FLA_Obj_show( "A", A, "%9.2e ", "" );

    FLA_Copy( A, A_ref );
    *diff = 0.0;
  }
  else
  {
    //FLA_Obj_show( "A", A, "%9.2e + %9.2e ", "" );

//if ( variant == 7 && type == FLA_ALG_UNB_ASM )
    //FLA_Obj_show( "A", A, "%9.2e", "" );

    //if ( variant == 9 ) FLASH_Obj_flatten( AH, A );

    FLA_Axpy( FLA_MINUS_ONE, A_ref, A );
    FLA_Norm_frob( A, norm );
    FLA_Obj_extract_real_scalar( norm, diff );

    //*diff = FLA_Max_elemwise_diff( A_ref, A );
  }


  *gflops = 6.0 * k * m * ( n - 1 ) /
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 2.0;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( G_save, G );

  //FLASH_Obj_free( &AH );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &G_save );
  FLA_Obj_free( &norm );
}

