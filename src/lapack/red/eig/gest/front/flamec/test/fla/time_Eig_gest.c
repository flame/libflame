/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Eig_gest( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B );

void time_Eig_gest(
                integer param_combo, integer type, integer n_repeats, integer n,
                FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj norm,
                double *dtime, double *diff, double *gflops );


void time_Eig_gest(
                integer param_combo, integer type, integer nrepeats, integer m,
                FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj norm,
                double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_save = 1.0e9;

  FLA_Obj
    A_save, B_save;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B, &B_save );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( B, B_save );


  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLA_Copy_external( A_save, A );
    FLA_Copy_external( B_save, B );

    *dtime = FLA_Clock();

    switch( param_combo ){

    case 0:
    {
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Eig_gest( FLA_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
        break;
      case FLA_ALG_FRONT:
        FLA_Eig_gest( FLA_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
        break;
      }
      break;
    }

    case 1:
    {
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Eig_gest( FLA_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
        break;
      case FLA_ALG_FRONT:
        FLA_Eig_gest( FLA_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
        break;
      }
      break;
    }

    case 2:
    {
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Eig_gest( FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
        break;
      case FLA_ALG_FRONT:
        FLA_Eig_gest( FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
        break;
      }
      break;
    }

    case 3:
    {
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Eig_gest( FLA_NO_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
        break;
      case FLA_ALG_FRONT:
        FLA_Eig_gest( FLA_NO_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
        break;
      }
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_save = fla_min( *dtime, dtime_save );
  }

  FLA_Hermitianize( uplo, A );

  // Recover A.
  if ( inv == FLA_NO_INVERSE )
  {
    if ( uplo == FLA_LOWER_TRIANGULAR )
    {
      // A = L' * A_orig * L
      // A_orig = inv(L') * A * inv(L)
      FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A );
      FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
      FLA_Trsm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
    }
    else // if ( uplo == FLA_UPPER_TRIANGULAR )
    {
      // A = U * A_orig * U'
      // A_orig = inv(U) * A * inv(U')
      FLA_Hermitianize( FLA_UPPER_TRIANGULAR, A );
      FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
      FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
    }
  }
  else // if ( inv == FLA_INVERSE )
  {
    if ( uplo == FLA_LOWER_TRIANGULAR )
    {
      // A = inv(L) * A_orig * inv(L')
      // A_orig = L * A * L'
      FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A );
      FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
      FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
    }
    else // if ( uplo == FLA_UPPER_TRIANGULAR )
    {
      // A = inv(U') * A_orig * inv(U)
      // A_orig = U' * A * U
      FLA_Hermitianize( FLA_UPPER_TRIANGULAR, A );
      FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
      FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B, A );
    }
  }

  *diff = FLA_Max_elemwise_diff( A, A_save );
/*
  FLA_Axpy_external( FLA_MINUS_ONE, A_save, A );
  FLA_Norm1( A, norm );
  FLA_Obj_extract_real_scalar( norm, diff );
*/
  *gflops = 1.0 * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) / 
            dtime_save / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_save;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( B_save, B );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &B_save );
}

