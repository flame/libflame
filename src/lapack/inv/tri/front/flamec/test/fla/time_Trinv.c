/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Trinv( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );

void time_Trinv(
                 integer param_combo, integer type, integer nrepeats, integer m, FLA_Uplo uplo, FLA_Diag diag,
                 FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj norm,
                 double *dtime, double *diff, double *gflops );


void time_Trinv(
                 integer param_combo, integer type, integer nrepeats, integer m, FLA_Uplo uplo, FLA_Diag diag,
                 FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj norm,
                 double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, b_save, b_orig_save;

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

    switch( param_combo ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trinv( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
        break;
      case FLA_ALG_FRONT:
        FLA_Trinv( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 1:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trinv( FLA_LOWER_TRIANGULAR, FLA_UNIT_DIAG, A );
        break;
      case FLA_ALG_FRONT:
        FLA_Trinv( FLA_LOWER_TRIANGULAR, FLA_UNIT_DIAG, A );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trinv( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
        break;
      case FLA_ALG_FRONT:
        FLA_Trinv( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 3:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trinv( FLA_UPPER_TRIANGULAR, FLA_UNIT_DIAG, A );
        break;
      case FLA_ALG_FRONT:
        FLA_Trinv( FLA_UPPER_TRIANGULAR, FLA_UNIT_DIAG, A );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = fla_min( *dtime, dtime_old );
  }

  {
    FLA_Trmv_external( uplo, FLA_NO_TRANSPOSE,
                       diag, A, b );

    FLA_Trmv_external( uplo, FLA_NO_TRANSPOSE,
                       diag, A_save, b );

    FLA_Axpy_external( FLA_MINUS_ONE, b_orig, b );

    FLA_Nrm2_external( b, norm );
    FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, norm,
                               1, 1, diff, 1, 1 );
  }

  *gflops = 1.0 / 4.0 *
            m * m * m /
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

