/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_LU_piv( FLA_Obj A, FLA_Obj p );
void time_LU(
              integer is_pivoting, integer type, integer nrepeats, integer m, integer n,
              FLA_Obj A, FLA_Obj p, FLA_Obj b, FLA_Obj b_ref, FLA_Obj b_norm,
              double *dtime, double *diff, double *gflops );


void time_LU(
              integer is_pivoting, integer type, integer nrepeats, integer m, integer n,
              FLA_Obj A, FLA_Obj p, FLA_Obj b, FLA_Obj b_ref, FLA_Obj b_norm,
              double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, p_save, b_save, b_ref_save;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, p, &p_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b, &b_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b_ref, &b_ref_save );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( p, p_save );
  FLA_Copy_external( b, b_save );
  FLA_Copy_external( b_ref, b_ref_save );


  for ( irep = 0 ; irep < nrepeats; irep++ ){

    FLA_Copy_external( A_save, A );
    FLA_Copy_external( p_save, p );
    FLA_Copy_external( b_save, b );
    FLA_Copy_external( b_ref_save, b_ref );

    *dtime = FLA_Clock();

    switch( is_pivoting ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_LU_piv( A, p );
        break;
      case FLA_ALG_FRONT:
        FLA_LU_piv( A, p );
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

  if ( type == FLA_ALG_REFERENCE )
  {
    FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, b );

    FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_UNIT_DIAG, A, b );
    FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, A, b );

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE,
                       A_save, b, FLA_MINUS_ONE, b_ref );

    FLA_Nrm2_external( b_ref, b_norm );
    FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, b_norm,
                               1, 1, diff, 1, 1 );
  }
  else
  {
    FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, b );

    FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_UNIT_DIAG, A, b );
    FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, A, b );

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE,
                       A_save, b, FLA_MINUS_ONE, b_ref );

    FLA_Nrm2_external( b_ref, b_norm );
    FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, b_norm,
                               1, 1, diff, 1, 1 );
  }

  *gflops = 2.0 / 3.0 * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) / 
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( p_save, p );
  FLA_Copy_external( b_save, b );
  FLA_Copy_external( b_ref_save, b_ref );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &p_save );
  FLA_Obj_free( &b_save );
  FLA_Obj_free( &b_ref_save );
}

