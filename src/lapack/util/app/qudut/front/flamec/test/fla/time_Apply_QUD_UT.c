/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


void time_Apply_QUD_UT(
                 integer n_repeats, integer mB, integer mC, integer mD, integer n, integer n_rhs, integer b_alg,
                 FLA_Obj R_BC, FLA_Obj R_BD, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W,
                 FLA_Obj bR_BC, FLA_Obj bR_BD, FLA_Obj bC, FLA_Obj bD,
                 double *dtime, double *diff, double *gflops );

void time_Apply_QUD_UT(
                 integer n_repeats, integer mB, integer mC, integer mD, integer n, integer n_rhs, integer b_alg,
                 FLA_Obj R_BC, FLA_Obj R_BD, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W,
                 FLA_Obj bR_BC, FLA_Obj bR_BD, FLA_Obj bC, FLA_Obj bD,
                 double *dtime, double *diff, double *gflops )
{
  integer
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    bR_BD_save, bC_save, bD_save;

  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, bR_BD, &bR_BD_save );
  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, bC, &bC_save );
  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, bD, &bD_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( bR_BD_save, bR_BD );
    FLA_Copy_external( bC_save, bC );
    FLA_Copy_external( bD_save, bD );

    *dtime = FLA_Clock();

    FLA_Apply_QUD_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                      T, W,
                         bR_BD,
                      C, bC,
                      D, bD );

    *dtime = FLA_Clock() - *dtime;
    dtime_old = fla_min( *dtime, dtime_old );

  }

  {

    FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, R_BD, bR_BD );

    FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, R_BC, bR_BC );

    *diff = FLA_Max_elemwise_diff( bR_BD, bR_BC );
  }

  *gflops = n * n_rhs * ( 2.0 * mC + 2.0 * mD + 0.5 * b_alg + 0.5 ) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( R_BD ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Obj_free( &bR_BD_save );
  FLA_Obj_free( &bC_save );
  FLA_Obj_free( &bD_save );
}

