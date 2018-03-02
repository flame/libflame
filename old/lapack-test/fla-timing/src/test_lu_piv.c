/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "test_lapack2flame.h"

#define FMULS FMULS_GETRF(m,n)
#define FADDS FADDS_GETRF(m,n)

int test_lu_piv( FILE* stream, param_t param, result_t *result)
{
    FLA_Datatype datatype = param.datatype;
    FLA_Obj      A, B, p, x, y;
    double       time, time_min =  MAX_TIME_VALUE;
    unsigned int i,
             m = param.dims[0],
             n = param.dims[1],
             repeat = param.repeat;
    int          is_complex;

    // Create matrices.
    FLA_Obj_create( datatype, m, n, 0, 0, &A );
    FLA_Obj_create( datatype, m, n, 0, 0, &B );
    FLA_Obj_create( datatype, n, 1, 0, 0, &x );
    FLA_Obj_create( datatype, m, 1, 0, 0, &y );
    FLA_Obj_create( FLA_INT , min(m,n), 1, 0, 0, &p );

    // Initialize the test matrices.
    FLA_Random_matrix( B );
    FLA_Random_matrix( x );

    // Repeat the experiment repeat times and record results.
    for ( i = 0; i < repeat; ++i )
    {
        FLA_Copy( B, A );
        time = FLA_Clock();

        FLA_LU_piv_blk_external( A, p );
        time = FLA_Clock() - time;

        time_min = min( time_min, time );
    }

    is_complex = FLA_Obj_is_complex( A );
    result->performance = ( FMULS * FP_PER_MUL(is_complex) +
                            FADDS * FP_PER_ADD(is_complex) )/time_min/FLOPS_PER_UNIT_PERF;

    // y := beta y + alpha B x   (beta = 0, alapha =1)
    if (m == n)
    {
        FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, B, x, FLA_ZERO, y );

        FLA_Trmv_external( FLA_UPPER_TRIANGULAR , FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A, x );
        FLA_Trmv_external( FLA_LOWER_TRIANGULAR , FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,    A, x );
        FLA_Apply_pivots_unb_external( FLA_LEFT, FLA_NO_TRANSPOSE, p, y );

        result->residual = FLA_Max_elemwise_diff( x, y );

    }
    else
    {
        result->residual = -1.0;
    }

    FLA_Obj_free( &A );
    FLA_Obj_free( &B );
    FLA_Obj_free( &x );
    FLA_Obj_free( &y );
    FLA_Obj_free( &p );

    return 0;
}

