/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "test_lapack2flame.h"

int test_appq( FILE* stream, param_t param, result_t *result)
{
    FLA_Datatype datatype = param.datatype;
    FLA_Obj      A, B, C, t, x, y, z;
    double       time, time_min =  MAX_TIME_VALUE;
    unsigned int i,
             m = param.dims[0],
             n = param.dims[1],
             repeat = param.repeat;
    int          is_complex;

    // Create matrices.
    FLA_Obj_create( datatype, m, m, 0, 0, &A );
    FLA_Obj_create( datatype, m, n, 0, 0, &B );
    FLA_Obj_create( datatype, m, n, 0, 0, &C );
    FLA_Obj_create( datatype, m, 1, 0, 0, &t );

    FLA_Obj_create( datatype, n, 1, 0, 0, &x );
    FLA_Obj_create( datatype, m, 1, 0, 0, &y );
    FLA_Obj_create( datatype, m, 1, 0, 0, &z );

    // Initialize the test matrices.
    FLA_Random_matrix( A );
    FLA_Random_matrix( C );
    FLA_Random_matrix( x );
    FLA_QR_blk_external( A, t );

    // Repeat the experiment repeat times and record results.
    for ( i = 0; i < repeat; ++i )
    {
        FLA_Copy( C, B );
        time = FLA_Clock();
        FLA_Apply_Q_blk_external( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_COLUMNWISE, A, t, B );
        time = FLA_Clock() - time;

        time_min = min( time_min, time );
    }

    is_complex = FLA_Obj_is_complex( A );
    result->performance = -time_min;

    // y := beta y + alpha B x   (beta = 0, alapha =1)
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, B, x, FLA_ZERO, y );

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, C, x, FLA_ZERO, z );
    FLA_Apply_Q_blk_external( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_COLUMNWISE, A, t, z );

    result->residual = FLA_Max_elemwise_diff( y, z );

    FLA_Obj_free( &A );
    FLA_Obj_free( &B );
    FLA_Obj_free( &C );
    FLA_Obj_free( &x );
    FLA_Obj_free( &y );
    FLA_Obj_free( &z );
    FLA_Obj_free( &t );

    return 0;
}

