
#include "test_lapack2flame.h"

int test_svd( FILE* stream, param_t param, result_t *result)
{
    FLA_Datatype datatype = param.datatype;
    FLA_Obj      A, B, s, U, V, x, y, w, z;
    double       time, time_min =  MAX_TIME_VALUE;
    unsigned int i,
             m = param.dims[0],
             n = param.dims[1],
             repeat = param.repeat;
    int          is_complex;

    // Create matrices.
    FLA_Obj_create( datatype, m, n, 0, 0, &A );
    FLA_Obj_create( datatype, m, n, 0, 0, &B );

    FLA_Obj_create( datatype,        m, m, 0, 0, &U );
    FLA_Obj_create( datatype,        n, n, 0, 0, &V );
    FLA_Obj_create( datatype, min(m,n), 1, 0, 0, &s );

    FLA_Obj_create( datatype,        n, 1, 0, 0, &x );
    FLA_Obj_create( datatype,        m, 1, 0, 0, &y );
    FLA_Obj_create( datatype, min(m,n), 1, 0, 0, &w );
    FLA_Obj_create( datatype,        m, 1, 0, 0, &z );

    // Initialize the test matrices.
    FLA_Random_matrix( B );
    FLA_Random_matrix( x );

    // Repeat the experiment repeat times and record results.
    for ( i = 0; i < repeat; ++i )
    {
        FLA_Copy( B, A );
        time = FLA_Clock();
        FLA_Svd_external( FLA_SVD_VECTORS_ALL, FLA_SVD_VECTORS_ALL,
                          A, s, U, V );
        time = FLA_Clock() - time;

        time_min = min( time_min, time );
    }

    is_complex = FLA_Obj_is_complex( A );
    result->performance = -time_min;

    // y := beta y + alpha B x   (beta = 0, alapha =1)
    if (m == n)
    {
        FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, B, x, FLA_ZERO, y );

        FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, V, x, FLA_ZERO, w );
        FLA_Apply_diag_matrix( FLA_LEFT,  FLA_NO_CONJUGATE, s, w );
        FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, U, w, FLA_ZERO, z );

        result->residual = FLA_Max_elemwise_diff( y, z );

    }
    else
    {
        result->residual = -1.0;
    }

    FLA_Obj_free( &A );
    FLA_Obj_free( &B );
    FLA_Obj_free( &U );
    FLA_Obj_free( &V );
    FLA_Obj_free( &s );
    FLA_Obj_free( &x );
    FLA_Obj_free( &y );
    FLA_Obj_free( &w );
    FLA_Obj_free( &z );

    return 0;
}

