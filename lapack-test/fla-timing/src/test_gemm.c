
#include "test_lapack2flame.h"

#define FMULS FMULS_GEMM(m, n, k)
#define FADDS FADDS_GEMM(m, n, k)

int test_gemm( FILE* stream, param_t param, result_t *result)
{
    FLA_Datatype datatype = param.datatype;
    FLA_Trans    transa = param.trans[0], transb = param.trans[1];
    FLA_Obj      A, B, C, x, y, z, w;
    FLA_Obj      alpha, beta;
    double       time, time_min =  MAX_TIME_VALUE;
    unsigned int i,
             m = param.dims[0],
             n = param.dims[1],
             k = param.dims[2],
             repeat = param.repeat;
    int          is_trans, is_complex;

    // Create matrices.
    is_trans = (transa == FLA_NO_TRANSPOSE);
    FLA_Obj_create( datatype, (is_trans ? m:k), (is_trans ? k:m), 0,0, &A );
    is_trans = (transb == FLA_NO_TRANSPOSE);
    FLA_Obj_create( datatype, (is_trans ? k:n), (is_trans ? n:k), 0,0, &B );
    FLA_Obj_create( datatype,                m,                n, 0,0, &C );

    FLA_Obj_create( datatype, n, 1, 0, 0, &x );
    FLA_Obj_create( datatype, m, 1, 0, 0, &y );
    FLA_Obj_create( datatype, m, 1, 0, 0, &z );
    FLA_Obj_create( datatype, k, 1, 0, 0, &w );

    // Initialize the test matrices.
    FLA_Random_matrix( A );
    FLA_Random_matrix( B );
    FLA_Random_matrix( C );

    FLA_Random_matrix( x );
    FLA_Set( FLA_ZERO, y );
    FLA_Set( FLA_ZERO, w );
    FLA_Set( FLA_ZERO, z );

    // Constants.
    alpha = FLA_MINUS_ONE;
    beta  = FLA_ZERO;

    // Repeat the experiment repeat times and record results.
    for ( i = 0; i < repeat; ++i )
    {
        time = FLA_Clock();
        FLA_Gemm_external( transa, transb, alpha, A, B, beta, C );
        time = FLA_Clock() - time;

        time_min = min( time_min, time );
    }

    is_complex = FLA_Obj_is_complex( C );
    result->performance = ( FMULS * FP_PER_MUL(is_complex) +
                            FADDS * FP_PER_ADD(is_complex) )/time_min/FLOPS_PER_UNIT_PERF;

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, C, x, FLA_ZERO, y );
    FLA_Gemv_external( transb,           FLA_ONE, B, x, FLA_ZERO, w );
    FLA_Gemv_external( transa,           alpha,   A, w, FLA_ZERO, z );

    result->residual = FLA_Max_elemwise_diff( y, z );

    FLA_Obj_free( &A );
    FLA_Obj_free( &B );
    FLA_Obj_free( &C );
    FLA_Obj_free( &x );
    FLA_Obj_free( &y );
    FLA_Obj_free( &z );
    FLA_Obj_free( &w );

    return 0;
}

