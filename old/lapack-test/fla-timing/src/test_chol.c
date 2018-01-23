/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "test_lapack2flame.h"

#define FMULS FMULS_POTRF(m)
#define FADDS FADDS_POTRF(m)

int test_chol( FILE* stream, param_t param, result_t *result)
{
    FLA_Datatype datatype = param.datatype;
    FLA_Trans    transa, transb;
    FLA_Uplo     uplo = param.uplos[0];
    FLA_Obj      A, B, x, y;
    double       time, time_min =  MAX_TIME_VALUE;
    unsigned int i,
             m = param.dims[0],
             repeat = param.repeat;
    int          is_complex;

    // Create matrices.
    FLA_Obj_create( datatype, m, m, 0, 0, &A );
    FLA_Obj_create( datatype, m, m, 0, 0, &B );
    FLA_Obj_create( datatype, m, 1, 0, 0, &x );
    FLA_Obj_create( datatype, m, 1, 0, 0, &y );

    // Initialize the test matrices.
    FLA_Random_spd_matrix( uplo, B );
    FLA_Triangularize( uplo, FLA_NONUNIT_DIAG, B );
    FLA_Random_matrix( x );

    // Repeat the experiment repeat times and record results.
    for ( i = 0; i < repeat; ++i )
    {
        FLA_Copy( B, A );
        time = FLA_Clock();
        FLA_Chol_blk_external( uplo, A );
        time = FLA_Clock() - time;

        time_min = min( time_min, time );
    }

    is_complex = FLA_Obj_is_complex( A );
    result->performance = ( FMULS * FP_PER_MUL(is_complex) +
                            FADDS * FP_PER_ADD(is_complex) )/time_min/FLOPS_PER_UNIT_PERF;


    // y := beta y + alpha B x   (beta = 0, alapha =1)
    FLA_Hermitianize( uplo, B );
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, B, x, FLA_ZERO, y );

    if (uplo == FLA_LOWER_TRIANGULAR)
    {
        // x := L ( L^T x )
        transa = FLA_CONJ_TRANSPOSE;
        transb = FLA_NO_TRANSPOSE;
    }
    else if (uplo == FLA_UPPER_TRIANGULAR)
    {
        // x := U^T ( U x )
        transa = FLA_NO_TRANSPOSE;
        transb = FLA_CONJ_TRANSPOSE;
    }
    else
    {
        fprintf(stderr, "Error not a valid uplo\n");
        CHKERR(-1);
    }
    FLA_Trmv_external( uplo , transa, FLA_NONUNIT_DIAG, A, x );
    FLA_Trmv_external( uplo , transb, FLA_NONUNIT_DIAG, A, x );

    result->residual = FLA_Max_elemwise_diff( x, y );

    FLA_Obj_free( &A );
    FLA_Obj_free( &B );
    FLA_Obj_free( &x );
    FLA_Obj_free( &y );

    return 0;
}

