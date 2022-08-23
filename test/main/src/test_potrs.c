/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

// Local prototypes.
void fla_test_potrs_experiment(test_params_t *params, integer datatype, integer  p_cur, integer  q_cur, integer  pci, integer  n_repeats,double* perf, double* time_min,double* residual);
void prepare_potrs_run(char* uplo, integer m, integer nrhs, void *A, integer datatype, void *b, integer n_repeats, double* time_min_);
void invoke_potrs(char* uplo, integer datatype, integer* m, void* A, integer* lda, integer *nrhs, void* b, integer* info);

void fla_test_potrs(test_params_t *params)
{
    char* op_str = "Cholesky factorization";
    char* front_str = "POTRS";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_potrs_experiment);
}

void fla_test_potrs_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double* perf,
    double* t,
    double* residual)
{
    integer n, info = 0, nrhs;
    void *A = NULL, *A_test = NULL;
    void *B = NULL, *X = NULL;
    void *B_test = NULL;
    double time_min = 1e9;
    char uplo = params->lin_solver_paramslist[pci].Uplo;
    nrhs = params->lin_solver_paramslist[pci].nrhs;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    /* Get input matrix dimensions. */
    n = p_cur;

    /* Create input matrix parameters */
    create_matrix(datatype, &A, n, n);
    create_matrix(datatype, &A_test, n, n);
    create_matrix(datatype, &B, n, nrhs);
    create_matrix(datatype, &X, n, nrhs);
    create_matrix(datatype, &B_test, n, nrhs);

    /* Initialize input symmetric positive definite matrix A */
    reset_matrix(datatype, n, n, A, n);
    rand_spd_matrix(datatype, &uplo, &A, n, n);
    copy_matrix(datatype, "full", n, n, A, n, A_test, n);
    /* cholesky factorisation of A as input to potrs */
    invoke_potrf(&uplo, datatype, &n, A, &n, &info);

    /* Generate random matrix B */
    rand_matrix(datatype, B, n, nrhs, n);
    copy_matrix(datatype, "full", n, nrhs, B, n, B_test, n);

    /* Invoke potrs API to find x using Ax-b */
    prepare_potrs_run(&uplo, n, nrhs, A, datatype, B_test, n_repeats, &time_min);
    copy_matrix(datatype, "full", n, nrhs, B_test, n, X, n);
    /* execution time */
    *t = time_min;
    /* Compute the performance of the best experiment repeat. */
    /* (2.0)m^2 flops for Ax=b computation. */
    *perf = (double)((2.0 * n * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* Validate potrs call by computing Ax-b */
    validate_potrs(n, nrhs, A_test, X, B, datatype, residual);

    free_matrix(A);
    free_matrix(A_test);
    free_matrix(B_test);
    free_matrix(B);
    free_matrix(X);

}

void prepare_potrs_run(char* uplo, 
    integer n,
    integer nrhs,
    void *A,
    integer datatype,
    void *B,
    integer n_repeats,
    double* time_min_)
{
    void *A_save = NULL, *B_test = NULL;
    double time_min = 1e9, exe_time;
    integer i, info = 0;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, n, n);
    create_matrix(datatype, &B_test, n, nrhs);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
        for each iteration */
        copy_matrix(datatype, "full", n, n, A, n, A_save, n);
        copy_matrix(datatype, "full", n, nrhs, B, n, B_test, n);
        exe_time = fla_test_clock();
        invoke_potrs(uplo, datatype, &n, A_save, &n, &nrhs, B_test, &info);
        exe_time = fla_test_clock() - exe_time;
        /* Get the best execution time */
        time_min = min(time_min, exe_time);
    }
    copy_matrix(datatype, "full", n, nrhs, B_test, n, B, n);
    *time_min_ = time_min;
    free_matrix(A_save);
    free_vector(B_test);
}

void invoke_potrs(char* uplo, integer datatype,
    integer* n,
    void* A,
    integer* lda,
    integer *nrhs,
    void* B,
    integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            spotrs_(uplo, n, nrhs, A, n, B, n, info);
            break;
        }
        case DOUBLE:
        {
            dpotrs_(uplo, n, nrhs, A, n, B, n, info);
            break;
        }
        case COMPLEX:
        {
            cpotrs_(uplo, n, nrhs, A, n, B, n, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zpotrs_(uplo, n, nrhs, A, n, B, n, info);
            break;
        }
    }
}
