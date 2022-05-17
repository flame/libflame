/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

#define NUM_PARAM_COMBOS 2
#define NUM_MATRIX_ARGS  1


// Local prototypes.
void fla_test_potrs_experiment(test_params_t *params, integer datatype, integer  p_cur, integer  q_cur, integer  pci, integer  n_repeats,double* perf, double* time_min,double* residual);
void prepare_potrs_run(char* uplo, integer m, void *A, integer datatype, void *b, integer n_repeats, double* time_min_);
void invoke_potrs(char* uplo, integer datatype, integer* m, void* A, integer* lda, integer *nrhs, void* b, integer* info);

void fla_test_potrs(test_params_t *params)
{
    char* op_str = "Cholesky factorization";
    char* front_str = "POTRS";
    char* lapack_str = "LAPACK";
    char* pc_str[NUM_PARAM_COMBOS] = {"U","L"};

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS,
                            params, LIN, fla_test_potrs_experiment);
}

void fla_test_potrs_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double* perf,
    double* time_min,
    double* residual)
{
    integer m, info = 0;
    void *A = NULL, *A_test = NULL;
    void *b = NULL, *x = NULL;
    void *b_test = NULL;
    char uplo = params->lin_solver_paramslist[pci].Uplo;

    /* Get input matrix dimensions. */
    m = p_cur;

    /* Create input matrix parameters */
    create_matrix(datatype, &A, m, m);

    /* Initialize input symmetric positive definite matrix A */
    rand_spd_matrix(datatype, &uplo, &A, m, m);

    /* Make a copy of input matrix A. This is required to validate the API functionality.*/
    create_matrix(datatype, &A_test, m, m);

    /* Create vectors to compute Ax-b */
    create_vector(datatype, &b, m);
    create_vector(datatype, &x, m);
    create_vector(datatype, &b_test, m);

    copy_matrix(datatype, "full", m, m, A, m, A_test, m);

    /* cholesky factorisation of A as input to potrs */
    invoke_potrf(&uplo, datatype, &m, A_test, &m, &info);

    /* Generate random vector b */
    rand_vector(datatype, b, m, 1);
    copy_vector(datatype, m, b, 1, b_test, 1);

    /* Invoke potrs API to find x using Ax-b */
    prepare_potrs_run(&uplo, m, A_test, datatype, b_test, n_repeats, time_min);
    copy_vector(datatype, m, b_test, 1, x, 1);

    /* Compute the performance of the best experiment repeat. */
    /* (2.0)m^2 flops for Ax=b computation. */
    *perf = (double)(2.0 * m * m) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* Validate potrs call by computing Ax-b */
    validate_potrs(&uplo, m, A, A_test, datatype, x, b, residual);

    free_matrix(A);
    free_matrix(A_test);
    free_matrix(b_test);
    free_matrix(b);
    free_matrix(x);
}

void prepare_potrs_run(char* uplo, integer m,
    void *A,
    integer datatype,
    void *b,
    integer n_repeats,
    double* time_min_)
{
    void *A_save = NULL;
    double time_min = 1e9, exe_time;
    integer i, nrhs = 1, info = 0;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, m, m);
    copy_matrix(datatype, "full", m, m, A, m, A_save, m);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
        for each iteration */
        copy_matrix(datatype, "full", m, m, A_save, m, A, m);
        exe_time = fla_test_clock();
        invoke_potrs(uplo, datatype, &m, A, &m, &nrhs, b, &info);
        exe_time = fla_test_clock() - exe_time;
        /* Get the best execution time */
        time_min = min(time_min, exe_time);
    }

    *time_min_ = time_min;
    free_matrix(A_save);
}

void invoke_potrs(char* uplo, integer datatype,
    integer* m,
    void* A,
    integer* lda,
    integer *nrhs,
    void* b,
    integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            spotrs_(uplo, m, nrhs, A, m, b, m, info);
            break;
        }
        case DOUBLE:
        {
            dpotrs_(uplo, m, nrhs, A, m, b, m, info);
            break;
        }
        case COMPLEX:
        {
            cpotrs_(uplo, m, nrhs, A, m, b, m, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zpotrs_(uplo, m, nrhs, A, m, b, m, info);
            break;
        }
    }
}
