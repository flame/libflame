/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_potrf_experiment(test_params_t *params, integer datatype, integer  p_cur, integer  q_cur, integer  pci, integer  n_repeats,double* perf, double* time_min, double* residual);
void prepare_potrf_run(char* uplo, integer m, void *A, integer datatype, integer n_repeats, double* time_min_);

void fla_test_potrf(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Cholesky factorization";
    char* front_str = "POTRF";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_potrf_experiment);
}

void fla_test_potrf_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double* perf,
    double* time_min,
    double* residual)
{
    integer m;
    void *A = NULL, *A_test = NULL;
    char uplo = params->lin_solver_paramslist[pci].Uplo;

    /* Get input matrix dimensions */
    m = p_cur;

    /* Create input matrix parameters */
    create_matrix(datatype, &A, m, m);

    /* Initialize input symmetric positive definite matrix A */
    rand_spd_matrix(datatype, &uplo, &A, m, m);

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, &A_test, m, m);
    copy_matrix(datatype, "full", m, m, A, m, A_test, m);

    prepare_potrf_run(&uplo, m, A_test, datatype, n_repeats, time_min);

    /* Compute the performance of the best experiment repeat */
    /* (1/3)m^3 for real and (4/3)m^3 for complex*/
    *perf = (double)(1.0 / 3.0 * m * m * m) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    validate_potrf(&uplo, m, A, A_test, datatype, residual);

    free_matrix(A);
    free_matrix(A_test);
}

void prepare_potrf_run(char* uplo, integer m,
    void *A,
    integer datatype,
    integer n_repeats,
    double* time_min_)
{
    void *A_save = NULL;
    double time_min = 1e9, exe_time;
    integer i, info = 0;

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
        invoke_potrf(uplo, datatype, &m, A, &m, &info);
        exe_time = fla_test_clock() - exe_time;
        /* Get the best execution time */
        time_min = min(time_min, exe_time);
    }

    *time_min_ = time_min;
    free_matrix(A_save);
}

void invoke_potrf(char* uplo, integer datatype, integer* m, void* a, integer* lda, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            spotrf_(uplo, m, a, lda, info);
            break;
        }
        case DOUBLE:
        {
            dpotrf_(uplo, m, a, lda, info);
            break;
        }
        case COMPLEX:
        {
            cpotrf_(uplo, m, a, lda, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zpotrf_(uplo, m, a, lda, info);
            break;
        }
    }
}
