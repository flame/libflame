/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_gesv_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_gesv_run(integer n_A, integer nrhs, void *A, void *B, integer* ipiv, integer datatype, integer n_repeats, double* time_min_);
void invoke_gesv(integer datatype, integer *nrhs, integer *n, void *a, integer *lda, integer *ipiv, void *b, integer *ldb, integer *info);

void fla_test_gesv(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Linear Solve using LU";
    char* front_str = "GESV";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_gesv_experiment);

}


void fla_test_gesv_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    double* perf,
    double* t,
    double* residual)
{
    integer n, cs_A, NRHS;
    void* IPIV;
    void *A, *A_save, *B, *B_save;
    double time_min = 1e9;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    NRHS = params->lin_solver_paramslist[pci].nrhs;
    /* Determine the dimensions*/
    n = p_cur;
    cs_A = n;
    /* Create the matrices for the current operation*/
    create_matrix(datatype, &A, n, n);
    create_matrix(datatype, &A_save, n, n);
    create_vector(INTEGER, &IPIV, n);
    create_matrix(datatype, &B, n, NRHS);
    create_matrix(datatype, &B_save, n, NRHS);
    /* Initialize the test matrices*/
    rand_matrix(datatype, A, n, n, cs_A);
    rand_matrix(datatype, B, n, NRHS, cs_A);
    /* Save the original matrix*/
    copy_matrix(datatype, "full", n, n, A, cs_A, A_save, cs_A);
    copy_matrix(datatype, "full", n, NRHS, B, cs_A, B_save, cs_A);
    /* call to API */
    prepare_gesv_run(n, NRHS, A_save, B_save, IPIV, datatype, n_repeats, &time_min);
    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2mn^2 - (2/3)n^3 flops */
    *perf = (double)((2.0 * n * n *n) - ((2.0 / 3.0) *n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    validate_gesv(n, NRHS, A, B, B_save, datatype, residual);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_save);
    free_vector(IPIV);
    free_matrix(B);
    free_matrix(B_save);
}


void prepare_gesv_run(integer n_A,
    integer nrhs,
    void* A,
    void* B,
    integer* IPIV,
    integer datatype,
    integer n_repeats,
    double* time_min_)
{
    integer cs_A;
    integer i;
    void *A_test, *B_test;
    integer info = 0;
    double time_min = 1e9, exe_time;
    /* Get column stride */
    cs_A = n_A;
    /* Save the original matrix */
    create_matrix(datatype, &A_test, n_A, n_A);
    create_matrix(datatype, &B_test, n_A, nrhs);


    for (i = 0; i < n_repeats; ++i)
    {

        /* Copy original input data */
        copy_matrix(datatype, "full", n_A, n_A, A, cs_A, A_test, cs_A);
        copy_matrix(datatype, "full", n_A, nrhs, B, cs_A, B_test, cs_A);

        exe_time = fla_test_clock();

        /*  call  gesv API  */
        invoke_gesv(datatype, &n_A, &nrhs, A_test, &n_A, IPIV, B_test, &n_A, &info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

    }

    *time_min_ = time_min;
    /*  Save the final result to B matrix*/
    copy_matrix(datatype, "full", n_A, nrhs, B_test, cs_A, B, cs_A);

    free_matrix(A_test);
    free_matrix(B_test);

}


/*
 *  gesv_API calls LAPACK interface of
 *  Singular value decomposition - gesvd
 *  */
void invoke_gesv(integer datatype, integer *n, integer *nrhs, void *a, integer *lda, integer *ipiv, void* b, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }
        
        case DOUBLE:
        {
            dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            cgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }
    }
}

