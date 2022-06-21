/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_getrs_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_getrs_run(char *trans, integer m_A, integer n_A, void *A, void *B, integer* ipiv, integer datatype, integer n_repeats, double* time_min_);
void invoke_getrs(integer datatype, char *trans, integer *nrhs, integer *n, void *a, integer *lda, integer *ipiv, void *b, integer *ldb, integer *info);

void fla_test_getrs(test_params_t *params)
{
    char* op_str = "LU factorization";
    char* front_str = "GETRS";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_getrs_experiment);

}


void fla_test_getrs_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    double* perf,
    double* t,
    double* residual)
{
    integer m, cs_A;
    void* IPIV;
    void *A, *A_test, *B, *B_save, *X;
    double time_min = 1e9;
    char TRANS = params->lin_solver_paramslist[pci].transr;

    /* Determine the dimensions*/
    m = p_cur;
    cs_A = m;
    /* Create the matrices for the current operation*/
    create_matrix(datatype, &A, m, m);
    create_vector(INTEGER, &IPIV, m);
    create_vector(datatype, &B, m);
    create_vector(datatype, &B_save, m);
    create_vector(datatype, &X, m);
    /* Initialize the test matrices*/
    rand_matrix(datatype, A, m, m, cs_A);
    rand_vector(datatype, B, m, 1);
    /* Save the original matrix*/
    create_matrix(datatype, &A_test, m, m);
    copy_matrix(datatype, "full", m, m, A, cs_A, A_test, cs_A);
    copy_vector(datatype, m, B, 1, B_save, 1);
    /* call to API */
    prepare_getrs_run(&TRANS, m, m, A_test, B, IPIV, datatype, n_repeats, &time_min);
    copy_vector(datatype, m, B, 1, X, 1);
    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2mn^2 - (2/3)n^3 flops */
    *perf = (double)((2.0 * m * m * m) - ((2.0 / 3.0) * m * m * m)) / time_min / FLOPS_PER_UNIT_PERF;
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    validate_getrs(&TRANS, m, m, A, B_save, X, datatype, residual);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(IPIV);
    free_vector(B);
    free_vector(X);
    free_vector(B_save);
}


void prepare_getrs_run(char *TRANS,
    integer m_A,
    integer n_A,
    void* A,
    void* B,
    integer* IPIV,
    integer datatype,
    integer n_repeats,
    double* time_min_)
{
    integer cs_A;
    integer i;
    void *A_save, *B_test;
    integer info = 0, nrhs=1;
    double time_min = 1e9, exe_time;
    /* Get column stride */
    cs_A = m_A;
    /* Save the original matrix */
    create_matrix(datatype, &A_save, m_A, m_A);
    copy_matrix(datatype, "full", m_A, m_A, A, cs_A, A_save, cs_A);
    create_vector(datatype, &B_test, m_A);


    for (i = 0; i < n_repeats; ++i)
    {
        /* Copy original input data */
        copy_matrix(datatype, "full", m_A, m_A, A, cs_A, A_save, cs_A);
        copy_vector(datatype, m_A, B, 1, B_test, 1);

        exe_time = fla_test_clock();

        /*  call to API getrf to get AFACT */
        invoke_getrf(datatype, &m_A, &m_A, A_save, &cs_A, IPIV, &info);
        /*  call  getrs API with AFACT */
        invoke_getrs(datatype, TRANS, &m_A, &nrhs, A_save, &n_A, IPIV, B_test, &m_A, &info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);
        /*  Save the final result to B matrix*/
        copy_vector(datatype, m_A, B_test, 1, B, 1);
    }

    *time_min_ = time_min;

    free_matrix(A_save);
    free_vector(B_test);

}


/*
 *  GETRS_API calls LAPACK interface of
 *  Singular value decomposition - gesvd
 *  */
void invoke_getrs(integer datatype, char* trans, integer *n, integer *nrhs, void *a, integer *lda, integer *ipiv, void* b, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            dgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            cgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }
    }
}

