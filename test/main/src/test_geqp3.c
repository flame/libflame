/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_geqp3_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                               integer pci, integer n_repeats, double* perf, double* t,double* residual);
void prepare_geqp3_run(integer m_A, integer n_A, void *A, integer *jpvt, void *T, integer datatype,
                       integer n_repeats, double* time_min_);
void invoke_geqp3(integer datatype, integer* m, integer* n, void* a, integer* lda, integer *jpvt,
                  void* tau, void* work, integer* lwork, void* rwork, integer* info);

void fla_test_geqp3(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "QR factorization with column pivoting";
    char* front_str = "GEQP3";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_geqp3_experiment);

}

void fla_test_geqp3_experiment(test_params_t *params,
                               integer datatype,
                               integer p_cur,
                               integer q_cur,
                               integer pci,
                               integer n_repeats,
                               double* perf,
                               double* t,
                               double* residual)
{
    integer m, n, cs_A;
    void *A = NULL, *A_test = NULL, *T = NULL;
    integer *jpvt;
    double time_min = 1e9;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;

    /* Get input matrix dimensions */
    m = p_cur;
    n = q_cur;
    cs_A = m;

    /* Create input matrix parameters */
    create_matrix(datatype, &A, m, n);
    create_vector(datatype, &T, min(m,n));

    /* Initialize input matrix A with random numbers */
    rand_matrix(datatype, A, m, n, cs_A);

    /* Make a copy of input matrix A,required for validation. */
    create_matrix(datatype, &A_test, m, n);
    copy_matrix(datatype, "full", m, n, A, cs_A, A_test, cs_A);

    /* Create pivot array */
    create_vector(INTEGER, (void **) &jpvt, n);

    prepare_geqp3_run(m, n, A_test, jpvt, T, datatype, n_repeats, &time_min);

    /* execution time */
    *t = time_min;

    /* performance computation
     * 2mn^2 - (2/3)n^3 flops
     */
    if(m >= n)
        *perf = (double)((2.0 * m * n * n) - (( 2.0 / 3.0 ) * n * n * n )) / time_min / FLOPS_PER_UNIT_PERF;
    else
        *perf = (double)((2.0 * n * m * m) - (( 2.0 / 3.0 ) * m * m * m )) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    validate_geqp3(m, n, A, A_test, jpvt, T, datatype, residual);


    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(T);
    free_vector(jpvt);
}


void prepare_geqp3_run(integer m_A, integer n_A,
    void *A,
    integer *jpvt,
    void *T,
    integer datatype,
    integer n_repeats,
    double* time_min_)
{
    integer cs_A, min_A, i;
    void *A_save = NULL, *T_test = NULL, *work = NULL;
    void *rwork;
    integer lwork = -1, info = 0;
    double time_min = 1e9, exe_time;

    cs_A = m_A;
    min_A = min(m_A, n_A);

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion. */
    create_matrix(datatype, &A_save, m_A, n_A);
    copy_matrix(datatype, "full", m_A, n_A, A, cs_A, A_save, cs_A);

    /* Make a workspace query the first time. This will provide us with
       and ideal workspace size based on internal block size. */
    create_vector(datatype, &work, 1);

    /* rwork for complex types */
    if (datatype >= COMPLEX)
        create_realtype_vector(datatype, &rwork, 2 * n_A);

    /* call to  geqp3 API */
    invoke_geqp3(datatype, &m_A, &n_A, NULL, &cs_A, NULL, NULL, work, &lwork, rwork, &info);

    /* Get work size */
    lwork = get_work_value( datatype, work );

    /* Output buffers will be freshly allocated for each iterations, free up 
       the current output buffers. */ 
    free_vector(work);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration */
        copy_matrix(datatype, "full", m_A, n_A, A_save, cs_A, A, cs_A);

        /* Reset pivot buffer */
        reset_vector(INTEGER, jpvt, n_A, 1);

        /* T_test vector will hold the scalar factors of the elementary reflectors. */
        create_vector(datatype, &T_test, min_A);

        /* Create work buffer */
        create_matrix(datatype, &work, lwork, 1);

        exe_time = fla_test_clock();

        /* Call to  gerqf API */
        invoke_geqp3(datatype, &m_A, &n_A, A, &cs_A, jpvt, T_test, work, &lwork, rwork, &info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

        /* Make a copy of the output buffers, for validation. */
        copy_vector(datatype, min_A, T_test, 1, T, 1);

        /* Free up the output buffers */
        free_matrix(work);
        free_vector(T_test);
    }

    *time_min_ = time_min;

    free_matrix(A_save);
    if (datatype >= COMPLEX)
        free_vector(rwork);
}


void invoke_geqp3(integer datatype, integer* m, integer* n, void* a, integer* lda, integer *jpvt,
                  void* tau, void* work, integer* lwork, void *rwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            dgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            cgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
            break;
        }
    }
}
