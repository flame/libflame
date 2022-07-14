/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_orgqr_experiment(test_params_t *params, integer datatype, integer  p_cur, integer  q_cur, integer  pci, integer  n_repeats, double* perf, double* t,double* residual);
void prepare_orgqr_run(integer m, integer n, void *A, void *T, void* work, integer *lwork, integer datatype, integer n_repeats, double* time_min_);
void invoke_orgqr(integer datatype, integer* m, integer* n, integer *min_A, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);

void fla_test_orgqr(test_params_t *params)
{
    char* op_str = "QR factorization";
    char* front_str = "ORGQR";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_orgqr_experiment);
}

void fla_test_orgqr_experiment(test_params_t *params,
    integer datatype,
    integer p_cur,
    integer q_cur,
    integer pci,
    integer n_repeats,
    double* perf,
    double* time_min,
    double* residual)
{
    integer m, n, cs_A;
    void *A = NULL, *A_test = NULL, *T_test = NULL;
    void *work = NULL, *work_test = NULL;
    void *Q = NULL, *R = NULL, *R_test = NULL;
    integer lwork = -1, info = 0;

    /* Get input matrix dimensions.*/
    m = p_cur;
    n = q_cur;
    cs_A = m;
    *time_min = 0.;
    *perf = 0.;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;

    if(m >= n)
    {
        /* Create input matrix parameters */
        create_matrix(datatype, &A, m, n);

        /* create tau vector */
        create_vector(datatype, &T_test, min(m,n));

        /* Initialize input matrix A with random numbers */
        rand_matrix(datatype, A, m, n, cs_A);

        /* Make a copy of input matrix A. This is required to validate the API functionality.*/
        create_matrix(datatype, &A_test, m, n);
        copy_matrix(datatype, "full", m, n, A, cs_A, A_test, cs_A);

        /* create Q matrix to check orthogonality */
        create_matrix(datatype, &Q, m, n);
        reset_matrix(datatype, m, n, Q, m);
 
        /* Make a workspace query the first time. This will provide us with
           and ideal workspace size based on internal block size.*/
        create_vector(datatype, &work, 1);

        /* call to  geqrf API */
        invoke_geqrf(datatype, &m, &n, NULL, &cs_A, NULL, work, &lwork, &info);

        /* Get work size */
        lwork = get_work_value( datatype, work );

        /* Output buffers will be freshly allocated for each iterations, free up
           the current output buffers.*/
        free_vector(work);

        /* create work buffer */
        create_matrix(datatype, &work, lwork, 1);
        create_vector(datatype, &work_test, lwork);

        /* QR Factorisation on matrix A to generate Q and R */
        invoke_geqrf(datatype, &m, &n, A_test, &cs_A, T_test, work, &lwork, &info);

        if(m == n)
        {
            create_matrix(datatype, &R, m, n);
            reset_matrix(datatype, m, n, R, cs_A);
            copy_matrix(datatype, "Upper", m, n, A_test, m, R, m);
        }
        else
        {
            create_matrix(datatype, &R, n, n);
            create_matrix(datatype, &R_test, n, n);
            reset_matrix(datatype, n, n, R, n);
            reset_matrix(datatype, n, n, R_test, n);
            copy_submatrix(datatype, A_test, m, n, R_test, n, n, 0, 0);
            copy_matrix(datatype, "Upper", n, n, R_test, n, R, n);
        } 
        copy_matrix(datatype, "full", m, n, A_test, m, Q, m);

        /*invoke orgqr API */
        prepare_orgqr_run(m, n, Q, T_test, work_test, &lwork, datatype, n_repeats, time_min);

        /* performance computation
           (2/3)*n2*(3m - n) */
        *perf = (double)((2.0 * m * n * n) - (( 2.0 / 3.0 ) * n * n * n )) / *time_min / FLOPS_PER_UNIT_PERF;
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            *perf *= 4.0;

        /* output validation */
        validate_orgqr(m, n, A, Q, R, work_test, datatype, residual);

        /* Free up the buffers */
        free_matrix(A);
        free_matrix(A_test);
        free_matrix(work);
        free_vector(work_test);
        free_vector(T_test);
        free_matrix(Q);
        free_matrix(R);
        if(m > n)
            free_matrix(R_test);
    }
}

void prepare_orgqr_run(integer m, integer n,
    void* A,
    void* T,
    void* work,
    integer *lwork,
    integer datatype,
    integer n_repeats,
    double* time_min_)
{
    integer cs_A, i;
    void *A_save = NULL;
    integer info = 0;
    double time_min = 1e9, exe_time;

    cs_A = m;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, m, n);
    copy_matrix(datatype, "full", m, n, A, cs_A, A_save, cs_A);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m, n, A_save, cs_A, A, cs_A);

        exe_time = fla_test_clock();

        /* Call to  orgqr API */
        invoke_orgqr(datatype, &m, &n, &n, A, &cs_A, T, work, lwork, &info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

    }

    *time_min_ = time_min;

    free_matrix(A_save);
}


void invoke_orgqr(integer datatype, integer* m, integer* n, integer *min_A, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sorgqr_(m, n, n, a, m, tau, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            dorgqr_(m, n, n, a, m, tau, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            cungqr_(m, n, n, a, m, tau, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zungqr_(m, n, n, a, m, tau, work, lwork, info);
            break;
        }
    }
}
