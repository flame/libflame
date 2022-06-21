/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_getri_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_getri_run(integer m_A, integer n_A, void *A, integer* ipiv, integer datatype, integer n_repeats, double* time_min_);
void invoke_getri(integer datatype, integer *n, void *a, integer *lda, integer *ipiv, void *work, integer *lwork, integer *info);

void fla_test_getri(test_params_t *params)
{
    char* op_str = "Inverse through LU factorization";
    char* front_str = "GETRI";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_getri_experiment);

}


void fla_test_getri_experiment(test_params_t *params,
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
    void *A, *A_test;
    double time_min = 1e9;


    /* Determine the dimensions*/
    m = p_cur;
    cs_A = m;
    /* Create the matrices for the current operation*/
    create_matrix(datatype, &A, m, m);
    create_vector(INTEGER, &IPIV, m);
    /* Initialize the test matrices*/
    rand_matrix(datatype, A, m, m, cs_A);

    /* Save the original matrix*/
    create_matrix(datatype, &A_test, m, m);
    copy_matrix(datatype, "full", m, m, A, cs_A, A_test, cs_A);

    /* call to API */
    prepare_getri_run(m, m, A_test, IPIV, datatype, n_repeats, &time_min);

    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2mn^2 - (2/3)n^3 flops */
    *perf = (double)((2.0 * m * m * m) - ((2.0 / 3.0) * m * m * m)) / time_min / FLOPS_PER_UNIT_PERF;
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    validate_getri(m, m, A, A_test, IPIV, datatype, residual);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(IPIV);
}


void prepare_getri_run(integer m_A,
    integer n_A,
    void* A,
    integer* IPIV,
    integer datatype,
    integer n_repeats,
    double* time_min_)
{
    integer cs_A, lwork;
    integer i;
    void *A_save, *work;
    integer info = 0;
    double time_min = 1e9, exe_time;
    lwork = i_n_one;
    /* Get column stride */
    cs_A = m_A;
    /* Save the original matrix */
    create_matrix(datatype, &A_save, m_A, n_A);
    copy_matrix(datatype, "full", m_A, n_A, A, cs_A, A_save, cs_A);
    create_vector(datatype, &work, 1);

    // call to  getri API
    invoke_getri(datatype, &n_A, NULL, &n_A, NULL, work, &lwork, &info);

    // Get work size
    lwork = get_work_value(datatype, work);
    free_vector(work);

    for (i = 0; i < n_repeats; ++i)
    {

        /* Copy original input data */
        copy_matrix(datatype, "full", m_A, n_A, A, cs_A, A_save, cs_A);
        // create work buffer
        create_matrix(datatype, &work, lwork, 1);
        exe_time = fla_test_clock();

        /*  call to API getrf to get AFACT */
        invoke_getrf(datatype, &m_A, &n_A, A_save, &cs_A, IPIV, &info);
        /*  call  getri API with AFACT to get A INV */
        invoke_getri(datatype, &n_A, A_save, &n_A, IPIV, work, &lwork, &info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);
        // Free up the output buffers
        free_matrix(work);
    }

    *time_min_ = time_min;
    /*  Save the final result to A matrix*/
    copy_matrix(datatype, "full", m_A, n_A, A_save, cs_A, A, cs_A);
    free_matrix(A_save);

}


/*
 *  GETRI_API calls LAPACK interface of
 *  Singular value decomposition - gesvd
 *  */
void invoke_getri(integer datatype, integer *n, void *a, integer *lda, integer *ipiv, void* work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sgetri_(n, a, lda, ipiv, work, lwork, info);
            break;
        }
        
        case DOUBLE:
        {
            dgetri_(n, a, lda, ipiv, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            cgetri_(n, a, lda, ipiv, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zgetri_(n, a, lda, ipiv, work, lwork, info);
            break;
        }
    }
}

