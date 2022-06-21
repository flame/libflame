/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes */
void fla_test_gesvd_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci, integer n_repeats, double* perf, double* t, double* residual);
void prepare_gesvd_run(char *jobu, char *jobvt, integer m_A, integer n_A, void *A, void *s, void *U, void *V, integer datatype, integer n_repeats, double* time_min_);
void invoke_gesvd(integer datatype, char* jobu, char *jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* iwork, integer* info);

void fla_test_gesvd(test_params_t *params)
{
    char* op_str = "Singular value decomposition";
    char* front_str = "GESVD";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, RECT_INPUT, params, SVD, fla_test_gesvd_experiment);
}


void fla_test_gesvd_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    double* perf,
    double *time_min,
    double* residual)
{
    integer m, n, cs_A;
    char jobu, jobvt;
    void *A = NULL, *U = NULL, *V = NULL, *s = NULL, *A_test = NULL;
    *residual = params->svd_paramslist[pci].svd_threshold;

    /* Get input matrix dimensions. */
    jobu = params->svd_paramslist[pci].jobu_gesvd;
    jobvt = params->svd_paramslist[pci].jobvt_gesvd;

    if(jobu == 'S' || jobvt == 'S' || (jobu == 'O' && jobvt != 'O') || (jobu != 'O' && jobvt == 'O'))
    {
           m = p_cur;
           n = p_cur;
    }
    else
    {
           m = p_cur;
           n = q_cur;
    }

    cs_A = m;

    /* Create input matrix parameters. */
    create_matrix(datatype, &A, m, n);
    create_matrix(datatype, &U, m, m);
    create_matrix(datatype, &V, n, n);
    create_realtype_vector(datatype, &s, min(m, n));

    /* Initialize input matrix A with random numbers. */
    rand_matrix(datatype, A, m, n, cs_A);

    /* Make a copy of input matrix A. This is required to validate the API functionality. */
    create_matrix(datatype, &A_test, m, n);
    copy_matrix(datatype, "full", m, n, A, cs_A, A_test, cs_A);

    prepare_gesvd_run(&jobu, &jobvt, m, n, A_test, s, U, V, datatype, n_repeats, time_min);

    /* Performance Computation
     * Singular values only, 4mn^2 - 4n^3/3 flops
     * Singular values and some singular vectors U (m x n) and V (n x n), 14mn^2 + 8n^3 flops 
     * Link : http://icl.cs.utk.edu/magma/forum/viewtopic.php?f=2&t=921 */
    if(jobu == 'N' || jobvt == 'N')
    {
           if(m >= n)
              *perf = (double)((4.0 * m * n * n) - (( 4.0 * n * n * n ) / 3.0 )) / *time_min / FLOPS_PER_UNIT_PERF;
           else
              *perf = (double)((4.0 * n * m * m) - ((4.0 * m * m * m ) / 3.0)) / *time_min / FLOPS_PER_UNIT_PERF;
    }
    else
    {
       if(m >= n)
              *perf = (double)((14.0 * m * n * n) + ( 8.0  * n * n * n )) / *time_min / FLOPS_PER_UNIT_PERF;
           else
              *perf = (double)((14.0 * n * m * m) + (8.0 * m * m * m )) / *time_min / FLOPS_PER_UNIT_PERF;
    }
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
              *perf *= 4.0;

    /* output validation */
    if(jobu == 'A' && jobvt == 'A')
           validate_gesvd(&jobu, &jobvt, m, n, A, A_test, s, U, V, datatype, residual); 

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(U);
    free_matrix(V);
    free_vector(s);
}



void prepare_gesvd_run(char *jobu, char *jobvt,
    integer m_A, integer n_A,
    void *A,
    void *s,
    void *U,
    void *V,
    integer datatype,
    integer n_repeats,
    double* time_min_)
{
    integer cs_A, cs_U, cs_V;
    integer min_m_n, max_m_n;
    void *A_save, *s_test;
    void *work, *iwork, *rwork;
    void *U_test, *V_test;
    integer lwork, liwork, lrwork;
    integer i;
    integer info = 0;
    double time_min = 1e9, exe_time;


    cs_A = m_A;
    cs_U = m_A;
    cs_V = n_A;
    min_m_n = min(m_A, n_A);
    max_m_n = max(m_A, n_A);

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, m_A, n_A);
    copy_matrix(datatype, "full", m_A, n_A, A, m_A, A_save, m_A);

    /* Get rwork and iwork array size since it is not depedent on internal blocks*/
    lrwork = max( (5 * min_m_n * min_m_n + 5 * min_m_n) , ( 2 * max_m_n * min_m_n + 2 * min_m_n * min_m_n + min_m_n));
    liwork = 8 * min_m_n;

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size. */
    lwork = -1;
    create_vector(datatype, &work, 1);

    /* call gesvd API */
    invoke_gesvd(datatype, jobu, jobvt, &m_A, &n_A, NULL, &cs_A, NULL, NULL, &cs_U, NULL, &cs_V, work, &lwork, NULL, NULL, &info);

    /* Get the work size */
    lwork = get_work_value( datatype, work );

    /* Output buffers will be freshly allocated for each iterations, free up 
       the current output buffers.*/ 
    free_vector(work);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m_A, n_A, A_save, m_A, A, m_A);

        create_matrix(datatype, &U_test, m_A, m_A);
        create_matrix(datatype, &V_test, n_A, n_A);
        create_realtype_vector(datatype, &s_test, min_m_n);
        create_vector(datatype, &work, lwork);
        create_vector(INTEGER, &iwork, liwork);

        if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            create_realtype_vector(datatype, &rwork, lrwork);
        else
            rwork = NULL;


        exe_time = fla_test_clock();

        /* call to API */
        invoke_gesvd(datatype, jobu, jobvt, &m_A, &n_A, A, &cs_A, s_test, U_test, &cs_U, V_test, &cs_V, work, &lwork, rwork, iwork, &info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

        /* Make a copy of the output buffers. This is required to validate the API functionality. */
        copy_matrix(datatype, "full", m_A, m_A, U_test, m_A, U, m_A);
        copy_matrix(datatype, "full", n_A, n_A, V_test, n_A, V, n_A);
        copy_realtype_vector(datatype, min_m_n, s_test, 1, s, 1);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(iwork);
        if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            free_vector(rwork);
        free_matrix(U_test);
        free_matrix(V_test);
        free_vector(s_test);
    }

    *time_min_ = time_min;

    free(A_save);
}


void invoke_gesvd(integer datatype, char* jobu, char* jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* iwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
            break;
        }
        
        case DOUBLE:
        {
            dgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            cgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
            break;
        }
    }
}
