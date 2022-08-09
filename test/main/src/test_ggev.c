/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"


/* Local prototypes */
void fla_test_ggev_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci, integer n_repeats, double* perf, double* t, double* residual);
void prepare_ggev_run(char *jobvl, char *jobvr, integer n, void *a, void *b, void* alpha, void * alphar, void * alphai, void *beta,	void *vl, integer ldvl, 
                      void *vr, integer ldvr,	integer datatype, integer n_repeats, double* time_min_);
void invoke_ggev(integer datatype, char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* b, integer* ldb, integer* alpha, integer* alphar,
    integer* alphai, integer* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, void* rwork, integer* info);


void fla_test_ggev(test_params_t *params)
{
    char* op_str = "Computing Eigen value and Eigen vectors";
    char* front_str = "GGEV";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_NSYM, fla_test_ggev_experiment);
}


void fla_test_ggev_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double   *perf,
    double   *t,
    double   *residual)
{
    integer m, n, cs_A,ldvl, ldvr;
    void *A = NULL, *B = NULL, *VL = NULL, *VR = NULL;
    void * alpha = NULL, *alphar=NULL, *alphai=NULL, *beta, *A_test , *B_test;
    double time_min = 1e9;
    *residual = params->eig_non_sym_paramslist[pci].GenNonSymEigProblem_threshold;
    //char JOBVL = params->eig_non_sym_paramslist[pci].jobvsl;
    //char JOBVR = params->eig_non_sym_paramslist[pci].jobvsr;
    char JOBVL = 'V';
    char JOBVR = 'V';

    /* Get input matrix dimensions */
    m = p_cur;
    n = q_cur;
    cs_A = n;
    ldvl = n;
    ldvr = n;

    /* Create input matrix parameters */
    create_matrix(datatype, &A, n, n);
    create_matrix(datatype, &B, n, n);
    create_matrix(datatype, &VL, ldvl, n);
    create_matrix(datatype, &VR, ldvr, n);
    if (datatype == FLOAT || datatype == DOUBLE)
    {
        create_vector(datatype, &alphar, n);
        create_vector(datatype, &alphai, n);
    }
    else
    {
        create_vector(datatype, &alpha, n);
    }
    create_vector(datatype, &beta, n);

    /* Initialize input matrix A with random numbers */
    rand_matrix(datatype, A, n, n, cs_A);
    rand_matrix(datatype, B, n, n, cs_A);

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, &A_test, n, n);
    create_matrix(datatype, &B_test, n, n);
    copy_matrix(datatype, "full", n, n, A, cs_A, A_test, cs_A);
    copy_matrix(datatype, "full", n, n, B, cs_A, B_test, cs_A);

    prepare_ggev_run(&JOBVL,&JOBVR, n, A_test, B_test, alpha, alphar, alphai, beta,  VL, ldvl, VR, ldvr, datatype, n_repeats, &time_min);

    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2mn^2 - (2/3)n^3 flops */
    *perf = (double)((2.0 * m * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if (JOBVL == 'V' && JOBVR == 'V')
    {
        validate_ggev(&JOBVL, &JOBVR, n, A, B, alpha, alphar, alphai, beta, VL, ldvl, VR, ldvr, datatype, residual);
    }
    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(VL);
    free_matrix(VR);
    if (datatype == FLOAT || datatype == DOUBLE)
    {
        free_vector(alphar);
        free_vector(alphai);
    }
    else
    {
        free_vector(alpha);
    }

    free_vector(beta);

}

void prepare_ggev_run(char *jobvl, char *jobvr, integer n_A, void *A,
                        void *B, void* alpha, void * alphar, void * alphai, void* beta,
                        void *VL, integer ldvl, void *VR, integer ldvr,
                        integer datatype, integer n_repeats, double* time_min_)
{
    void *A_save = NULL, *B_save = NULL , *work = NULL, * rwork = NULL;
    integer i, lda, ldb;
    integer lwork, info = 0;
    double time_min = 1e9, exe_time;

    lda = n_A;
    ldb = n_A;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, n_A, n_A);
    copy_matrix(datatype, "full", n_A, n_A, A, n_A, A_save, n_A);
    create_matrix(datatype, &B_save, n_A, n_A);
    copy_matrix(datatype, "full", n_A, n_A, B, n_A, B_save, n_A);

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size. */
    lwork = -1;
    create_vector(datatype, &work, 8*n_A);

    /* call to  ggev API to get work query */
    invoke_ggev(datatype, jobvl, jobvr, &n_A, NULL, &lda, NULL, &ldb, NULL, NULL, NULL, NULL, NULL, &ldvl, NULL, &ldvr, work, &lwork, rwork, &info);

    /* Get work size */
    lwork = get_work_value( datatype, work);

    /* Output buffers will be freshly allocated for each iterations, free up 
       the current output buffers.*/ 
    free_vector(work);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers for each iteration */
        copy_matrix(datatype, "full", n_A, n_A, A_save, n_A, A, n_A);
        copy_matrix(datatype, "full", n_A, n_A, B_save, n_A, B, n_A);
        create_vector(datatype, &work, lwork);
        if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            create_realtype_vector(datatype, &rwork, 8*n_A);
        }
        exe_time = fla_test_clock();

        /* call to API */

        invoke_ggev(datatype, jobvl, jobvr, &n_A, A, &lda, B, &ldb, alpha, alphar, alphai, beta, VL, &ldvl, VR, &ldvr, work, &lwork, rwork , &info);
        
        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
        if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            free_vector(rwork);
        }
    }

    *time_min_ = time_min;
    copy_matrix(datatype, "full", n_A, n_A, A_save, n_A, A, n_A);
    copy_matrix(datatype, "full", n_A, n_A, B_save, n_A, B, n_A);

    free(A_save);
    free(B_save);


}


void invoke_ggev(integer datatype, char *jobvl, char *jobvr,integer *n, void *a, integer *lda, void *b, integer *ldb, integer* alpha, integer * alphar,
    integer* alphai, integer* beta, void *vl, integer *ldvl, void *vr, integer *ldvr, void* work, integer* lwork, void* rwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sggev_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
            break;
        }
        
        case DOUBLE:
        {
            dggev_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            cggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
            break;
        }
    }
}