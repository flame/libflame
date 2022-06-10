/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1

/* Local prototypes.*/
void fla_test_syevd_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
integer n_repeats, double* perf, double* t, double* residual);
void prepare_syevd_run(char* jobz, char* uplo, integer n, void* A, void* w, integer datatype, integer n_repeats, double* time_min_);
void invoke_syevd(integer datatype, char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, void* rwork, integer* lrwork, void* iwork, integer* liwork, integer* info);

void fla_test_syevd(test_params_t *params)
{
    char* op_str = "Eigen Decomposition";
    char* front_str = "SYEVD";
    char* lapack_str = "LAPACK";
    char* pc_str[NUM_PARAM_COMBOS] = {""};

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS, params, EIG_SYM, fla_test_syevd_experiment);
}

void fla_test_syevd_experiment(test_params_t *params,
                               integer  datatype,
                               integer  p_cur,
                               integer  q_cur,
                               integer pci,
                               integer n_repeats,
                               double* perf,
                               double *time_min,
                               double* residual)
{
    integer n, lda;
    char jobz, uplo;
    void *A = NULL, *w = NULL, *A_test = NULL;

    /* Get input matrix dimensions.*/
    jobz = params->eig_sym_paramslist[pci].jobz;
    uplo = params->eig_sym_paramslist[pci].uplo;
    *residual = params->eig_sym_paramslist[pci].threshold_value;

    n = p_cur;
    lda = max(1,n);

    /* Create input matrix parameters */
    create_matrix(datatype, &A, n, n);
    create_realtype_vector(datatype, &w, n);

    /* input matrix A with random symmetric numbers or complex hermitian matrix */
    if(datatype == FLOAT || datatype == DOUBLE)
        rand_sym_matrix(datatype, A, n, n, lda);
    else
        rand_hermitian_matrix(datatype, n, &A, lda);

    /* Make a copy of input matrix A. This is required to validate the API functionality.*/
    create_matrix(datatype, &A_test, n, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    prepare_syevd_run(&jobz, &uplo, n, A_test, w, datatype, n_repeats, time_min);

    /* performance computation
       (8/3)n^3 flops for eigen vectors
       (4/3)n^3 flops for eigen values */
    if( jobz == 'V')
        *perf = (double)((8.0 / 3.0) * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    else
        *perf = (double)((4.0 / 3.0) * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if( jobz == 'V')
        validate_syevd(&jobz, &uplo, n, A, A_test, w, datatype, residual);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(w);
}

void prepare_syevd_run(char *jobz,
		       char *uplo,
                       integer n,
                       void *A,
                       void *w,
                       integer datatype,
                       integer n_repeats,
                       double* time_min_)
{
    integer lda;
    void *A_save, *w_test, *work, *iwork, *rwork=NULL;
    integer lwork, liwork, lrwork;
    integer i;
    integer info = 0;
    double time_min = 1e9, exe_time;

    lda = max(1,n);

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, n, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size.*/
    lwork = -1;
    liwork = -1;

    create_vector(datatype, &work, 1);
    create_vector(INTEGER, &iwork, 1);

    if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX )
    {
	lrwork = -1;
        create_realtype_vector(datatype, &rwork, 1);
    }
    else
    {
        rwork = NULL;
    }

    /* call to  gesdd API */
    invoke_syevd(datatype, jobz, uplo, &n, NULL, &lda, NULL, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

    /* Get work size */
    lwork = get_work_value( datatype, work );
    liwork = get_work_value( INTEGER, iwork );
    lrwork = get_work_value( datatype, rwork );

    /* Output buffers will be freshly allocated for each iterations, free up
       the current output buffers.*/
    free_vector(work);
    free_vector(iwork);
    if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        free_vector(rwork);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);

        create_realtype_vector(datatype, &w_test, n);
        create_vector(datatype, &work, lwork);
        create_vector(INTEGER, &iwork, liwork);

	if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX )
            create_realtype_vector(datatype, &rwork, lrwork);
        else
            rwork = NULL;

        exe_time = fla_test_clock();

        /* call to API */
        invoke_syevd(datatype, jobz, uplo, &n, A, &lda, w_test, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

        /* Make a copy of the output buffers. This is required to validate the API functionality.*/
        copy_realtype_vector(datatype, n, w_test, 1, w, 1);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(iwork);

	if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            free_vector(rwork);

	free_vector(w_test);
    }

    *time_min_ = time_min;

    free(A_save);
}

void invoke_syevd(integer datatype, char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, void* rwork, integer* lrwork, void* iwork, integer* liwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            ssyevd_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
            break;
        }
        case DOUBLE:
        {
            dsyevd_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
            break;
        }
        case COMPLEX:
        {
            cheevd_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zheevd_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
            break;
        }
    }
}
