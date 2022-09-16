/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_stevd_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
integer n_repeats, double* perf, double* t, double* residual);
void prepare_stevd_run(char* jobz, integer n, void* Z, void* D, void* E, integer datatype, integer n_repeats, double* time_min_);
void invoke_stevd(integer datatype, char* jobz, integer* n, void* z, integer* ldz, void* d, void* e, void* work, integer* lwork, void* iwork, integer* liwork, integer* info);

void fla_test_stevd(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Eigen Decomposition of symmetrix tridiagonal matrix";
    char* front_str = "STEVD";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_stevd_experiment);
}

void fla_test_stevd_experiment(test_params_t *params,
                               integer  datatype,
                               integer  p_cur,
                               integer  q_cur,
                               integer pci,
                               integer n_repeats,
                               double* perf,
                               double *time_min,
                               double* residual)
{
    integer n, ldz;
    char jobz;
    void *Z = NULL, *Z_test = NULL;
    void *D = NULL, *D_test = NULL;
    void *E = NULL, *E_test = NULL;

    /* Get input matrix dimensions.*/
    jobz = params->eig_sym_paramslist[pci].jobz;
    *residual = params->eig_sym_paramslist[pci].threshold_value;

    if(datatype == FLOAT || datatype == DOUBLE)
    {
        n = p_cur;
        ldz = max(1,n);

        /* Create input matrix parameters */
        create_matrix(datatype, &Z, n, n);
        reset_matrix(datatype, n, n, Z, ldz);
        create_vector(datatype, &D, n);
        create_vector(datatype, &E, n-1);

        /* input matrix Z with random symmetric numbers and D,E matrix with diagonal and subdiagonal values */
        rand_sym_tridiag_matrix(datatype, Z, n, n, ldz);
        get_diagonal(datatype, Z, n, n, ldz, D);
        get_subdiagonal(datatype, Z, n, n, ldz, E);

        /* Make a copy of input matrix A. This is required to validate the API functionality.*/
        create_matrix(datatype, &Z_test, n, n);
        reset_matrix(datatype, n, n, Z_test, ldz);
        create_vector(datatype, &D_test, n);
        create_vector(datatype, &E_test, n-1);
        copy_vector(datatype, n, D, 1, D_test, 1);
        copy_vector(datatype, n-1, E, 1, E_test, 1);

        prepare_stevd_run(&jobz, n, Z_test, D_test, E_test, datatype, n_repeats, time_min);

        /* performance computation
           6 * n^3 + n^2 flops for eigen vectors
           6 * n^2 flops for eigen values */
        if( jobz == 'V')
            *perf = (double)((6.0 * n * n * n) + (n * n)) / *time_min / FLOPS_PER_UNIT_PERF;
        else
            *perf = (double)(6.0 * n * n) / *time_min / FLOPS_PER_UNIT_PERF;

        /* output validation */
        validate_syevd(&jobz, n, Z, Z_test, D_test, datatype, residual);

        /* Free up the buffers */
        free_matrix(Z);
        free_vector(D);
        free_vector(E);
        free_matrix(Z_test);
        free_vector(D_test);
        free_vector(E_test);
    }
}

void prepare_stevd_run(char *jobz,
                       integer n,
                       void *Z,
                       void *D,
                       void *E,
                       integer datatype,
                       integer n_repeats,
                       double* time_min_)
{
    integer ldz;
    void *Z_save, *D_save, *E_save, *work, *iwork;
    integer lwork, liwork;
    integer i;
    integer info = 0;
    double time_min = 1e9, exe_time;

    ldz = max(1,n);

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &Z_save, n, n);
    create_vector(datatype, &D_save, n);
    create_vector(datatype, &E_save, n-1);

    copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);
    copy_vector(datatype, n, D, 1, D_save, 1);
    copy_vector(datatype, n-1, E, 1, E_save, 1);

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size.*/
    lwork = -1;
    liwork = -1;

    create_vector(datatype, &work, 1);
    create_vector(INTEGER, &iwork, 1);

    /* call to  gesdd API */
    invoke_stevd(datatype, jobz, &n, NULL, &ldz, NULL, NULL, work, &lwork, iwork, &liwork, &info);

    /* Get work size */
    lwork = get_work_value( datatype, work );
    liwork = get_work_value( INTEGER, iwork );

    /* Output buffers will be freshly allocated for each iterations, free up
       the current output buffers.*/
    free_vector(work);
    free_vector(iwork);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);
        copy_vector(datatype, n, D_save, 1, D, 1);
        copy_vector(datatype, n-1, E_save, 1, E, 1);

        create_vector(datatype, &work, lwork);
        create_vector(INTEGER, &iwork, liwork);

        exe_time = fla_test_clock();
        /* call to API */
        invoke_stevd(datatype, jobz, &n, Z, &ldz, D, E, work, &lwork, iwork, &liwork, &info);
        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(iwork);

    }

    *time_min_ = time_min;

    free(Z_save);
    free(D_save);
    free(E_save);
}

void invoke_stevd(integer datatype, char* jobz, integer* n, void* z, integer* ldz, void* d, void* e, void* work, integer* lwork, void* iwork, integer* liwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sstevd_(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
            break;
        }
        case DOUBLE:
        {
            dstevd_(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
            break;
        }
    }
}
