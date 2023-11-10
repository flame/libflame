/*
    Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_syevx_experiment(test_params_t *params, integer datatype,
                              integer p_cur, integer  q_cur, integer pci,
                              integer n_repeats, integer einfo, double* perf,
                              double* t, double* residual);
void prepare_syevx_run(char* jobz, char* range, char* uplo, integer n, void* A,
                       integer lda, void *vl, void *vu, integer il,
                       integer iu, void *abstol, void* w, integer ldz,
                       integer datatype, integer n_repeats,
                       double* time_min_, integer* info);
void invoke_syevx(integer datatype, char* jobz, char* range, char* uplo,
                  integer* n, void* a, integer* lda, void* vl, void* vu,
                  integer* il, integer* iu, void* abstol, integer* m, void* w,
                  void* z, integer* ldz, void* work, integer* lwork,
                  void* rwork, void* iwork, void* ifail, integer* info);

void fla_test_syevx(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Eigen Values and Vectors in specified range";
    char* front_str = "SYEVX";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_syevx_experiment);
        tests_not_run = 0;
    }
    if (argc == 17)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[16]);
    }
    if (argc >= 16 && argc <= 17)
    {
        integer i, num_types,N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype,type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].jobz = argv[3][0];
        params->eig_sym_paramslist[0].range_x = argv[4][0];
        params->eig_sym_paramslist[0].uplo = argv[5][0];
        N = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].lda = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);

        params->eig_sym_paramslist[0].VL = atof(argv[8]);
        params->eig_sym_paramslist[0].VU = atof(argv[9]);

        if (params->eig_sym_paramslist[0].range_x == 'I')
        {
            /* 1 <= IL <= IU <= N, if N > 0;
               IL = 1 and IU = 0 if N = 0. */
            if (N == 0)
            {
                params->eig_sym_paramslist[0].IL = 1;
                params->eig_sym_paramslist[0].IU = 0;
                printf("\nIL = 1 and IU = 0 if N = 0\n");
            }
            else
            {
                params->eig_sym_paramslist[0].IL = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
                params->eig_sym_paramslist[0].IU = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
            }
        }

        params->eig_sym_paramslist[0].abstol = atof(argv[12]);

        params->eig_sym_paramslist[0].ldz = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);

        g_lwork = strtoimax(argv[14], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[15], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->eig_sym_paramslist[0].threshold_value = CLI_NORM_THRESH;

            for(i = 0; i < num_types; i++)
            {
                stype = argv[2][i];
                datatype = get_datatype(stype);

                /* Check for invalide dataype */
                if(datatype == INVALID_TYPE)
                {
                    invalid_dtype = 1;
                    continue;
                }

                /* Check for duplicate datatype presence */
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                /* Call the test code */
                fla_test_syevx_experiment(params, datatype,
                                          N, N,
                                          0,
                                          n_repeats, einfo,
                                          &perf, &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str,
                                      stype,
                                      SQUARE_INPUT,
                                      N, N,
                                      residual, params->eig_sym_paramslist[0].threshold_value,
                                      time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for syevx\n");
        printf("./<EXE> syevx <precisions - sdcz> <JOBZ> <RANGE> <UPLO> <N> <LDA> <VL> <VU> <IL> <IU> <ABSTOL> <LDZ> <LWORK> <repeats>\n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
    if (g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
        g_ext_fptr = NULL;
    }
    return;
}

void fla_test_syevx_experiment(test_params_t *params,
                               integer  datatype,
                               integer  p_cur,
                               integer  q_cur,
                               integer pci,
                               integer n_repeats,
                               integer einfo,
                               double* perf,
                               double *time_min,
                               double* residual)
{
    integer n, lda, ldz, il, iu, info = 0, vinfo = 0;
    char jobz, uplo, range;
    void *A = NULL, *w = NULL, *A_test = NULL;
    void *vl, *vu, *abstol;

    /* Get input matrix dimensions.*/
    jobz = params->eig_sym_paramslist[pci].jobz;
    uplo = params->eig_sym_paramslist[pci].uplo;
    range = params->eig_sym_paramslist[pci].range_x;
    *residual = params->eig_sym_paramslist[pci].threshold_value;

    n = p_cur;
    lda = params->eig_sym_paramslist[pci].lda;
    ldz = params->eig_sym_paramslist[pci].ldz;

    il = params->eig_sym_paramslist[pci].IL;
    iu = params->eig_sym_paramslist[pci].IU;

    create_realtype_vector(datatype, &vl, 1);
    create_realtype_vector(datatype, &vu, 1);
    create_realtype_vector(datatype, &abstol, 1);

    if (datatype == FLOAT || datatype == COMPLEX)
    {
        *(real*)vl = params->eig_sym_paramslist[pci].VL;
        *(real*)vu = params->eig_sym_paramslist[pci].VU;
        *(real*)abstol = params->eig_sym_paramslist[pci].abstol;

        /* When abstol value is set to -1, assign default value.
           NOTE: Eigenvalues will be computed most accurately
                 when ABSTOL is set to twice the underflow
                 threshold 2*SLAMCH('S') */
        if (*(real*)abstol == -1)
            *(real*)abstol = 2 * slamch_("S");
    }
    else
    {
        *(doublereal*)vl = params->eig_sym_paramslist[pci].VL;
        *(doublereal*)vu = params->eig_sym_paramslist[pci].VU;
        *(doublereal*)abstol = params->eig_sym_paramslist[pci].abstol;

        /* When abstol value is set to -1, assign default value.
           NOTE: Eigenvalues will be computed most accurately
                 when ABSTOL is set to twice the underflow
                 threshold 2*DLAMCH('S') */
        if (*(doublereal*)abstol == -1)
            *(doublereal*)abstol = 2 * dlamch_("S");
    }

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if (config_data)
    {
        if (lda == -1)
        {
            lda = fla_max(1,n);
        }
        /* LDZ >= 1;
           if JOBZ = 'V', LDZ >= max(1,N) */
        if (ldz == -1)
        {
            if (jobz == 'V')
            {
                ldz = fla_max(1,n);
            }
            else
            {
                ldz = 1;
            }
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, &A, lda, n);
    create_realtype_vector(datatype, &w, n);
    if (g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
    }
    else
    {
        /* input matrix A with random symmetric numbers
           or complex hermitian matrix */
        if (datatype == FLOAT || datatype == DOUBLE)
            rand_sym_matrix(datatype, A, n, n, lda);
        else
            rand_hermitian_matrix(datatype, n, &A, lda);
    }
    /* Make a copy of input matrix A.
       This is required to validate the API functionality.*/
    create_matrix(datatype, &A_test, lda, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    prepare_syevx_run(&jobz, &range, &uplo, n, A_test, lda, vl, vu, il, iu,
                      abstol, w, ldz, datatype, n_repeats, time_min,
                      &info);

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
    if (info == 0 && range == 'A')
        validate_syevd(&jobz, n, A, A_test, lda, w, datatype, residual, &vinfo);

    FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_vector(vl);
    free_vector(vu);
    free_vector(abstol);
    free_matrix(A);
    free_matrix(A_test);
    free_vector(w);
}

void prepare_syevx_run(char* jobz, char* range, char* uplo, integer n, void* A,
                       integer lda, void* vl, void* vu, integer il,
                       integer iu, void* abstol, void* w, integer ldz,
                       integer datatype, integer n_repeats,
                       double* time_min_, integer* info)
{
    void *A_save = NULL, *work = NULL, *rwork = NULL;
    void *w_test = NULL, *z__ = NULL;
    integer i, m, lwork;
    double time_min = 1e9, exe_time;
    void *iwork = NULL, *ifail = NULL;

    if(*range == 'I')
        m = iu - il + 1;
    else
        m = n;

    /* Make a copy of the input matrix A.
       Same input values will be passed in eaach itertaion.*/
    create_matrix(datatype, &A_save, lda, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);
    create_vector(INTEGER, &iwork, 5*n);
    create_vector(INTEGER, &ifail, n);

    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX )
        create_realtype_vector(datatype, &rwork, (7*n));
    else
        rwork = NULL;

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size.*/
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
        /* call to  syevx API */
        invoke_syevx(datatype, jobz, range, uplo, &n, NULL, &lda, vl, vu,
                     &il, &iu, abstol, &m, NULL, NULL, &ldz, work, &lwork,
                     rwork, iwork, ifail, info);
        /* Get work size */
        if(*info == 0)
        {
            lwork = get_work_value(datatype, work);
        }
        free_vector(work);
    }
    else
    {
        lwork = g_lwork;
    }

    *info = 0;
    for (i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);

        create_realtype_vector(datatype, &w_test, n);
        create_vector(datatype, &work, lwork);
        create_matrix(datatype, &z__, ldz, fla_max(1, m));

        exe_time = fla_test_clock();

        /* call to API */
        invoke_syevx(datatype, jobz, range, uplo, &n, A, &lda, vl, vu, &il,
                     &iu, abstol, &m, w_test, z__, &ldz, work, &lwork, rwork,
                     iwork, ifail, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Make a copy of the output buffers.
           This is required to validate the API functionality.*/
        copy_realtype_vector(datatype, n, w_test, 1, w, 1);

        /* If JOBZ = 'V', the first M columns of Z contain the
           orthonormal eigenvectors of the matrix A corresponding to
           the selected eigenvalues.
           Copy eigen vectors to A to validate API functionality */
        if(*jobz == 'V')
            copy_matrix(datatype, "full", m, m, z__, ldz, A, lda);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(w_test);
        free_matrix(z__);
    }

    *time_min_ = time_min;
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        free_vector(rwork);
    free_vector(iwork);
    free_vector(ifail);
    free_matrix(A_save);
}

void invoke_syevx(integer datatype, char* jobz, char* range, char* uplo,
                  integer* n, void* a, integer* lda, void* vl, void* vu,
                  integer* il, integer* iu, void* abstol, integer* m, void* w,
                  void* z, integer* ldz, void* work, integer* lwork,
                  void* rwork, void* iwork, void* ifail, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu,
                              abstol, m, w, z, ldz, work, lwork, iwork, ifail,
                              info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu,
                              abstol, m, w, z, ldz, work, lwork, iwork, ifail,
                              info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu,
                              abstol, m, w, z, ldz, work, lwork, rwork, iwork,
                              ifail, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu,
                              abstol, m, w, z, ldz, work, lwork, rwork, iwork,
                              ifail, info);
            break;
        }
    }
}
