/*
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_getrs_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, integer einfo, double* perf, double* t, double* residual);
void prepare_getrs_run(char *trans, integer m_A, integer n_A, void *A, integer lda, void *B, integer ldb, integer* ipiv, integer datatype, integer n_repeats, double* time_min_, integer *info);
void invoke_getrs(integer datatype, char *trans, integer *nrhs, integer *n, void *a, integer *lda, integer *ipiv, void *b, integer *ldb, integer *info);

static FILE* g_ext_fptr = NULL;

void fla_test_getrs(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "LU factorization";
    char* front_str = "GETRS";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_getrs_experiment);
        tests_not_run = 0;
    }
    if(argc == 10)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[9]);
    }
    if(argc >= 9 && argc <= 10)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].transr = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldb = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        
        n_repeats = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->lin_solver_paramslist[0].solver_threshold = CLI_NORM_THRESH;

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
                fla_test_getrs_experiment(params, datatype,
                                          N, N,
                                          0,
                                          n_repeats, einfo,
                                          &perf, &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str,
                                      stype,
                                      SQUARE_INPUT,
                                      N, N,
                                      residual, params->lin_solver_paramslist[0].solver_threshold,
                                      time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for getrs\n");
        printf("./<EXE> getrs <precisions - sdcz> <TRANS> <N> <NRHS> <LDA> <LDB> <repeats>\n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
    if (g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
    }
    return;

}


void fla_test_getrs_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    integer einfo,
    double* perf,
    double* t,
    double* residual)
{
    integer n, lda, ldb, NRHS;
    integer info = 0, vinfo = 0;
    void* IPIV;
    void *A, *A_test, *B, *B_save, *X;
    double time_min = 1e9;
    char TRANS = params->lin_solver_paramslist[pci].transr;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    NRHS = params->lin_solver_paramslist[pci].nrhs;

    /* Determine the dimensions*/
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    ldb = params->lin_solver_paramslist[pci].ldb;

    if(lda < n || ldb < n)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Create the matrices for the current operation*/
    create_matrix(datatype, &A, lda, n);
    create_vector(INTEGER, &IPIV, n);
    create_matrix(datatype, &B, ldb, NRHS);
    create_matrix(datatype, &B_save, ldb, NRHS);
    create_matrix(datatype, &X, n, NRHS);
    create_matrix(datatype, &A_test, lda, n);
    /* Initialize the test matrices*/
    if (g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
        init_matrix_from_file(datatype, B, n, NRHS, ldb, g_ext_fptr);
    }
    else
    {
        /* Initialize input matrix with random numbers */
        rand_matrix(datatype, A, n, n, lda);
        rand_matrix(datatype, B, n, NRHS, ldb);
    }

    /* Save the original matrix*/

    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);
    copy_matrix(datatype, "full", n, NRHS, B, ldb, B_save, ldb);

    /*  call to API getrf to get AFACT */
    invoke_getrf(datatype, &n, &n, A_test, &lda, IPIV, &info);

    /* call to API */
    prepare_getrs_run(&TRANS, n, NRHS, A_test, lda, B, ldb, IPIV, datatype, n_repeats, &time_min, &info);
    copy_matrix(datatype, "full", n, NRHS, B, ldb, X, n);
    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2mn^2 - (2/3)n^3 flops */
    *perf = (double)((2.0 * n * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if(info == 0)
        validate_getrs(&TRANS, n, NRHS, A, lda, B_save, ldb, X, datatype, residual, &vinfo);
    
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
        
    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(IPIV);
    free_matrix(B);
    free_matrix(X);
    free_matrix(B_save);
}


void prepare_getrs_run(char *TRANS,
    integer n_A,
    integer nrhs,
    void* A,
    integer lda,
    void* B,
    integer ldb,
    integer* IPIV,
    integer datatype,
    integer n_repeats,
    double* time_min_,
    integer* info)
{
    integer i;
    void *A_save, *B_test;
    double time_min = 1e9, exe_time;

    /* Save the original matrix */
    create_matrix(datatype, &A_save, lda, n_A);
    copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);
    create_matrix(datatype, &B_test, ldb, nrhs);


    *info = 0;
    for (i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Copy original input data */
        copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);
        copy_matrix(datatype, "full", n_A, nrhs, B, ldb, B_test, ldb);

        exe_time = fla_test_clock();

        /*  call  getrs API with AFACT */
        invoke_getrs(datatype, TRANS, &n_A, &nrhs, A_save, &lda, IPIV, B_test, &ldb, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

    }

    *time_min_ = time_min;
    /*  Save the final result to B matrix*/
    copy_matrix(datatype, "full", n_A, nrhs, B_test, ldb, B, ldb);

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
            fla_lapack_sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }
    }
}
