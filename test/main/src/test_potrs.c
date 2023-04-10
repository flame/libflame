/*
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

// Local prototypes.
void fla_test_potrs_experiment(test_params_t *params, integer datatype, integer  p_cur, integer  q_cur, integer  pci, integer  n_repeats, integer einfo, double* perf, double* time_min,double* residual);
void prepare_potrs_run(char* uplo, integer m, integer nrhs, void *A, integer lda, integer datatype, void *b, integer ldb, integer n_repeats, double* time_min_, integer *info);
void invoke_potrs(char* uplo, integer datatype, integer* m, void* A, integer* lda, integer *nrhs, void* b, integer* ldb, integer* info);
static FILE* g_ext_fptr = NULL;

void fla_test_potrs(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Cholesky factorization";
    char* front_str = "POTRS";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_potrs_experiment);
        tests_not_run = 0;
    }
    if (argc == 10)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[9]);
    }
    if (argc >= 9 && argc <= 10)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].Uplo = argv[3][0];
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
                fla_test_potrs_experiment(params, datatype,
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
        printf("\nIllegal arguments for potrs\n");
        printf("./<EXE> potrs <precisions - sdcz> <UPLO> <N> <NRHS> <LDA> <LDB> <repeats>\n");
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

void fla_test_potrs_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    integer  einfo,
    double* perf,
    double* t,
    double* residual)
{
    integer n, info = 0, nrhs, lda, ldb, vinfo = 0;
    void *A = NULL, *A_test = NULL;
    void *B = NULL, *X = NULL;
    void *B_test = NULL;
    double time_min = 1e9;
    char uplo = params->lin_solver_paramslist[pci].Uplo;
    nrhs = params->lin_solver_paramslist[pci].nrhs;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    /* Get input matrix dimensions. */
    n= p_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    ldb = params->lin_solver_paramslist[pci].ldb;

    if(lda < n || ldb < n)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Create input matrix parameters */
    create_matrix(datatype, &A, lda, n);
    create_matrix(datatype, &A_test, lda, n);
    create_matrix(datatype, &B, ldb, nrhs);
    create_matrix(datatype, &X, n, nrhs);
    create_matrix(datatype, &B_test, ldb, nrhs);

    /* Initialize input symmetric positive definite matrix A */
    reset_matrix(datatype, n, n, A, lda);
    if (g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
        init_matrix_from_file(datatype, B, n, nrhs, ldb, g_ext_fptr);
    }
    else
    {
        /* Initialize input matrix with random numbers */
        rand_spd_matrix(datatype, &uplo, &A, n, lda);
        rand_matrix(datatype, B, n, nrhs, ldb);
    }
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);
    /* cholesky factorisation of A as input to potrs */
    invoke_potrf(&uplo, datatype, &n, A, &lda, &info);

    copy_matrix(datatype, "full", n, nrhs, B, ldb, B_test, ldb);

    /* Invoke potrs API to find x using Ax-b */
    prepare_potrs_run(&uplo, n, nrhs, A, lda, datatype, B_test, ldb, n_repeats, &time_min, &info);
    copy_matrix(datatype, "full", n, nrhs, B_test, ldb, X, n);
    /* execution time */
    *t = time_min;
    /* Compute the performance of the best experiment repeat. */
    /* (2.0)m^2 flops for Ax=b computation. */
    *perf = (double)((2.0 * n * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* Validate potrs call by computing Ax-b */
    if(info == 0)
        validate_potrs(n, nrhs, A_test, lda, X, B, ldb, datatype, residual, &vinfo);
    
    FLA_TEST_CHECK_EINFO(residual, info, einfo);

    free_matrix(A);
    free_matrix(A_test);
    free_matrix(B_test);
    free_matrix(B);
    free_matrix(X);

}

void prepare_potrs_run(char* uplo, 
    integer n,
    integer nrhs,
    void *A,
    integer lda,
    integer datatype,
    void *B,
    integer ldb,
    integer n_repeats,
    double* time_min_,
    integer* info)
{
    void *A_save = NULL, *B_test = NULL;
    double time_min = 1e9, exe_time;
    integer i;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, lda, n);
    create_matrix(datatype, &B_test, ldb, nrhs);

    *info = 0;
    for (i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
        for each iteration */
        copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);
        copy_matrix(datatype, "full", n, nrhs, B, ldb, B_test, ldb);
        exe_time = fla_test_clock();
        invoke_potrs(uplo, datatype, &n, A_save, &lda, &nrhs, B_test, &ldb, info);
        exe_time = fla_test_clock() - exe_time;
        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }
    copy_matrix(datatype, "full", n, nrhs, B_test, ldb, B, ldb);
    *time_min_ = time_min;
    free_matrix(A_save);
    free_vector(B_test);
}

void invoke_potrs(char* uplo, integer datatype,
    integer* n,
    void* A,
    integer* lda,
    integer *nrhs,
    void* B,
    integer* ldb,
    integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_spotrs(uplo, n, nrhs, A, lda, B, ldb, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dpotrs(uplo, n, nrhs, A, lda, B, ldb, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cpotrs(uplo, n, nrhs, A, lda, B, ldb, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zpotrs(uplo, n, nrhs, A, lda, B, ldb, info);
            break;
        }
    }
}
