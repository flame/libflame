/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_getrf_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_getrf_run(integer m_A, integer n_A, void *A, integer lda, integer* ipiv, integer datatype, integer n_repeats, double* time_min_, integer *info);
void invoke_getrf(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *ipiv, integer *info);
static FILE* g_ext_fptr = NULL;

void fla_test_getrf(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "LU factorization";
    char* front_str = "GETRF";
    integer tests_not_run = 1, invalid_dtype = 0;
    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_getrf_experiment);
        tests_not_run = 0;
    }
    if (argc == 8)
    {
        /* Read matrix input data from a file */
        g_ext_fptr = fopen(argv[7], "r");
        if (g_ext_fptr == NULL)
        {
            printf("\n Invalid input file argument \n");
            return;
        }
    }
    if (argc >= 7 && argc <= 8)
    {
        integer i, num_types, M,N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype,type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_getrf_experiment(params, datatype,
                                          M, N,
                                          0,
                                          n_repeats,
                                          &perf, &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str,
                                      stype,
                                     RECT_INPUT,
                                      M, N,
                                      residual, params->lin_solver_paramslist[0].solver_threshold,
                                      time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for getrf\n");
        printf("./<EXE> getrf <precisions - sdcz> <M> <N> <LDA> <repeats>\n");
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

void fla_test_getrf_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    double* perf,
    double* t,
    double* residual)
{
    integer m, n, lda;
    integer info = 0, vinfo = 0;
    void* IPIV;
    void *A, *A_test;
    double time_min = 1e9;

    /* Determine the dimensions*/
    m = p_cur;
    n = q_cur;
    lda = params->lin_solver_paramslist[pci].lda;

    if(lda < m)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Create the matrices for the current operation*/
    create_matrix(datatype, &A, lda, n);
    create_vector(INTEGER, &IPIV, fla_min(m, n));

    /* Initialize the test matrices*/
    if (g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, m, n, lda, g_ext_fptr);
    }
    else
    {
        /* Initialize input matrix with random numbers */
        rand_matrix(datatype, A, m, n, lda);
    }

    /* Save the original matrix*/
    create_matrix(datatype, &A_test, lda, n);
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    /* call to API */
    prepare_getrf_run(m, n, A_test, lda, IPIV, datatype, n_repeats, &time_min, &info);

    /* execution time */
    *t = time_min;

    /* performance computation */
    if(m == n)
    {
        *perf = (2.0 / 3.0) * n * n * n / time_min / FLOPS_PER_UNIT_PERF;
    }
    else if(m > n)
    {
        *perf = (1.0 / 3.0) * n * n * (3 * m - n) / time_min / FLOPS_PER_UNIT_PERF;
    }
    else
    {
        *perf = (1.0 / 3.0) * m * m * (3 * n - m) / time_min / FLOPS_PER_UNIT_PERF;
    }
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if (info == 0)
       validate_getrf(m, n, A, A_test, lda, IPIV, datatype, residual, &vinfo);

    /* Assigning bigger value to residual as execution fails */
    if (info < 0 || vinfo < 0)
        *residual = DBL_MAX;

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(IPIV);
}

void prepare_getrf_run(integer m_A,
    integer n_A,
    void* A,
    integer lda,
    integer* IPIV,
    integer datatype,
    integer n_repeats,
    double* time_min_,
    integer* info)
{
    integer i;
    void *A_save;
    double time_min = 1e9, exe_time;

    /* Save the original matrix */
    create_matrix(datatype, &A_save, lda, n_A);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    for (i = 0; i < n_repeats && *info == 0; ++i)
    {

        /* Copy original input data */
        copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

        exe_time = fla_test_clock();

        /*  call to API */
        invoke_getrf(datatype, &m_A, &n_A, A_save, &lda, IPIV, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    /*  Save the AFACT to matrix A */
    copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);
    free_matrix(A_save);

}

/*
 *  GETRF_API calls LAPACK interface of
 *  Singular value decomposition - gesvd
 *  */
void invoke_getrf(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *ipiv, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sgetrf_(m, n, a, lda, ipiv, info);
            break;
        }

        case DOUBLE:
        {
            dgetrf_(m, n, a, lda, ipiv, info);
            break;
        }

        case COMPLEX:
        {
            cgetrf_(m, n, a, lda, ipiv, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zgetrf_(m, n, a, lda, ipiv, info);
            break;
        }
    }
}
