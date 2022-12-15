/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_potrf_experiment(test_params_t *params, integer datatype, integer  p_cur, integer  q_cur, integer  pci, integer  n_repeats,double* perf, double* time_min, double* residual);
void prepare_potrf_run(char* uplo, integer m, void *A, integer lda, integer datatype, integer n_repeats, double* time_min_, integer *info);
static FILE* g_ext_fptr = NULL;

void fla_test_potrf(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Cholesky factorization";
    char* front_str = "POTRF";
    integer tests_not_run = 1, invalid_dtype = 0;
    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_potrf_experiment);
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
        integer i, num_types,N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype,type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].Uplo = argv[3][0];
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
                fla_test_potrf_experiment(params, datatype,
                                          N, N,
                                          0,
                                          n_repeats,
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
        printf("\nIllegal arguments for potrf\n");
        printf("./<EXE> potrf <precisions - sdcz> <Uplo> <N> <LDA> <repeats>\n");
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

void fla_test_potrf_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double* perf,
    double* time_min,
    double* residual)
{
    integer m, lda;
    integer info = 0, vinfo = 0;
    void *A = NULL, *A_test = NULL;
    char uplo = params->lin_solver_paramslist[pci].Uplo;

    /* Get input matrix dimensions */
    m = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;

    if(lda < m)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Create input matrix parameters */
    create_matrix(datatype, &A, lda, m);

    if (g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, m, m, lda, g_ext_fptr);
    }
    else
    {
        /* Initialize input matrix with random numbers */
        rand_spd_matrix(datatype, &uplo, &A, m, lda);
    }

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, &A_test, lda, m);
    copy_matrix(datatype, "full", m, m, A, lda, A_test, lda);

    prepare_potrf_run(&uplo, m, A_test, lda, datatype, n_repeats, time_min, &info);

    /* Compute the performance of the best experiment repeat */
    /* (1/3)m^3 for real and (4/3)m^3 for complex*/
    *perf = (double)(1.0 / 3.0 * m * m * m) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    if(info == 0)
        validate_potrf(&uplo, m, A, A_test, lda, datatype, residual, &vinfo);

    /* Assigning bigger value to residual as execution fails */
    if (info < 0 || vinfo < 0)
        *residual = DBL_MAX;

    free_matrix(A);
    free_matrix(A_test);
}

void prepare_potrf_run(char* uplo, integer m,
    void *A,
    integer lda,
    integer datatype,
    integer n_repeats,
    double* time_min_,
    integer* info)
{
    void *A_save = NULL;
    double time_min = 1e9, exe_time;
    integer i;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, lda, m);
    copy_matrix(datatype, "full", m, m, A, lda, A_save, lda);

    for (i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
        for each iteration */
        copy_matrix(datatype, "full", m, m, A_save, lda, A, lda);
        exe_time = fla_test_clock();
        invoke_potrf(uplo, datatype, &m, A, &lda, info);
        exe_time = fla_test_clock() - exe_time;
        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    free_matrix(A_save);
}

void invoke_potrf(char* uplo, integer datatype, integer* m, void* a, integer* lda, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            spotrf_(uplo, m, a, lda, info);
            break;
        }
        case DOUBLE:
        {
            dpotrf_(uplo, m, a, lda, info);
            break;
        }
        case COMPLEX:
        {
            cpotrf_(uplo, m, a, lda, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zpotrf_(uplo, m, a, lda, info);
            break;
        }
    }
}
