/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

// Local prototypes.
void fla_test_geqrf_experiment(test_params_t *params, integer datatype, integer  p_cur, integer  q_cur, integer  pci, integer  n_repeats,
                                    double* perf, double* t,double* residual);
void prepare_geqrf_run(integer m_A, integer n_A, void *A, integer lda, void *T, integer datatype, integer n_repeats, double* time_min_);

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
static FILE* g_ext_fptr = NULL;

void fla_test_geqrf(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "RQ factorization";
    char* front_str = "GEQRF";
    integer tests_not_run = 1, invalid_dtype = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_geqrf_experiment);
        tests_not_run = 0;
    }
    if(argc == 9)
    {
        /* Read matrix input data from a file */
        g_ext_fptr = fopen(argv[8], "r");
        if (g_ext_fptr == NULL)
        {
            printf("\n Invalid input file argument \n");
            return;
        }
    }
    if(argc >= 8 && argc <= 9)
    {
        integer i, num_types, M,N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        g_lwork = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_geqrf_experiment(params, datatype,
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
        printf("\nIllegal arguments for geqrf\n");
        printf("./<EXE> geqrf <precisions - sdcz> <M> <N> <LDA> <LWORK> <repeats>\n");
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

void fla_test_geqrf_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double* perf,
    double* t,
    double* residual)
{
    integer m, n, lda;
    void *A = NULL, *A_test = NULL, *T = NULL;
    double time_min = 1e9;

    // Get input matrix dimensions.
    m = p_cur;
    n = q_cur;
    lda = params->lin_solver_paramslist[pci].lda;

    // Create input matrix parameters
    create_matrix(datatype, &A, lda, n);
    create_vector(datatype, &T, fla_min(m,n));

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

    // Make a copy of input matrix A. This is required to validate the API functionality.
    create_matrix(datatype, &A_test, lda, n);
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    prepare_geqrf_run(m, n, A_test, lda, T, datatype, n_repeats, &time_min);

    // execution time
    *t = time_min;

    // performance computation
    // 2mn^2 - (2/3)n^3 flops
    if(m >= n)
        *perf = (double)((2.0 * m * n * n) - (( 2.0 / 3.0 ) * n * n * n )) / time_min / FLOPS_PER_UNIT_PERF;
    else
        *perf = (double)((2.0 * n * m * m) - (( 2.0 / 3.0 ) * m * m * m )) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    // output validation
    validate_geqrf(m, n, A, A_test, lda, T, datatype, residual);

    // Free up the buffers
    free_matrix(A);
    free_matrix(A_test);
    free_vector(T);
}


void prepare_geqrf_run(integer m_A, integer n_A,
    void *A,
    integer lda,
    void *T,
    integer datatype,
    integer n_repeats,
    double* time_min_)
{
    integer min_A, i;
    void *A_save = NULL, *T_test = NULL, *work = NULL;
    integer lwork = -1, info = 0;
    double time_min = 1e9, exe_time;

    min_A = fla_min(m_A, n_A);

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, lda, n_A);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    /* Make a workspace query the first time. This will provide us with
       and ideal workspace size based on internal block size.*/
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        // call to  geqrf API
        invoke_geqrf(datatype, &m_A, &n_A, NULL, &lda, NULL, work, &lwork, &info);

        // Get work size
        lwork = get_work_value( datatype, work );

        /* Output buffers will be freshly allocated for each iterations, free up 
       the current output buffers.*/ 
        free_vector(work);
    }
    else
    {
        lwork = g_lwork;
    }

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);

        // T_test vector will hold the scalar factors of the elementary reflectors.
        create_vector(datatype, &T_test, min_A);

        // Create work buffer
        create_matrix(datatype, &work, lwork, 1);

        exe_time = fla_test_clock();

        // Call to  gerqf API
        invoke_geqrf(datatype, &m_A, &n_A, A, &lda, T_test, work, &lwork, &info);

        exe_time = fla_test_clock() - exe_time;

        // Get the best execution time
        time_min = fla_min(time_min, exe_time);

        // Make a copy of the output buffers. This is required to validate the API functionality.
        copy_vector(datatype, min_A, T_test, 1, T, 1);

        // Free up the output buffers
        free_vector(work);
        free_vector(T_test);
    }

    *time_min_ = time_min;

    free_matrix(A_save);
}


void invoke_geqrf(integer datatype, integer* m, integer* n, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sgeqrf_(m, n, a, lda, tau, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            dgeqrf_(m, n, a, lda, tau, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            cgeqrf_(m, n, a, lda, tau, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zgeqrf_(m, n, a, lda, tau, work, lwork, info);
            break;
        }
    }
}
