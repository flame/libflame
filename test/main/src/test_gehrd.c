/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"


/* Local prototypes */
void fla_test_gehrd_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, integer einfo, double* perf, double* t, double* residual);
void prepare_gehrd_run(integer n, integer* ilo, integer* ihi, void* A, integer lda, void *tau, integer datatype,
                        integer n_repeats, double* time_min_, integer* info);
void invoke_gehrd(integer datatype, integer* n, integer* ilo, integer* ihi, void* a, integer* lda, void *tau, void* work,
                    integer* lwork, integer* info);
static FILE* g_ext_fptr = NULL;

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
void fla_test_gehrd(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Reduces matrix to upper hessenberg from";
    char* front_str = "GEHRD";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_gehrd_experiment);
        tests_not_run = 0;
    }
    if(argc == 10)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[9]);
    }
    if(argc >=9 && argc <= 10)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype,type_flag[4] = {0};
        char *endptr;

        /* Prase the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ilo = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ihi = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        g_lwork = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_gehrd_experiment(params, datatype,
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
        printf("\nIllegal arguments for GEHRD\n");
        printf("./<EXE> gehrd <precisions - sdcz> <N> <ILO> <IHI> <LDA> <LWORK> <repeats>\n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
    if (g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
    }
}

void fla_test_gehrd_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    integer  einfo,
    double   *perf,
    double   *time_min,
    double   *residual)
{
    integer n, lda;
    integer ilo, ihi, info = 0, vinfo = 0;
    void *A = NULL, *A_Test = NULL, *tau = NULL;

    /* Get input matrix dimensions. */
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;

    if(lda < n)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Initialize parameter needed for gehrd() call. */
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    ilo = params->lin_solver_paramslist[pci].ilo;
    ihi = params->lin_solver_paramslist[pci].ihi;

    /* Create input matrix parameters*/
    create_matrix(datatype, &A, lda, n);
    create_vector(datatype, &tau, n-1);

    if(g_ext_fptr != NULL)
    {
        init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
    }
    else
    {
        /* Intialize matrix H with ILO and IHI conditions to generate hessenberg matrix */
        get_generic_triangular_matrix(datatype, n, A, lda, ilo, ihi);
    }

    /* Make copy of matrix A. This is required to validate the API functionality */
    create_matrix(datatype, &A_Test, lda, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_Test, lda);

    prepare_gehrd_run(n, &ilo, &ihi, A_Test, lda, tau, datatype, n_repeats, time_min, &info);

    /* Performance computation
       (2/3)*(ihi - ilo)^2(2ihi + 2ilo + 3n) flops for real values
       4*((2/3)*(ihi - ilo)^2(2ihi + 2ilo + 3n)) flops for complex values */

    if(datatype == FLOAT || datatype == DOUBLE)
        *perf = (double)((2.0 / 3.0) * pow((ihi - ilo),2) * (2*ihi + 2*ilo + 3*n)) / *time_min / FLOPS_PER_UNIT_PERF;
    else
        *perf = (double)(4.0 * ((2.0 / 3.0) * pow((ihi - ilo),2) * (2*ihi + 2*ilo + 3*n))) / *time_min / FLOPS_PER_UNIT_PERF;

    /* Output Validation */
    if(info == 0)
        validate_gehrd(n, ilo, ihi, A, A_Test, lda, tau, datatype, residual, &vinfo);

    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    
    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_Test);
    free_vector(tau);
}

void prepare_gehrd_run(integer n, integer* ilo, integer* ihi, void* A, integer lda, void* tau,
                        integer datatype, integer n_repeats, double* time_min_, integer* info)
{
    void *A_save = NULL, *work = NULL, *tau_test = NULL;
    integer i, lwork;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in each itertaion.*/
    create_matrix(datatype, &A_save, lda, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.*/
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* call to  gehrd API */
        invoke_gehrd(datatype ,&n, ilo, ihi, NULL, &lda, NULL, work, &lwork, info);

        /* Output buffers will be freshly allocated for each iterations, free up
        the current output buffers.*/
        if(*info == 0)
        {
            /* Get work size */
            lwork = get_work_value( datatype, work );
        }

        free_vector(work);
    }
    else
    {
        lwork = g_lwork;
    }

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix H and Z value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
        create_vector(datatype, &work, lwork);
        create_vector(datatype, &tau_test, n-1);
        exe_time = fla_test_clock();

        /* Call to gehrd API */
        invoke_gehrd(datatype, &n, ilo, ihi, A, &lda, tau_test, work, &lwork, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        copy_vector(datatype, n-1, tau_test, 1, tau, 1);

        /* Free up the output buffers */
        free_vector(work);
    }
    *time_min_ = time_min;

    free_matrix(A_save);
    free_vector(tau_test);
}

void invoke_gehrd(integer datatype, integer* n, integer* ilo, integer* ihi, void* A, integer* lda, void *tau, void* work, integer* lwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgehrd(n, ilo, ihi, A, lda, tau, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgehrd(n, ilo, ihi, A, lda, tau, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgehrd(n, ilo, ihi, A, lda, tau, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgehrd(n, ilo, ihi, A, lda, tau, work, lwork, info);
            break;
        }
    }
}
