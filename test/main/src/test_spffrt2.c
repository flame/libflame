/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_spffrt2_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_spffrt2_run(integer n_A, integer ncolm, integer pn, void *A, integer datatype, integer n_repeats, double* time_min_);
void invoke_spffrt2(integer datatype, void *a, integer* n, integer * ncolm, void *work, void *work2);
static FILE* g_ext_fptr = NULL;

void fla_test_spffrt2(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Computes LDLT partial factorization";
    char* front_str = "SPFFRT2";
    integer tests_not_run = 1, invalid_dtype = 0;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_spffrt2_experiment);
        tests_not_run = 0;
    }
    if (argc == 7)
    {
        /* Read matrix input data from a file */
        g_ext_fptr = fopen(argv[6], "r");
        if (g_ext_fptr == NULL)
        {
            printf("\n Invalid input file argument \n");
            return;
        }
    }
    if (argc >= 6 && argc <= 7)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ncolm = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_spffrt2_experiment(params, datatype,
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
        printf("\nIllegal arguments for spffrt2 \n");
        printf("./<EXE> spffrt2 <precisions - sdcz>  <N> <ncolm> <repeats> [file]\n");
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


void fla_test_spffrt2_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    double* perf,
    double* t,
    double* residual)
{
    integer n, ncolm, pn;
    void *A, *AP;
    double time_min = 1e9;

    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    ncolm = params->lin_solver_paramslist[pci].ncolm;

    /* Determine the dimensions*/    
    n = p_cur;
    pn = n * (n + 1) / 2;
    /* Create the matrices for the current operation*/
    create_matrix(datatype, &A, n, n);
    create_vector(datatype, &AP, pn);
    if (g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, n, n, n, g_ext_fptr);
    }
    else
    {
        /* Initialize input matrix with random numbers */
        rand_sym_matrix(datatype, A, n, n, n);
    }
    /* Pack a symmetric matrix in column first order */
    pack_matrix_lt(datatype, A, AP, n, n);

    /* call to API */
    prepare_spffrt2_run(n, ncolm, pn, AP, datatype, n_repeats, &time_min);

    /* execution time */
    *t = time_min;

    /* performance computation */
    if (datatype == FLOAT || datatype == DOUBLE)
        *perf = (double)ncolm / 6.0f * (2.0 * ncolm * ncolm - 6.0 * ncolm * n + 3.0 * ncolm + 6.0 * n * n - 6.0 * n + 7) / time_min / FLOPS_PER_UNIT_PERF;
    else if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf = (double)ncolm / 3.0f * (4.0 * ncolm * ncolm - 12.0 * ncolm * n + 9.0 * ncolm + 12.0 * n * n - 18.0 * n + 8) / time_min / FLOPS_PER_UNIT_PERF;

    /* output validation */
    if(ncolm <= n && n > 0 && ncolm > 0)
    {
        validate_spffrt2(n, ncolm, A, AP, datatype, residual);
    }
    else
    {   /* Assigning bigger value to residual as execution fails */
        *residual = DBL_MAX;
    }

    /* Free up the buffers */
    free_matrix(A);
    free_vector(AP);
}


void prepare_spffrt2_run(integer n_A,
    integer ncolm,
    integer pn,
    void* AP,
    integer datatype,
    integer n_repeats,
    double* time_min_)
{
    integer i;
    void *AP_save, *work = NULL, *work2 = NULL;
    double time_min = 1e9, exe_time;

    create_vector(datatype, &AP_save, pn);
    create_vector(datatype, &work, 2*n_A);
    create_vector(datatype, &work2, 2*n_A);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Copy original input data */
        copy_vector(datatype, pn, AP, i_one, AP_save, i_one);

        exe_time = fla_test_clock();

        /*  call  spffrt2 API with AFACT to get A INV */
        invoke_spffrt2(datatype, AP_save, &n_A, &ncolm, work, work2);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    /*  Save the final result to A matrix*/
    copy_vector(datatype, pn, AP_save, i_one, AP, i_one);
    free_vector(AP_save);
    free_vector(work);
    free_vector(work2);
}


/*
 *  spffrt2_API calls LAPACK interface
 *  */
void invoke_spffrt2(integer datatype, void *ap, integer *n, integer *ncolm, void* work, void* work2)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sspffrt2(ap, n, ncolm, work, work2);
            break;
        }
        
        case DOUBLE:
        {
            fla_lapack_dspffrt2(ap, n, ncolm, work, work2);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cspffrt2(ap, n, ncolm, work, work2);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zspffrt2(ap, n, ncolm, work, work2);
            break;
        }
    }
}

