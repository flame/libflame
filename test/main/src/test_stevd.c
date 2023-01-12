/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_stevd_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
integer n_repeats, double* perf, double* t, double* residual);
void prepare_stevd_run(char* jobz, integer n, void* Z, integer ldz, void* D, void* E, integer datatype, integer n_repeats, double* time_min_);
void invoke_stevd(integer datatype, char* jobz, integer* n, void* z, integer* ldz, void* d, void* e, void* work, integer* lwork, void* iwork, integer* liwork, integer* info);

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
static integer g_liwork;
static FILE* g_ext_fptr = NULL;

void fla_test_stevd(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Eigen Decomposition of symmetrix tridiagonal matrix";
    char* front_str = "STEVD";
    integer tests_not_run = 1, invalid_dtype = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        g_liwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_stevd_experiment);
        tests_not_run = 0;
    }
    if(argc == 10)
    {
        /* Read matrix input data from a file */
        g_ext_fptr = fopen(argv[9], "r");
        if (g_ext_fptr == NULL)
        {
            printf("\n Invalid input file argument \n");
            return;
        }
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
        params->eig_sym_paramslist[0].jobz = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);

        g_lwork = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        g_liwork = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_stevd_experiment(params, datatype,
                                          N, N,
                                          0,
                                          n_repeats,
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
        printf("Invalid arguments for stevd\n");
        printf("Usage: ./<EXE> stevd <precisions - sdcz> <JOBZ> <N> <LDZ> <LWORk> <LIWORK> <repeats>\n");
    }
    else if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n");
    }
    if(g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
    }
    return;
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
    ldz = params->eig_sym_paramslist[pci].ldz;
    *residual = params->eig_sym_paramslist[pci].threshold_value;

    if(datatype == FLOAT || datatype == DOUBLE)
    {
        n = p_cur;

        /* Create input matrix parameters */
        create_matrix(datatype, &Z, ldz, n);
        reset_matrix(datatype, n, n, Z, ldz);
        create_vector(datatype, &D, n);
        create_vector(datatype, &E, n-1);

        if (g_ext_fptr != NULL)
        {
            /* Initialize input matrix with custom data */
            init_matrix_from_file(datatype, Z, n, n, ldz, g_ext_fptr);
        }
        else
        {
            /* input matrix Z with random symmetric numbers and D,E matrix with diagonal and subdiagonal values */
            rand_sym_tridiag_matrix(datatype, Z, n, n, ldz);
        }

        get_diagonal(datatype, Z, n, n, ldz, D);
        get_subdiagonal(datatype, Z, n, n, ldz, E);

        /* Make a copy of input matrix A. This is required to validate the API functionality.*/
        create_matrix(datatype, &Z_test, ldz, n);
        reset_matrix(datatype, n, n, Z_test, ldz);
        create_vector(datatype, &D_test, n);
        create_vector(datatype, &E_test, n-1);
        copy_vector(datatype, n, D, 1, D_test, 1);
        copy_vector(datatype, n-1, E, 1, E_test, 1);

        prepare_stevd_run(&jobz, n, Z_test, ldz, D_test, E_test, datatype, n_repeats, time_min);

        /* performance computation
           6 * n^3 + n^2 flops for eigen vectors
           6 * n^2 flops for eigen values */
        if( jobz == 'V')
            *perf = (double)((6.0 * n * n * n) + (n * n)) / *time_min / FLOPS_PER_UNIT_PERF;
        else
            *perf = (double)(6.0 * n * n) / *time_min / FLOPS_PER_UNIT_PERF;

        /* output validation */
        validate_syevd(&jobz, n, Z, Z_test, ldz, D_test, datatype, residual);

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
                       integer ldz,
                       void *D,
                       void *E,
                       integer datatype,
                       integer n_repeats,
                       double* time_min_)
{
    void *Z_save, *D_save, *E_save, *work, *iwork;
    integer lwork, liwork;
    integer i;
    integer info = 0;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &Z_save, ldz, n);
    create_vector(datatype, &D_save, n);
    create_vector(datatype, &E_save, n-1);

    copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);
    copy_vector(datatype, n, D, 1, D_save, 1);
    copy_vector(datatype, n-1, E, 1, E_save, 1);

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size.*/
    if(g_lwork <= 0 || g_liwork <= 0)
    {
        lwork = -1;
        liwork = -1;
        create_vector(INTEGER, &iwork, 1);
        create_vector(datatype, &work, 1);
        /* call to  stevd API */
        invoke_stevd(datatype, jobz, &n, NULL, &ldz, NULL, NULL, work, &lwork, iwork, &liwork, &info);

        /* Get work size */
        lwork = get_work_value(datatype, work );
        liwork = get_work_value(INTEGER, iwork );
        /* Output buffers will be freshly allocated for each iterations, free up
        the current output buffers.*/
        free_vector(work);
        free_vector(iwork);
    }
    else
    {
        lwork = g_lwork;
        liwork = g_liwork;
    }

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
        time_min = fla_min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(iwork);

    }

    *time_min_ = time_min;

    free_matrix(Z_save);
    free_matrix(D_save);
    free_matrix(E_save);
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
