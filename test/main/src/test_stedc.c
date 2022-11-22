/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes. */
void fla_test_stedc_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer  pci,
                                    integer  n_repeats, double* perf, double* t, double* residual);
void prepare_stedc_run(char* compz, integer n, void* D, void* E, void* Z,
                      integer ldz, integer datatype, integer n_repeats, double* time_min_);
void invoke_stedc(integer datatype, char* compz, integer* n, void* D, void* E, void* Z,
                  integer* ldz, void* work, integer* lwork, void* rwork,
                  integer* lrwork, integer* iwork, integer* liwork, integer *info);

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
static integer g_liwork;
static integer g_lrwork;
void fla_test_stedc(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Eigenvalues/eigenvectors of symmetric tridiagonal matrix";
    char* front_str = "STEDC";
    integer tests_not_run = 1, invalid_dtype = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        g_liwork = -1;
        g_lrwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_stedc_experiment);
        tests_not_run = 0;
    }
    else if(argc == 10)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].compz = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);

        g_lwork = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        g_liwork = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        g_lrwork = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_stedc_experiment(params, datatype,
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
        printf("Invalid arguments for stedc\n");
        printf("Usage: ./<EXE> stedc <precisions - sdcz> <COMPZ> <N> <LDZ> <LWORk> <LIWORK> <LRWORK> <repeats>\n");
    }
    else if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n");
    }
    return;
}

void fla_test_stedc_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double* perf,
    double* t,
    double* residual)
{
    integer n, info = -1, realtype;
    void *D = NULL, *D_test = NULL, *E = NULL, *E_test = NULL, *Z_test = NULL;
    void *Z_input = NULL, *A = NULL;
    double time_min = 1e9;
    char compz, uplo;
    
    /* Get input matrix dimensions. */
    n = q_cur;
    
    /* Initialize parameter needed for STEDC() call. */
    compz = params->eig_sym_paramslist[pci].compz;
    create_matrix(datatype, &A, n, n);
    
    realtype = get_realtype(datatype);
    create_vector(realtype, &D, n);
    create_vector(realtype, &E, n-1);
    
    /* Create random symmetric/hermitian matrix if compz = V. */
    if (compz == 'V') {
        if ((datatype == FLOAT) || (datatype == DOUBLE)) {
            rand_sym_matrix(datatype, A, n, n, n);
        } else {
            rand_hermitian_matrix(datatype, n, &A, n);
        }
    } else { /* Create tridiagonal matrix using random Diagonal, subdiagonal elements if compz != V. */
        rand_vector(realtype, D, n, 1);
        rand_vector(realtype, E, n-1, 1);
        copy_sym_tridiag_matrix(datatype, D, E, n, n, A, n);
    }
    create_matrix(datatype, &Z_input, n, n);
    copy_matrix(datatype, "full", n, n, A, n, Z_input, n);

    /* Call SYTRD(), ORGTR() to get tridiagonal/orthogonal matrix when compz = V. */
    if (compz == 'V') {
        /* Initialize parameter needed for SYTRD() call. */
        uplo = 'U';
        /* Call SYTRD() orthogonal matrix and tridiagonal elements.
           invoke_sytrd() internally calls ORGTR() to get orthogonal matrix.*/
        invoke_sytrd(datatype, &uplo, compz, n, A, n, D, E, &info);
        if (info != 0) {
            free_matrix(A);
            free_vector(D);
            free_vector(E);
            free_matrix(Z_input);
            return;
        }
    }
    /* Make a copy of input matrices. This is required to validate the API functionality. */
    create_matrix(datatype, &Z_test, n, n);
    if (compz == 'V') {
        copy_matrix(datatype, "full", n, n, A, n, Z_test, n);
    }
    create_vector(realtype, &D_test, n);
    copy_vector(realtype, n, D, 1, D_test, 1);
    create_vector(realtype, &E_test, n-1);
    copy_vector(realtype, n-1, E, 1, E_test, 1);

    prepare_stedc_run(&compz, n, D_test, E_test, Z_test, n, datatype, n_repeats, &time_min);
    
    /* Execution time. */
    *t = time_min;

    /* Performance computation
       (6)n^3 flops for eigen vectors
       (4/3)n^3 flops for eigen values. */
    *perf = (double)((4.0 / 3.0) * n * n * n) / *t / FLOPS_PER_UNIT_PERF;
    if (compz != 'N') {
        *perf += (double)(6 * n * n * n) / *t / FLOPS_PER_UNIT_PERF;
    }
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX) {
        *perf *= 2.0;
    }

    /* Output validation. */
    if (compz != 'N') {
        validate_stedc(compz, n, D_test, Z_input, Z_test, datatype, residual);
    } else {
        *residual = 0.0;
    }

    /* Free up buffers. */
    free_matrix(Z_input);
    free_matrix(Z_test);
    free_vector(D_test);
    free_vector(E_test);
    free_matrix(A);
    free_vector(D);
    free_vector(E);
}

void prepare_stedc_run(char* compz, integer n, void* D, void* E, void* Z,
                      integer ldz, integer datatype, integer n_repeats, double* time_min_)
{
    integer index, info = 0, lwork, liwork, lrwork, realtype;
    void *D_save = NULL, *E_save = NULL, *E_test = NULL, *Z_save = NULL;
    void *work = NULL, *iwork = NULL, *rwork = NULL;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrices. Same input values will be passed in
       each itertaion.*/
    if (*compz == 'V') {
        create_matrix(datatype, &Z_save, ldz, n);
        copy_matrix(datatype, "full", ldz, n, Z, ldz, Z_save, ldz);
    }
    realtype = get_realtype(datatype);
    create_vector(realtype, &D_save, n);
    copy_vector(realtype, n, D, 1, D_save, 1);
    create_vector(realtype, &E_save, n-1);
    copy_vector(realtype, n-1, E, 1, E_save, 1);
    
    /* Call to STEDC() API to get work buffers size. */
    if(g_lwork <= 0 || g_liwork <= 0 || ((datatype == COMPLEX || datatype == DOUBLE_COMPLEX) && g_lrwork <= 0))
    {
        /* Make a workspace query the first time. This will provide us with
        and ideal workspace size based on internal block size.*/
        create_vector(datatype, &work, 1);
        create_vector(realtype, &rwork, 1);
        create_vector(INTEGER, &iwork, 1);
        lwork = -1;
        liwork = -1;
        lrwork = -1;
        invoke_stedc(datatype, compz, &n, D, E_test, Z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

        /* Get work buffers size. */
        if (info == 0) 
        {
            lwork = get_work_value(datatype, work );
            liwork = get_work_value(INTEGER, iwork );
            lrwork = get_work_value(realtype, rwork );      
            free_vector(work);
            free_vector(rwork);
            free_vector(iwork);
        } 
        else 
        {
            free_vector(work);
            free_vector(rwork);
            free_vector(iwork);
            return;
        }
    }
    else
    {   
        lwork = g_lwork;
        liwork = g_liwork;
        lrwork = g_lrwork;
    }
    create_vector(datatype, &work, lwork);
    if ((datatype == COMPLEX) || (datatype == DOUBLE_COMPLEX)) {
        create_vector(realtype, &rwork, lrwork);
    }
    create_vector(INTEGER, &iwork, liwork);
    create_vector(realtype, &E_test, n-1);

    for (index = 0; index < n_repeats; ++index)
    {
        /* Restore input matrices and allocate memory to output buffers
           for each iteration. */
        if (*compz == 'V') {
            copy_matrix(datatype, "full", ldz, n, Z_save, ldz, Z, ldz);
        }
        copy_vector(realtype, n, D_save, 1, D, 1);
        copy_vector(realtype, n-1, E_save, 1, E_test, 1);
        reset_vector(datatype, work, lwork, 1);
        if ((datatype == COMPLEX) || (datatype == DOUBLE_COMPLEX)) {
            reset_vector(realtype, rwork, lrwork, 1);
        }
        reset_vector(INTEGER, iwork, liwork, 1);
        
        exe_time = fla_test_clock();

        /* Call to STEDC() API. */
        invoke_stedc(datatype, compz, &n, D, E_test, Z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time. */
        time_min = min(time_min, exe_time);
    }
    *time_min_ = time_min;

    /* Free up buffers. */
    if (*compz == 'V') {
        free_vector(Z_save);
    }
    free_vector(D_save);
    free_vector(E_save);
    free_vector(E_test);
    free_vector(work);
    if ((datatype == COMPLEX) || (datatype == DOUBLE_COMPLEX)) {
        free_vector(rwork);
    }
    free_vector(iwork);
}

void invoke_stedc(integer datatype, char* compz, integer* n, void* D, void* E, void* Z, integer* ldz, void* work, integer* lwork, void* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sstedc_(compz, n, D, E, Z, ldz, work, lwork, iwork, liwork, info);
            break;
        }
        case DOUBLE:
        {
            dstedc_(compz, n, D, E, Z, ldz, work, lwork, iwork, liwork, info);
            break;
        }
        case COMPLEX:
        {
            cstedc_(compz, n, D, E, Z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zstedc_(compz, n, D, E, Z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
            break;
        }
    }
}
