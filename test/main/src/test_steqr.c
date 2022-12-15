/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_steqr_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
integer n_repeats, double* perf, double* t, double* residual);
void prepare_steqr_run(char* compz, integer n, void* Z, integer ldz, void* D, void* E, integer datatype, integer n_repeats, double* time_min_, integer* info);
void invoke_steqr(integer datatype, char* compz, integer* n, void* z, integer* ldz, void* d, void* e, void* work, integer* info);
static FILE* g_ext_fptr = NULL;

void fla_test_steqr(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Eigen Decomposition of symmetrix tridiagonal matrix";
    char* front_str = "STEQR";
    integer tests_not_run = 1, invalid_dtype = 0;

    if(argc == 1)
    {
        /* Test with parameters from config */
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_steqr_experiment);
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
        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->eig_sym_paramslist[0].threshold_value = CLI_NORM_THRESH;
            params->eig_sym_paramslist[0].uplo = 'L';

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
                fla_test_steqr_experiment(params, datatype,
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
        printf("Invalid arguments for STEQR\n");
        printf("Usage: ./<EXE> steqr <precisions - sdcz> <COMPZ> <N> <LDZ> <repeats>\n");
    }
    else if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n");
    }
    if (g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
    }
    return;
}

void fla_test_steqr_experiment(test_params_t *params,
                               integer  datatype,
                               integer  p_cur,
                               integer  q_cur,
                               integer pci,
                               integer n_repeats,
                               double* perf,
                               double *time_min,
                               double* residual)
{
    integer n, ldz, lda, info = 0, vinfo = 0;
    char compz, uplo;
    void *Z = NULL, *Z_test = NULL, *A = NULL, *Q = NULL;
    void *D = NULL, *D_test = NULL, *E = NULL, *E_test = NULL;

    /* Get input matrix dimensions.*/
    compz = params->eig_sym_paramslist[pci].compz;
    *residual = params->eig_sym_paramslist[pci].threshold_value;
    uplo = params->eig_sym_paramslist[pci].uplo;

    n = p_cur;
    ldz = params->eig_sym_paramslist[pci].ldz;
    lda = ldz;

    if(ldz < n || lda < n)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Create input matrix parameters */
    create_matrix(datatype, &Z, ldz, n);
    create_matrix(datatype, &A, lda, n);
    create_matrix(datatype, &Q, ldz, n);

    reset_matrix(datatype, n, n, Z, ldz);
    reset_matrix(datatype, n, n, A, lda);
    reset_matrix(datatype, n, n, Q, ldz);

    create_vector(get_realtype(datatype), &D, n);
    create_vector(get_realtype(datatype), &E, n-1);

    if (g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
    }
    else
    {
        /* input matrix Z with random symmetric numbers and D,E matrix with diagonal and subdiagonal values */
        if (datatype == FLOAT || datatype == DOUBLE)
        {
            rand_sym_matrix(datatype, A, n, n, lda);
        }
        else
        {
            rand_hermitian_matrix(datatype, n, &A, lda);
        }
    }

    copy_matrix(datatype, "full", n, n, A, lda, Q, ldz);
    /* Make a copy of input matrix Z. This is required to validate the API functionality.*/
    create_matrix(datatype, &Z_test, ldz, n);
    reset_matrix(datatype, n, n, Z_test, ldz);

    invoke_sytrd(datatype, &uplo, compz, n, Q, lda, D, E, &info);
    if(info < 0)
    {
        *residual = DBL_MAX;
        free_matrix(Z_test);
        free_matrix(A);
        free_matrix(Z);
        free_matrix(Q);
        free_vector(D);
        free_vector(E);
        return;
    }
    /*form tridiagonal matrix Z by copying from matrix*/
    copy_sym_tridiag_matrix(datatype, D, E, n, n, Z, ldz);

    if(compz == 'I')
        set_identity_matrix(datatype, n, n, Z_test, ldz);
    else
    {
        copy_matrix(datatype, "full", n, n, Q, lda, Z_test, ldz);
        copy_matrix(datatype, "full", n, n, A, lda, Z, ldz);
    }
    create_vector(get_realtype(datatype), &D_test, n);
    create_vector(get_realtype(datatype), &E_test, n-1);
    copy_vector(get_realtype(datatype), n, D, 1, D_test, 1);
    copy_vector(get_realtype(datatype), n-1, E, 1, E_test, 1);

    prepare_steqr_run(&compz, n, Z_test, ldz, D_test, E_test, datatype, n_repeats, time_min, &info);

    /* performance computation
       24 n^2 flops for eigen vectors of Z, compz = 'N'
       7 n^3 flops for eigen vectors of Z, compz = 'V' or 'I'
       14 n^3 flops for eigen vectors of Z for complex, compz = 'V' or 'I' */
 
    if( compz == 'I' || compz == 'V')
        *perf = (double)(7.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    else if( compz == 'N')
        *perf = (double)(24.0 * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf = (double)(14.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;

    /* output validation */
    if (info == 0)
        validate_syevd(&compz, n, Z, Z_test, ldz, D_test, datatype, residual, &vinfo);
    
    if (info < 0 || vinfo < 0)
        *residual = DBL_MAX;

    /* Free up the buffers */
    free_matrix(Z);
    free_vector(D);
    free_vector(E);
    free_matrix(A);
    free_matrix(Q);
    free_matrix(Z_test);
    free_vector(D_test);
    free_vector(E_test);
}

void prepare_steqr_run(char *compz,
                       integer n,
                       void *Z,
                       integer ldz,
                       void *D,
                       void *E,
                       integer datatype,
                       integer n_repeats,
                       double* time_min_,
                       integer* info)
{
    void *Z_save = NULL, *D_save = NULL, *E_save = NULL, *work = NULL;
    integer i;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &Z_save, ldz, n);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);

    create_vector(get_realtype(datatype), &D_save, n);
    create_vector(get_realtype(datatype), &E_save, n-1);
    copy_vector(get_realtype(datatype), n, D, 1, D_save, 1);
    copy_vector(get_realtype(datatype), n-1, E, 1, E_save, 1);

    for (i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);
        copy_vector(get_realtype(datatype), n, D_save, 1, D, 1);
        copy_vector(get_realtype(datatype), n-1, E_save, 1, E, 1);

        create_vector(get_realtype(datatype), &work, 2 * (n - 1));
        exe_time = fla_test_clock();

        /* call to API */
        invoke_steqr(datatype, compz, &n, Z, &ldz, D, E, work, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
    }

    *time_min_ = time_min;

    free_matrix(Z_save);
    free_matrix(D_save);
    free_matrix(E_save);
}

void invoke_steqr(integer datatype, char* compz, integer* n, void* z, integer* ldz, void* d, void* e, void* work, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            ssteqr_(compz, n, d, e, z, ldz, work, info);
            break;
        }
        case DOUBLE:
        {
            dsteqr_(compz, n, d, e, z, ldz, work, info);
            break;
        }
        case COMPLEX:
        {
            csteqr_(compz, n, d, e, z, ldz, work, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zsteqr_(compz, n, d, e, z, ldz, work, info);
            break;
        }
    }
}
