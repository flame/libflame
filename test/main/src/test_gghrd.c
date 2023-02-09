/*
    Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"


/* Local prototypes */
void fla_test_gghrd_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_gghrd_run(char* compq, char* compz, integer n, integer* ilo, integer* ihi, void* a, integer lda,
                            void* b, integer ldb, void* q, integer ldq, void* z, integer ldz, integer datatype,
                            integer n_repeats, double* time_min_, integer* info);
void invoke_gghrd(integer datatype, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, void* a,
                            integer* lda, void* b, integer* ldb, void* q, integer* ldq,
                            void* z, integer* ldz, integer* info);
static FILE* g_ext_fptr = NULL;

void fla_test_gghrd(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Reduces a pair matrices (A,B) to generalized upper Hessenberg form";
    char* front_str = "GGHRD";
    integer tests_not_run = 1, invalid_dtype = 0;
    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_gghrd_experiment);
        tests_not_run = 0;
    }
    if(argc == 14)
    {
        /* Read matrix input data from a file */
        g_ext_fptr = fopen(argv[13], "r");
        if (g_ext_fptr == NULL)
        {
            printf("\n Invalid input file argument \n");
            return;
        }
    }
    if(argc >= 13 && argc <=14)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype,type_flag[4] = {0};
        char *endptr;

        /* Prase the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].compq_gghrd = argv[3][0];
        params->lin_solver_paramslist[0].compz_gghrd = argv[4][0];
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ilo = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ihi = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldb = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldq = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldz = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_gghrd_experiment(params, datatype,
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
        printf("\nIllegal arguments for GGHRD\n");
        printf("./<EXE> gghrd <precisions - sdcz> <compq> <compz> <N> <ILO> <IHI> <LDA> <LDB> <LDQ> <LDZ> <repeats>\n");
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

void fla_test_gghrd_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double   *perf,
    double   *time_min,
    double   *residual)
{
    integer n, ldz, lda, ldb, ldq;
    integer ilo, ihi, info = 0, vinfo = 0;
    void *A = NULL, *Z = NULL, *Q = NULL, *B = NULL, *A_test = NULL, *B_test = NULL, *Q_test = NULL, *Z_test = NULL;
    char compz, compq;

    /* Get input matrix dimensions. */
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    ldb = params->lin_solver_paramslist[pci].ldb;
    ldq = params->lin_solver_paramslist[pci].ldq;
    ldz = params->lin_solver_paramslist[pci].ldz;

    if(lda < n || ldq < n || ldz < n || ldb < n)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Initialize parameter */
    compz = params->lin_solver_paramslist[pci].compz_gghrd;
    compq = params->lin_solver_paramslist[pci].compq_gghrd;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    ilo = params->lin_solver_paramslist[pci].ilo;
    ihi = params->lin_solver_paramslist[pci].ihi;

    /* Create input matrix parameters*/
    create_matrix(datatype, &A, lda, n);
    create_matrix(datatype, &B, ldb, n);
    create_matrix(datatype, &Q, ldq, n);
    create_matrix(datatype, &Z, ldz, n);

    if(g_ext_fptr != NULL)
    {
        init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
        init_matrix_from_file(datatype, B, n, n, ldb, g_ext_fptr);
        init_matrix_from_file(datatype, Q, n, n, ldq, g_ext_fptr);
        init_matrix_from_file(datatype, Z, n, n, ldz, g_ext_fptr);
    }
    else
    {
        rand_matrix(datatype, B, n, n, ldb);
        get_orthogonal_matrix_from_QR(datatype, n, B, ldb, Q, ldq, &info);
        if(info < 0)
        {
            *residual = DBL_MAX;
            free_matrix(A);
            free_matrix(B);
            free_matrix(Q);
            free_matrix(Z);
            return;
        }
        if(compq == 'I')
            set_identity_matrix(datatype, n, n, Q, ldq);
        get_generic_triangular_matrix(datatype, n, A, lda, ilo, ihi);
        set_identity_matrix(datatype, n, n, Z, ldz);
    }

    /* Make copy of matrix A,B,Q and Z. This is required to validate the API functionality */
    create_matrix(datatype, &A_test, lda, n);
    create_matrix(datatype, &B_test, ldb, n);
    create_matrix(datatype, &Q_test, ldq, n);
    create_matrix(datatype, &Z_test, ldq, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);
    copy_matrix(datatype, "full", n, n, B, ldb, B_test, ldb);
    copy_matrix(datatype, "full", n, n, Q, ldq, Q_test, ldq);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_test, ldz);

    prepare_gghrd_run(&compq, &compq, n, &ilo, &ihi, A_test, lda, B_test, ldb, Q_test, ldq, Z_test, ldz, datatype, n_repeats, time_min, &info);

    /* Performance computation
       (7)n^3 flops for eigen vectors for real
       (25)n^3 flops for eigen vectors for complex
       (10)n^3 flops for Schur form is computed for real
       (35)n^3 flops for Schur form is computed for complex
       (20)n^3 flops full Schur factorization is computed for real
       (70)n^3 flops full Schur factorization is computed for complex */

    if(compz == 'N')
    {
        if(datatype == FLOAT || datatype == DOUBLE)
            *perf = (double)(7.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
        else
            *perf = (double)(25.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    }
    else if(compz == 'I')
    {
        if(datatype == FLOAT || datatype == DOUBLE)
            *perf = (double)(10.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
        else
            *perf = (double)(35.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    }
    else
    {
        if(datatype == FLOAT || datatype == DOUBLE)
            *perf = (double)(20.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
        else
            *perf = (double)(70.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    }

    /* Output Validation */
    if(info == 0)
        validate_gghrd(&compq, &compz, n, A, A_test, lda, B, B_test, ldb, Q, Q_test, ldq, Z, Z_test, ldz, datatype, residual, &vinfo);
    if(info < 0 || vinfo < 0)
        *residual = DBL_MAX;

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(B);
    free_matrix(Q);
    free_matrix(Z);
    free_matrix(A_test);
    free_matrix(B_test);
    free_matrix(Q_test);
    free_matrix(Z_test);
}

void prepare_gghrd_run(char* compq, char* compz, integer n, integer* ilo, integer* ihi, void* A, integer lda,
                            void* B, integer ldb, void* Q, integer ldq, void* Z, integer ldz, integer datatype,
                            integer n_repeats, double* time_min_, integer* info)
{
    void *A_save = NULL, *B_save = NULL, *Q_save = NULL, *Z_save = NULL;
    integer i;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A,B,Q and Z. Same input values will be passed in each itertaion.*/
    create_matrix(datatype, &A_save, lda, n);
    create_matrix(datatype, &B_save, ldb, n);
    create_matrix(datatype, &Q_save, ldq, n);
    create_matrix(datatype, &Z_save, ldz, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);
    copy_matrix(datatype, "full", n, n, B, ldb, B_save, ldb);
    copy_matrix(datatype, "full", n, n, Q, ldq, Q_save, ldq);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);

    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A,B,Q and Z value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
        copy_matrix(datatype, "full", n, n, B_save, ldb, B, ldb);
        copy_matrix(datatype, "full", n, n, Q_save, ldq, Q, ldq);
        copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);

        exe_time = fla_test_clock();

        /* Call to gghrd API */
        invoke_gghrd(datatype, compq, compz, &n, ilo, ihi, A, &lda, B, &ldb, Q, &ldq, Z, &ldz, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }
    *time_min_ = time_min;

    free_matrix(A_save);
    free_matrix(B_save);
    free_matrix(Q_save);
    free_matrix(Z_save);
}

void invoke_gghrd(integer datatype, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, void* a, integer* lda, void* b, integer* ldb, void* q, integer* ldq, void* z, integer* ldz, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
            break;
        }
    }
}
