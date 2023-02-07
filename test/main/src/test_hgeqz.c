/*
    Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"


/* Local prototypes */
void fla_test_hgeqz_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_hgeqz_run(char* job, char* compq, char* compz, integer n, integer* ilo, integer* ihi, void* h, integer ldh,
                            void* t, integer ldt, void *alpha, void *alphar, void *alphai, void* beta, void* q, integer ldq, void* z,
                            integer ldz, integer datatype, integer n_repeats, double* time_min_, integer* info);
void invoke_hgeqz(integer datatype, char* job, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, void* h,
                            integer* ldh, void* t, integer* ldt, void *alpha, void *alphar, void *alphai, void* beta, void* q, integer* ldq,
                            void* z, integer* ldz, void* work, integer* lwork, void* rwork, integer* info);
static FILE* g_ext_fptr = NULL;

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
void fla_test_hgeqz(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Computing Eigen value of a real matrix pair (H,T)";
    char* front_str = "HGEQZ";
    integer tests_not_run = 1, invalid_dtype = 0;
    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_hgeqz_experiment);
        tests_not_run = 0;
    }
    if(argc == 16)
    {
        /* Read matrix input data from a file */
        g_ext_fptr = fopen(argv[15], "r");
        if (g_ext_fptr == NULL)
        {
            printf("\n Invalid input file argument \n");
            return;
        }
    }
    if(argc >= 15 && argc <=16)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype,type_flag[4] = {0};
        char *endptr;

        /* Prase the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].job_seqr = argv[3][0];
        params->eig_sym_paramslist[0].compq_hgeqz = argv[4][0];
        params->eig_sym_paramslist[0].compz_hgeqz = argv[5][0];
        N = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ilo = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ihi = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].lda = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ldb = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ldq = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ldz = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);
        g_lwork = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[14], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_hgeqz_experiment(params, datatype,
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
        printf("\nIllegal arguments for HGEQZ\n");
        printf("./<EXE> hgeqz <precisions - sdcz> <job> <compq> <compz> <N> <ILO> <IHI> <LDH> <LDT> <LDQ> <LDZ> <LWORK> <repeats>\n");
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

void fla_test_hgeqz_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double   *perf,
    double   *time_min,
    double   *residual)
{
    integer n, ldz, ldh, ldt, ldq;
    integer ilo, ihi, info = 0, vinfo = 0;
    void *H = NULL, *Z = NULL, *Q = NULL, *T = NULL, *H_test = NULL, *T_test = NULL, *A = NULL;
    void *B = NULL, *Z_A = NULL, *Q_test = NULL, *Z_test = NULL;
    void *alpha = NULL, *alphar = NULL, *alphai = NULL, *beta = NULL, *scale = NULL, *Q_A = NULL;
    char compz, compq, job;

    /* Get input matrix dimensions. */
    n = p_cur;
    ldh = params->eig_sym_paramslist[pci].lda;
    ldt = params->eig_sym_paramslist[pci].ldb;
    ldq = params->eig_sym_paramslist[pci].ldq;
    ldz = params->eig_sym_paramslist[pci].ldz;

    if(ldh < n || ldq < n || ldz < n || ldt < n)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Initialize parameters */
    job = params->eig_sym_paramslist[pci].job_seqr;
    compz = params->eig_sym_paramslist[pci].compz_hgeqz;
    compq = params->eig_sym_paramslist[pci].compq_hgeqz;
    *residual = params->eig_sym_paramslist[pci].threshold_value;
    ilo = params->eig_sym_paramslist[pci].ilo;
    ihi = params->eig_sym_paramslist[pci].ihi;

    /* Create input matrix */
    create_matrix(datatype, &H, ldh, n);
    create_matrix(datatype, &T, ldt, n);
    create_matrix(datatype, &Q, ldq, n);
    create_matrix(datatype, &Z, ldz, n);
    create_matrix(datatype, &A, ldh, n);
    create_matrix(datatype, &B, ldt, n);
    create_matrix(datatype, &Q_A, ldq, n);
    create_matrix(datatype, &Z_A, ldz, n);
    if (datatype == FLOAT || datatype == DOUBLE)
    {
        create_vector(datatype, &alphar, n);
        create_vector(datatype, &alphai, n);
    }
    else
    {
        create_vector(datatype, &alpha, n);
    }
    create_vector(datatype, &beta, n);

    if(g_ext_fptr != NULL)
    {
        init_matrix_from_file(datatype, H, n, n, ldh, g_ext_fptr);
        init_matrix_from_file(datatype, T, n, n, ldt, g_ext_fptr);
        init_matrix_from_file(datatype, Q, n, n, ldq, g_ext_fptr);
        init_matrix_from_file(datatype, Z, n, n, ldz, g_ext_fptr);
    }
    else
    {
        /* Convert matrix according to ILO and IHI values */
        get_generic_triangular_matrix(datatype, n, A, ldh, ilo, ihi);
        /* Initialize matrix with random values */
        rand_matrix(datatype, B, n, n, ldt);
        /* Decompose matrix B in to QR and store orthogonal matrix in Q and R in B */
        get_orthogonal_matrix_from_QR(datatype, n, B, ldt, Q, ldq, &info);
        /* Make copy of matrix A and B. This is required to validate the API functionality */
        copy_matrix(datatype, "full", n, n, A, ldh, H, ldh);
        copy_matrix(datatype, "full", n, n, B, ldt, T, ldt);
        if(compq == 'I')
            set_identity_matrix(datatype, n, n, Q, ldq);
        set_identity_matrix(datatype, n, n, Z, ldz);
        /* Make copy of matrix Q and Z. This is required to validate the API functionality */
        copy_matrix(datatype, "full", n, n, Q, ldq, Q_A, ldq);
        copy_matrix(datatype, "full", n, n, Z, ldz, Z_A, ldz);
        /* Call to GGHRD API */
        invoke_gghrd(datatype, &compq, &compz, &n, &ilo, &ihi, H, &ldh, T, &ldt, Q, &ldq, Z, &ldz, &info);
        if(info < 0)
        {
            *residual = DBL_MAX;
            free_matrix(H);
            free_matrix(T);
            free_matrix(Q);
            free_matrix(Z);
            free_matrix(A);
            free_matrix(B);
            free_matrix(Q_A);
            free_matrix(Z_A);
            if (datatype == FLOAT || datatype == DOUBLE)
            {
                free_vector(alphar);
                free_vector(alphai);
            }
            else
                free_vector(alpha);
            free_vector(beta);
            return;
        }
        if(compq == 'I')
            set_identity_matrix(datatype, n, n, Q, ldq);
        if(compz == 'I')
            set_identity_matrix(datatype, n, n, Z, ldz);
    }

    /* Make copy of matrix H,T,Q and Z. This is required to validate the API functionality */
    create_matrix(datatype, &H_test, ldh, n);
    create_matrix(datatype, &T_test, ldt, n);
    create_matrix(datatype, &Q_test, ldq, n);
    create_matrix(datatype, &Z_test, ldz, n);
    copy_matrix(datatype, "full", n, n, H, ldh, H_test, ldh);
    copy_matrix(datatype, "full", n, n, T, ldt, T_test, ldt);
    copy_matrix(datatype, "full", n, n, Q, ldq, Q_test, ldq);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_test, ldz);

    prepare_hgeqz_run(&job, &compq, &compq, n, &ilo, &ihi, H_test, ldh, T_test, ldt, alpha, alphar, alphai, beta, Q_test, ldq, Z_test,
                        ldz, datatype, n_repeats, time_min, &info);

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
        validate_hgeqz(&job, &compq, &compz, n, H, H_test, A, ldh, T, T_test, B, ldt, Q, Q_test, Q_A, ldq, Z, Z_test, Z_A, ldz, datatype, residual, &vinfo);
    if(info < 0 || vinfo < 0)
        *residual = DBL_MAX;

    /* Free up the buffers */
    free_vector(scale);
    free_vector(beta);
    free_matrix(H);
    free_matrix(T);
    free_matrix(Q);
    free_matrix(Z);
    free_matrix(H_test);
    free_matrix(T_test);
    if (datatype == FLOAT || datatype == DOUBLE)
    {
        free_vector(alphar);
        free_vector(alphai);
    }
    else
        free_vector(alpha);
}

void prepare_hgeqz_run(char* job, char* compq, char* compz, integer n, integer* ilo, integer* ihi, void* H, integer ldh,
                            void* T, integer ldt, void *alpha, void *alphar, void *alphai, void* beta, void* Q, integer ldq, void* Z,
                            integer ldz, integer datatype, integer n_repeats, double* time_min_, integer* info)
{
    void *H_save = NULL, *T_save = NULL, *Q_save = NULL, *Z_save = NULL, *work = NULL, *rwork;
    integer i, lwork;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix H,T,Q and Z. Same input values will be passed in each itertaion.*/
    create_matrix(datatype, &H_save, ldh, n);
    create_matrix(datatype, &T_save, ldt, n);
    create_matrix(datatype, &Q_save, ldq, n);
    create_matrix(datatype, &Z_save, ldz, n);
    create_realtype_vector(datatype, &rwork, n);
    copy_matrix(datatype, "full", n, n, H, ldh, H_save, ldh);
    copy_matrix(datatype, "full", n, n, T, ldt, T_save, ldt);
    copy_matrix(datatype, "full", n, n, Q, ldq, Q_save, ldq);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.*/
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* call to  hgeqz API */
        invoke_hgeqz(datatype, job, compq, compz, &n, ilo, ihi, NULL, &ldh, NULL, &ldt, alpha, alphar, alphai,
                         beta, NULL, &ldq, NULL, &ldz, work, &lwork, rwork, info);

        /* Output buffers will be freshly allocated for each iterations, free up
        the current output buffers.*/
        if(*info == 0)
        {
            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
        }
        else
        {
            free_vector(work);
            free_matrix(H_save);
            free_matrix(T_save);
            free_matrix(Q_save);
            free_matrix(Z_save);
            return;
        }
    }
    else
    {
        lwork = g_lwork;
    }

    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix H and Z value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, H_save, ldh, H, ldh);
        copy_matrix(datatype, "full", n, n, T_save, ldt, T, ldt);
        copy_matrix(datatype, "full", n, n, Q_save, ldq, Q, ldq);
        copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);
        create_vector(datatype, &work, lwork);

        exe_time = fla_test_clock();

        /* Call to hgeqz API */
        invoke_hgeqz(datatype, job, compq, compz, &n, ilo, ihi, H, &ldh, T, &ldt, alpha, alphar, alphai,
                        beta, Q, &ldq, Z, &ldz, work, &lwork, rwork, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
    }
    *time_min_ = time_min;

    free(H_save);
    free(Z_save);
    free(Q_save);
    free(T_save);
    free_vector(rwork);
}

void invoke_hgeqz(integer datatype, char* job, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, void* h, integer* ldh, void* t, integer* ldt, void *alpha, void *alphar, void *alphai, void* beta, void* q, integer* ldq, void* z, integer* ldz, void* work, integer* lwork, void* rwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_shgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dhgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_chgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zhgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, info);
            break;
        }
    }
}
