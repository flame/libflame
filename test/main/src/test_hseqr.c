/*
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"


/* Local prototypes */
void fla_test_hseqr_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, integer einfo, double* perf, double* t, double* residual);
void prepare_hseqr_run(char* job, char* compz, integer n, integer* ilo, integer* ihi, void* h, integer ldh, void *w, void *wr, void* wi,
                            void* z, integer ldz, integer datatype, integer n_repeats, double* time_min_, integer* info);
void invoke_hseqr(integer datatype,char* job, char* compz, integer* n, integer* ilo, integer* ihi, void* h, integer* ldh, void *w,
                    void *wr, void* wi, void* z, integer* ldz, void* work, integer* lwork, integer* info);
static FILE* g_ext_fptr = NULL;

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
void fla_test_hseqr(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Computing Eigen value of a Hessenberg matrix";
    char* front_str = "HSEQR";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_hseqr_experiment);
        tests_not_run = 0;
    }
    if(argc == 13)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[12]);
    }
    if(argc >= 12 && argc <=13)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype,type_flag[4] = {0};
        char *endptr;

        /* Prase the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].job_seqr = argv[3][0];
        params->eig_sym_paramslist[0].compz_hseqr = argv[4][0];
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ilo = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ihi = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ldz = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        g_lwork = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_hseqr_experiment(params, datatype,
                                          N, N,
                                          0,
                                          n_repeats, einfo,
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
        printf("\nIllegal arguments for HSEQR\n");
        printf("./<EXE> hseqr <precisions - sdcz> <job> <compz> <N> <ILO> <IHI> <LDH> <LDZ> <LWORK> <repeats>\n");
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

void fla_test_hseqr_experiment(test_params_t *params,
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
    integer n, ldz, ldh;
    integer ilo, ihi, info = 0, vinfo = 0;
    void *H = NULL, *w = NULL, *wr = NULL, *wi = NULL, *Z = NULL, *H_test = NULL, *Z_Test = NULL;
    void *scale = NULL;
    char compz, job;

    /* Get input matrix dimensions. */
    n = p_cur;
    ldz = params->eig_sym_paramslist[pci].ldz;
    ldh = params->eig_sym_paramslist[pci].lda;

    if(ldz < n || ldh < n)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Initialize parameter needed for HSEQR() call. */
    job = params->eig_sym_paramslist[pci].job_seqr;
    compz = params->eig_sym_paramslist[pci].compz_hseqr;
    *residual = params->eig_sym_paramslist[pci].threshold_value;
    ilo = params->eig_sym_paramslist[pci].ilo;
    ihi = params->eig_sym_paramslist[pci].ihi;

    /* Create input matrix parameters*/
    create_matrix(datatype, &H, ldh, n);
    create_matrix(datatype, &Z, ldz, n);
    create_vector(datatype, &scale, n);

    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_vector(datatype, &w, n);
    }
    else
    {
        create_vector(datatype, &wr, n);
        create_vector(datatype, &wi, n);
    }

    if(g_ext_fptr != NULL)
    {
        init_matrix_from_file(datatype, H, n, n, ldh, g_ext_fptr);
        init_matrix_from_file(datatype, Z, n, n, ldz, g_ext_fptr);
    }
    else
    {
        /* Generate Hessenberg matrix H */
        get_hessenberg_matrix(datatype, n, H, ldh, Z, ldz, &ilo, &ihi, scale, &info);
        if(compz == 'I')
            set_identity_matrix(datatype, n, n, Z, ldz);
    }

    /* Make copy of matrix H and Z. This is required to validate the API functionality */
    create_matrix(datatype, &H_test, ldh, n);
    create_matrix(datatype, &Z_Test, ldz, n);
    copy_matrix(datatype, "full", n, n, H, ldh, H_test, ldh);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_Test, ldz);

    prepare_hseqr_run(&job, &compz, n, &ilo, &ihi, H_test, ldh, w, wr, wi, Z_Test, ldz, datatype, n_repeats, time_min, &info);

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
        validate_hseqr(&job, &compz, n, H, H_test, ldh, Z, Z_Test, ldz, datatype, residual, &vinfo);
    
    FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_vector(scale);
    free_matrix(H);
    free_matrix(Z);
    free_matrix(H_test);
    free_matrix(Z_Test);
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        free_vector(w);
    }
    else
    {
        free_vector(wr);
        free_vector(wi);
    }
}

void prepare_hseqr_run(char* job,
    char* compz,
    integer n,
    integer* ilo,
    integer* ihi,
    void* H,
    integer ldh,
    void *w,
    void *wr,
    void* wi,
    void* Z,
    integer ldz,
    integer datatype,
    integer n_repeats,
    double* time_min_,
    integer* info)
{
    void *H_save = NULL, *work = NULL, *Z_save = NULL;
    integer i, lwork;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix H and Z. Same input values will be passed in each itertaion.*/
    create_matrix(datatype, &H_save, ldh, n);
    create_matrix(datatype, &Z_save, ldz, n);
    copy_matrix(datatype, "full", n, n, H, ldh, H_save, ldh);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.*/
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* call to  hseqr API */
        invoke_hseqr(datatype, job, compz, &n, ilo, ihi, NULL, &ldh,
                    NULL, NULL, NULL, NULL, &ldz, work, &lwork, info);

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
        copy_matrix(datatype, "full", n, n, H_save, ldh, H, ldh);
        copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);
        create_vector(datatype, &work, lwork);

        exe_time = fla_test_clock();

        /* Call to hseqr API */
        invoke_hseqr(datatype, job, compz, &n, ilo, ihi, H, &ldh, w, wr, wi, Z, &ldz, work, &lwork, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
    }
    *time_min_ = time_min;

    free(H_save);
    free(Z_save);
}

void invoke_hseqr(integer datatype,char* job, char* compz, integer* n, integer* ilo, integer* ihi, void* h, integer* ldh, void *w, void *wr, void* wi, void* z, integer* ldz, void* work, integer* lwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_shseqr(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dhseqr(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_chseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zhseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info);
            break;
        }
    }
}
