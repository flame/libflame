/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"


/* Local prototypes */
void fla_test_ggev_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci, integer n_repeats, double* perf, double* t, double* residual);
void prepare_ggev_run(char *jobvl, char *jobvr, integer n, void *a, integer lda, void *b, integer ldb, void* alpha, void * alphar, void * alphai, void *beta,	void *vl, integer ldvl, 
                      void *vr, integer ldvr,	integer datatype, integer n_repeats, double* time_min_);
void invoke_ggev(integer datatype, char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* b, integer* ldb, integer* alpha, integer* alphar,
    integer* alphai, integer* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, void* rwork, integer* info);

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
void fla_test_ggev(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Computing Eigen value and Eigen vectors";
    char* front_str = "GGEV";
    integer tests_not_run = 1, invalid_dtype = 0;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_NSYM, fla_test_ggev_experiment);
        tests_not_run = 0;
    }
    else if(argc == 12)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_non_sym_paramslist[0].jobvsl = argv[3][0];
        params->eig_non_sym_paramslist[0].jobvsr = argv[4][0];
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->eig_non_sym_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->eig_non_sym_paramslist[0].ldb = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->eig_non_sym_paramslist[0].ldvl = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        params->eig_non_sym_paramslist[0].ldvr = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);

        g_lwork = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->eig_non_sym_paramslist[0].GenNonSymEigProblem_threshold = CLI_NORM_THRESH;

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
                fla_test_ggev_experiment(params, datatype,
                                          N, N,
                                          0,
                                          n_repeats,
                                          &perf, &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str,
                                      stype,
                                      SQUARE_INPUT,
                                      N, N,
                                      residual, params->eig_non_sym_paramslist[0].GenNonSymEigProblem_threshold,
                                      time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for ggev\n");
        printf("./<EXE> ggev <precisions - sdcz> <jobvl> <jobvr> <N> <LDA> <LDB> <LDVL> <LDVR> <LWORK> <repeats>\n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
    return;
}


void fla_test_ggev_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double   *perf,
    double   *t,
    double   *residual)
{
    integer m, lda, ldvl, ldvr, ldb;
    void *A = NULL, *B = NULL, *VL = NULL, *VR = NULL;
    void * alpha = NULL, *alphar=NULL, *alphai=NULL, *beta, *A_test , *B_test;
    double time_min = 1e9;
    *residual = params->eig_non_sym_paramslist[pci].GenNonSymEigProblem_threshold;
    //char JOBVL = params->eig_non_sym_paramslist[pci].jobvsl;
    //char JOBVR = params->eig_non_sym_paramslist[pci].jobvsr;
    char JOBVL = 'V';
    char JOBVR = 'V';

    /* Get input matrix dimensions */
    m = p_cur;

    lda = params->eig_non_sym_paramslist[pci].lda;
    ldb = params->eig_non_sym_paramslist[pci].ldb;
    ldvl = params->eig_non_sym_paramslist[pci].ldvl;
    ldvr = params->eig_non_sym_paramslist[pci].ldvr;

    /* Create input matrix parameters */
    create_matrix(datatype, &A, lda, m);
    create_matrix(datatype, &B, ldb, m);
    create_matrix(datatype, &VL, ldvl, m);
    create_matrix(datatype, &VR, ldvr, m);
    if (datatype == FLOAT || datatype == DOUBLE)
    {
        create_vector(datatype, &alphar, m);
        create_vector(datatype, &alphai, m);
    }
    else
    {
        create_vector(datatype, &alpha, m);
    }
    create_vector(datatype, &beta, m);

    /* Initialize input matrix A with random numbers */
    rand_matrix(datatype, A, m, m, lda);
    rand_matrix(datatype, B, m, m, ldb);

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, &A_test, lda, m);
    create_matrix(datatype, &B_test, ldb, m);
    copy_matrix(datatype, "full", m, m, A, lda, A_test, lda);
    copy_matrix(datatype, "full", m, m, B, ldb, B_test, ldb);

    prepare_ggev_run(&JOBVL,&JOBVR, m, A_test, lda, B_test, ldb, alpha, alphar, alphai, beta,  VL, ldvl, VR, ldvr, datatype, n_repeats, &time_min);

    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2m^3 - (2/3)m^3 flops */
    *perf = (double)((2.0 * m * m * m) - ((2.0 / 3.0) * m * m * m)) / time_min / FLOPS_PER_UNIT_PERF;
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if (JOBVL == 'V' && JOBVR == 'V')
    {
        validate_ggev(&JOBVL, &JOBVR, m, A, lda, B, ldb, alpha, alphar, alphai, beta, VL, ldvl, VR, ldvr, datatype, residual);
    }
    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(VL);
    free_matrix(VR);
    free_matrix(B);
    free_matrix(B_test);
    if (datatype == FLOAT || datatype == DOUBLE)
    {
        free_vector(alphar);
        free_vector(alphai);
    }
    else
    {
        free_vector(alpha);
    }

    free_vector(beta);

}

void prepare_ggev_run(char *jobvl, char *jobvr, integer n_A, void *A, integer lda,
                        void *B, integer ldb, void* alpha, void * alphar, void * alphai, void* beta,
                        void *VL, integer ldvl, void *VR, integer ldvr,
                        integer datatype, integer n_repeats, double* time_min_)
{
    void *A_save = NULL, *B_save = NULL , *work = NULL, * rwork = NULL;
    integer i;
    integer lwork, info = 0;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, lda, n_A);
    copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);
    create_matrix(datatype, &B_save, ldb, n_A);
    copy_matrix(datatype, "full", n_A, n_A, B, ldb, B_save, ldb);

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size. */
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 8*n_A);

        /* call to  ggev API to get work query */
        invoke_ggev(datatype, jobvl, jobvr, &n_A, NULL, &lda, NULL, &ldb, NULL, NULL, NULL, NULL, NULL, &ldvl, NULL, &ldvr, work, &lwork, rwork, &info);

        /* Get work size */
        lwork = get_work_value( datatype, work);

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
        /* Restore input matrix A value and allocate memory to output buffers for each iteration */
        copy_matrix(datatype, "full", n_A, n_A, A_save, lda, A, lda);
        copy_matrix(datatype, "full", n_A, n_A, B_save, ldb, B, ldb);
        create_vector(datatype, &work, lwork);
        if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            create_realtype_vector(datatype, &rwork, 8*n_A);
        }
        exe_time = fla_test_clock();

        /* call to API */

        invoke_ggev(datatype, jobvl, jobvr, &n_A, A, &lda, B, &ldb, alpha, alphar, alphai, beta, VL, &ldvl, VR, &ldvr, work, &lwork, rwork , &info);
        
        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
        if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            free_vector(rwork);
        }
    }

    *time_min_ = time_min;
    copy_matrix(datatype, "full", n_A, n_A, A_save, lda, A, lda);
    copy_matrix(datatype, "full", n_A, n_A, B_save, ldb, B, ldb);

    free(A_save);
    free(B_save);


}


void invoke_ggev(integer datatype, char *jobvl, char *jobvr,integer *n, void *a, integer *lda, void *b, integer *ldb, integer* alpha, integer * alphar,
    integer* alphai, integer* beta, void *vl, integer *ldvl, void *vr, integer *ldvr, void* work, integer* lwork, void* rwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sggev_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
            break;
        }
        
        case DOUBLE:
        {
            dggev_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            cggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
            break;
        }
    }
}
