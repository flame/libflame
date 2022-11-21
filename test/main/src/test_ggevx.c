/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"


/* Local prototypes */
void fla_test_ggevx_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_ggevx_run(char* balanc, char* jobvl, char* jobvr, char* sense, integer n_A, 
                        void* A, integer lda, void* B, integer ldb, void* alpha, void* alphar, void* alphai, void* beta,
                        void* VL, integer ldvl, void* VR, integer ldvr,
                        integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, 
                        void* rcondv, integer datatype, integer n_repeats, double* time_min_);
void invoke_ggevx(integer datatype, char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, 
                        void* b, integer* ldb, void* alpha, void* alphar, void* alphai, void* beta, void* vl, integer* ldvl, 
                        void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, 
                        void* rconde, void* rcondv, void* work, integer* lwork, void* rwork, integer* iwork, integer* bwork, integer* info);

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
FILE* g_ext_fptr = NULL;

void fla_test_ggevx(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Computing Eigen value and Eigen vectors with condition numbers";
    char* front_str = "GGEVX";
    integer tests_not_run = 1, invalid_dtype = 0;

    if(argc == 1)
    {
        /* Test with parameters from config */
        g_lwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_NSYM, fla_test_ggevx_experiment);
        tests_not_run = 0;
    }
    if(argc == 15)
    {
        /* Read matrix input data from a file */
        g_ext_fptr = fopen(argv[14], "r");
        if (g_ext_fptr == NULL)
        {
            printf("\n Invalid input file argument \n");
            return;
        }
    }
    if(argc >= 14 && argc <= 15)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_non_sym_paramslist[0].balance_ggevx = argv[3][0];
        params->eig_non_sym_paramslist[0].jobvsl = argv[4][0];
        params->eig_non_sym_paramslist[0].jobvsr = argv[5][0];
        params->eig_non_sym_paramslist[0].sense_ggevx = argv[6][0];
        N = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->eig_non_sym_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        params->eig_non_sym_paramslist[0].ldb = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        params->eig_non_sym_paramslist[0].ldvl = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        params->eig_non_sym_paramslist[0].ldvr = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
        g_lwork = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_ggevx_experiment(params, datatype,
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
        printf("\nIllegal arguments for GGEVX\n");
        printf("./<EXE> ggevx <precisions - sdcz> <balanc> <jobvl> <jobvr> <sense> <N> <LDA> <LDB> <LDVL> <LDVR> <LWORK> <repeats>\n");
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


void fla_test_ggevx_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer  pci,
    integer  n_repeats,
    double   *perf,
    double   *t,
    double   *residual)
{
    integer n, lda, ldb, ldvl, ldvr;
    integer ilo, ihi;
    void *A = NULL, *B = NULL, *VL = NULL, *VR = NULL;
    void *rscale=NULL, *lscale=NULL, *alpha = NULL, *alphar=NULL, *alphai=NULL, *beta=NULL, *A_test=NULL , *B_test=NULL;
    void *abnrm = NULL, *bbnrm = NULL, *rconde=NULL, *rcondv=NULL;
    double time_min = 1e9;
    char JOBVL = params->eig_non_sym_paramslist[pci].jobvsl;
    char JOBVR = params->eig_non_sym_paramslist[pci].jobvsr;
    char BALANC = params->eig_non_sym_paramslist[pci].balance_ggevx;
    char SENSE = params->eig_non_sym_paramslist[pci].sense_ggevx;
    *residual = params->eig_non_sym_paramslist[pci].GenNonSymEigProblem_threshold;
    /* Get input matrix dimensions */
    n = p_cur;
    lda = n;
    ldb = n;
    ldvl = n;
    ldvr = n;

    /* Create input matrix parameters */
    create_matrix(datatype, &A, lda, n);
    create_matrix(datatype, &B, ldb, n);
    create_matrix(datatype, &VL, ldvl, n);
    create_matrix(datatype, &VR, ldvr, n);
    create_vector(datatype, &beta, n);

    if(datatype == FLOAT || datatype == DOUBLE)
    {
        create_vector(datatype, &alphar, n);
        create_vector(datatype, &alphai, n);
    }
    else
    {
        create_vector(datatype, &alpha, n);
    }

    create_realtype_vector(datatype, &lscale, n);
    create_realtype_vector(datatype, &rscale, n);
    create_realtype_vector(datatype, &rconde, n);
    create_realtype_vector(datatype, &rcondv, n);
    create_realtype_vector(datatype, &abnrm, 1);
    create_realtype_vector(datatype, &bbnrm, 1);


    if(g_ext_fptr != NULL)
    {
        init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
        init_matrix_from_file(datatype, B, n, n, lda, g_ext_fptr);
    }
    else
    {
        /* Initialize input matrix A with random numbers */
        rand_matrix(datatype, A, n, n, lda);
        rand_matrix(datatype, B, n, n, ldb);
    }
    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, &A_test, n, n);
    create_matrix(datatype, &B_test, n, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, n);
    copy_matrix(datatype, "full", n, n, B, ldb, B_test, n);

    prepare_ggevx_run(&BALANC, &JOBVL, &JOBVR, &SENSE, n, A, lda, B, ldb, alpha, alphar, alphai, beta, VL, ldvl, VR, ldvr,
                        &ilo, &ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, datatype, n_repeats, &time_min);

    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2mn^2 - (2/3)n^3 flops */
    *perf = (double)((2.0 * n * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if(JOBVL == 'V' || JOBVR == 'V')
    {
        validate_ggevx(&BALANC, &JOBVL, &JOBVR, &SENSE, n, A_test, n, B_test, n, alpha, alphar, alphai, beta, VL, ldvl, VR, ldvr, datatype, residual);
    }

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(B);
    free_matrix(B_test);
    free_matrix(VL);
    free_matrix(VR);
    if(datatype == FLOAT || datatype == DOUBLE)
    {
        free_vector(alphar);
        free_vector(alphai);
    }
    else
    {
        free_vector(alpha);
    }

    free_vector(beta);
    free_vector(lscale);
    free_vector(rscale);
    free_vector(rconde);
    free_vector(rcondv);
    free_vector(abnrm);
    free_vector(bbnrm);

}

void prepare_ggevx_run(char* balanc, char* jobvl, char* jobvr, char* sense, integer n_A, 
                        void* A, integer lda, void* B, integer ldb, void* alpha, void* alphar, void* alphai, void* beta,
                        void* VL, integer ldvl, void* VR, integer ldvr,
                        integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, 
                        void* rcondv, integer datatype, integer n_repeats, double* time_min_)
{
    void *A_save = NULL, *B_save = NULL , *work = NULL, *rwork = NULL, *iwork = NULL, *bwork=NULL;;
    integer i;
    integer lwork, info = 0;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, lda, n_A);
    create_matrix(datatype, &B_save, ldb, n_A);
    copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);
    copy_matrix(datatype, "full", n_A, n_A, B, ldb, B_save, ldb);

    if(*sense != 'E')
    {
        create_vector(INTEGER, &iwork, 6 + n_A);
    }
    if(*sense != 'N')
    {
        create_vector(INTEGER, &bwork, n_A);
    }

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size. */
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* call to  ggevx API */
        invoke_ggevx(datatype, balanc, jobvl, jobvr, sense, &n_A, A, &lda, B, &ldb, alpha, alphar, alphai, beta, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork, &info);

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

    create_realtype_vector(datatype, &rwork, 8 * n_A);


    for(i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers for each iteration */
        copy_matrix(datatype, "full", n_A, n_A, A_save, n_A, A, lda);
        copy_matrix(datatype, "full", n_A, n_A, B_save, n_A, B, ldb);
        create_vector(datatype, &work, lwork);

        exe_time = fla_test_clock();

        /* call to API */

        invoke_ggevx(datatype, balanc, jobvl, jobvr, sense, &n_A, A, &lda, B, &ldb, alpha, alphar, alphai, beta, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork, &info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
    }

    *time_min_ = time_min;

    free_matrix(A_save);
    free_matrix(B_save);
    if(*sense != 'E')
    {
        free_vector(iwork);
    }
    if(*sense != 'N')
    {
        free_vector(bwork);
    }
    free_vector(rwork);

}


void invoke_ggevx(integer datatype, char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alpha, void* alphar, void* alphai, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, void* rcondv, void* work, integer* lwork, void *rwork, integer* iwork, integer* bwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork, bwork, info);
            break;
        }

        case DOUBLE:
        {
            dggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork, bwork, info);
            break;
        }

        case COMPLEX:
        {
            cggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, rwork, iwork, bwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, rwork, iwork, bwork, info);
            break;
        }
    }
}
