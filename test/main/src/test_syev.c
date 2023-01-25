/*
    Copygitright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_syev_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                integer n_repeats, double* perf, double* t, double* residual);
void prepare_syev_run(char* jobz, char* uplo, integer n, void* A, integer lda, void* w, integer datatype, integer n_repeats,
                        double* time_min_, integer* info);
void invoke_syev(integer datatype, char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork,
                     void *rwork, integer* info);

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
static FILE* g_ext_fptr = NULL;

void fla_test_syev(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Eigen Values and Vectors";
    char* front_str = "SYEV";
    integer tests_not_run = 1, invalid_dtype = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_syev_experiment);
        tests_not_run = 0;
    }
    if (argc == 10)
    {
        /* Read matrix input data from a file */
        g_ext_fptr = fopen(argv[9], "r");
        if (g_ext_fptr == NULL)
        {
            printf("\n Invalid input file argument \n");
            return;
        }
    }
    if (argc >= 9 && argc <= 10)
    {
        integer i, num_types,N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype,type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].jobz = argv[3][0];
        params->eig_sym_paramslist[0].uplo = argv[4][0];
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        g_lwork = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_syev_experiment(params, datatype,
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
        printf("\nIllegal arguments for syev\n");
        printf("./<EXE> syev <precisions - sdcz> <JOBZ> <UPLO> <N> <LDA> <LWORK> <repeats>\n");
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

void fla_test_syev_experiment(test_params_t *params,
                               integer  datatype,
                               integer  p_cur,
                               integer  q_cur,
                               integer pci,
                               integer n_repeats,
                               double* perf,
                               double *time_min,
                               double* residual)
{
    integer n, lda, info = 0, vinfo = 0;
    char jobz, uplo;
    void *A = NULL, *w = NULL, *A_test = NULL;

    /* Get input matrix dimensions.*/
    jobz = params->eig_sym_paramslist[pci].jobz;
    uplo = params->eig_sym_paramslist[pci].uplo;
    *residual = params->eig_sym_paramslist[pci].threshold_value;

    n = p_cur;
    lda = params->eig_sym_paramslist[pci].lda;

    if(lda < n)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Create input matrix parameters */
    create_matrix(datatype, &A, lda, n);
    create_realtype_vector(datatype, &w, n);
    if (g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
    }
    else
    {
        /* input matrix A with random symmetric numbers or complex hermitian matrix */
        if (datatype == FLOAT || datatype == DOUBLE)
            rand_sym_matrix(datatype, A, n, n, lda);
        else
            rand_hermitian_matrix(datatype, n, &A, lda);
    }
    /* Make a copy of input matrix A. This is required to validate the API functionality.*/
    create_matrix(datatype, &A_test, lda, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    prepare_syev_run(&jobz, &uplo, n, A_test, lda, w, datatype, n_repeats, time_min, &info);

    /* performance computation
       (8/3)n^3 flops for eigen vectors
       (4/3)n^3 flops for eigen values */
    if( jobz == 'V')
        *perf = (double)((8.0 / 3.0) * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    else
        *perf = (double)((4.0 / 3.0) * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if (info == 0)
        validate_syevd(&jobz, n, A, A_test, lda, w, datatype, residual, &vinfo);

    /* Assigning bigger value to residual as execution fails */
    if (info < 0 || vinfo < 0)
        *residual = DBL_MAX;
        
    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(w);
}

void prepare_syev_run(char *jobz,
                       char *uplo,
                       integer n,
                       void *A,
                       integer lda,
                       void *w,
                       integer datatype,
                       integer n_repeats,
                       double* time_min_,
                       integer* info)
{
    void *A_save = NULL, *work = NULL, *rwork = NULL, *w_test = NULL;
    integer i, lwork;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, lda, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size.*/

    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
        /* call to  syev API */
        invoke_syev(datatype, jobz, uplo, &n, NULL, &lda, NULL, work, &lwork, rwork, info);
        /* Get work size */
        if(*info == 0)
        {
            lwork = get_work_value(datatype, work);
            free_vector(work);
        }
        else
        {
            free_vector(work);
            free_matrix(A_save);
            return;
        }
    }
    else
    {
        lwork = g_lwork;
    }

    for (i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
        create_realtype_vector(datatype, &w_test, n);
        create_vector(datatype, &work, lwork);

        if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX )
                create_realtype_vector(datatype, &rwork, fla_max(1, 3*n -2));
        else
            rwork = NULL;

        exe_time = fla_test_clock();

        /* call to API */
        invoke_syev(datatype, jobz, uplo, &n, A, &lda, w_test, work, &lwork, rwork, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Make a copy of the output buffers. This is required to validate the API functionality.*/
        copy_realtype_vector(datatype, n, w_test, 1, w, 1);

        /* Free up the output buffers */
        free_vector(work);

        if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            free_vector(rwork);

        free_vector(w_test);
    }

    *time_min_ = time_min;

    free_matrix(A_save);
}

void invoke_syev(integer datatype, char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, void* rwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssyev(jobz, uplo, n, a, lda, w, work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dsyev(jobz, uplo, n, a, lda, w, work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
            break;
        }
    }
}
