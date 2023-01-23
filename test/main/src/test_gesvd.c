/*
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes */
void fla_test_gesvd_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci, integer n_repeats, double* perf, double* t, double* residual);
void prepare_gesvd_run(char *jobu, char *jobvt, integer m_A, integer n_A, void *A, integer lda, void *s, void *U, integer ldu, void *V, integer ldvt, integer datatype, integer n_repeats, double* time_min_, integer *info);
void invoke_gesvd(integer datatype, char* jobu, char *jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* iwork, integer* info);

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
static FILE* g_ext_fptr = NULL;

void fla_test_gesvd(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Singular value decomposition";
    char* front_str = "GESVD";
    integer tests_not_run = 1, invalid_dtype = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, SVD, fla_test_gesvd_experiment);
        tests_not_run = 0;
    }
    if (argc == 13)
    {
        /* Read matrix input data from a file */
        g_ext_fptr = fopen(argv[12], "r");
        if (g_ext_fptr == NULL)
        {
            printf("\n Invalid input file argument \n");
            return;
        }
    }
    if (argc >= 12 && argc <= 13)
    {
        integer i, num_types,N,M;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype,type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->svd_paramslist[0].jobu_gesvd = argv[3][0];
        params->svd_paramslist[0].jobvt_gesvd = argv[4][0];
        M = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->svd_paramslist[0].lda = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->svd_paramslist[0].ldu = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        params->svd_paramslist[0].ldvt = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);

        g_lwork = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        
        n_repeats = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->svd_paramslist[0].svd_threshold = CLI_NORM_THRESH;

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
                fla_test_gesvd_experiment(params, datatype,
                                          M, N,
                                          0,
                                          n_repeats,
                                          &perf, &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str,
                                      stype,
                                      RECT_INPUT,
                                      M, N,
                                      residual, params->svd_paramslist[0].svd_threshold,
                                      time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gesdd\n");
        printf("./<EXE> gesvd<precisions - sdcz> <JOBU> <JOBVT> <M> <N> <LDA> <LDU> <LDVT> <LWORK> <repeats>\n");
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


void fla_test_gesvd_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    double* perf,
    double *time_min,
    double* residual)
{
    integer m, n, lda, ldu, ldvt;
    integer info = 0, vinfo = 0;
    char jobu, jobvt;
    void *A = NULL, *U = NULL, *V = NULL, *s = NULL, *A_test = NULL;
    *residual = params->svd_paramslist[pci].svd_threshold;

    /* Get input matrix dimensions. */
    jobu = params->svd_paramslist[pci].jobu_gesvd;
    jobvt = params->svd_paramslist[pci].jobvt_gesvd;

    if(jobu == 'S' || jobvt == 'S' || (jobu == 'O' && jobvt != 'O') || (jobu != 'O' && jobvt == 'O'))
    {
           m = p_cur;
           n = p_cur;
    }
    else
    {
           m = p_cur;
           n = q_cur;
    }

    lda = params->svd_paramslist[pci].lda;
    ldu = params->svd_paramslist[pci].ldu;
    ldvt = params->svd_paramslist[pci].ldvt;

    if(lda < m || ldu < m || ldvt < n)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Create input matrix parameters. */
    create_matrix(datatype, &A, lda, n);
    create_matrix(datatype, &U, ldu, m);
    create_matrix(datatype, &V, ldvt, n);
    create_realtype_vector(datatype, &s, fla_min(m, n));

    if (g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, m, n, lda, g_ext_fptr);
    }
    else
    {
        /* Initialize input matrix with random numbers */
        rand_matrix(datatype, A, m, n, lda);
    }

    /* Make a copy of input matrix A. This is required to validate the API functionality. */
    create_matrix(datatype, &A_test, lda, n);
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    prepare_gesvd_run(&jobu, &jobvt, m, n, A_test, lda, s, U, ldu, V, ldvt, datatype, n_repeats, time_min, &info);

    /* Performance Computation
     * Singular values only, 4mn^2 - 4n^3/3 flops
     * Singular values and some singular vectors U (m x n) and V (n x n), 14mn^2 + 8n^3 flops 
     * Link : http://icl.cs.utk.edu/magma/forum/viewtopic.php?f=2&t=921 */
    if(jobu == 'N' || jobvt == 'N')
    {
           if(m >= n)
              *perf = (double)((4.0 * m * n * n) - (( 4.0 * n * n * n ) / 3.0 )) / *time_min / FLOPS_PER_UNIT_PERF;
           else
              *perf = (double)((4.0 * n * m * m) - ((4.0 * m * m * m ) / 3.0)) / *time_min / FLOPS_PER_UNIT_PERF;
    }
    else
    {
       if(m >= n)
              *perf = (double)((14.0 * m * n * n) + ( 8.0  * n * n * n )) / *time_min / FLOPS_PER_UNIT_PERF;
           else
              *perf = (double)((14.0 * n * m * m) + (8.0 * m * m * m )) / *time_min / FLOPS_PER_UNIT_PERF;
    }
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
              *perf *= 4.0;

    /* output validation */
    if((jobu == 'A' && jobvt == 'A') && info == 0)
           validate_gesvd(&jobu, &jobvt, m, n, A, A_test, lda, s, U, ldu, V, ldvt, datatype, residual, &vinfo); 
    
    /* Assigning bigger value to residual as execution fails*/
    if(info < 0 || vinfo < 0)
        *residual = DBL_MAX;

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(U);
    free_matrix(V);
    free_vector(s);
}



void prepare_gesvd_run(char *jobu, char *jobvt,
    integer m_A, integer n_A,
    void *A,
    integer lda,
    void *s,
    void *U,
    integer ldu,
    void *V,
    integer ldvt,
    integer datatype,
    integer n_repeats,
    double* time_min_,
    integer* info)
{
    integer min_m_n, max_m_n;
    void *A_save, *s_test;
    void *work, *iwork, *rwork;
    void *U_test, *V_test;
    integer lwork, liwork, lrwork;
    integer i;
    double time_min = 1e9, exe_time;

    min_m_n = fla_min(m_A, n_A);
    max_m_n = fla_max(m_A, n_A);

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, lda, n_A);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    /* Get rwork and iwork array size since it is not depedent on internal blocks*/
    lrwork = fla_max( (5 * min_m_n * min_m_n + 5 * min_m_n) , ( 2 * max_m_n * min_m_n + 2 * min_m_n * min_m_n + min_m_n));
    liwork = 8 * min_m_n;

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size. */
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* call gesvd API */
        invoke_gesvd(datatype, jobu, jobvt, &m_A, &n_A, NULL, &lda, NULL, NULL, &ldu, NULL, &ldvt, work, &lwork, NULL, NULL, info);
        if(*info < 0)
        {
            free_matrix(A_save);
            free_vector(work);
            return;
        }

        /* Get the work size */
        lwork = get_work_value( datatype, work );

        /* Output buffers will be freshly allocated for each iterations, free up 
        the current output buffers.*/ 
        free_vector(work);
    }
    else
    {
         lwork = g_lwork;
    }

    for (i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);

        create_matrix(datatype, &U_test, ldu, m_A);
        create_matrix(datatype, &V_test, ldvt, n_A);
        create_realtype_vector(datatype, &s_test, min_m_n);
        create_vector(datatype, &work, lwork);
        create_vector(INTEGER, &iwork, liwork);

        if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            create_realtype_vector(datatype, &rwork, lrwork);
        else
            rwork = NULL;


        exe_time = fla_test_clock();

        /* call to API */
        invoke_gesvd(datatype, jobu, jobvt, &m_A, &n_A, A, &lda, s_test, U_test, &ldu, V_test, &ldvt, work, &lwork, rwork, iwork, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Make a copy of the output buffers. This is required to validate the API functionality. */
        copy_matrix(datatype, "full", m_A, m_A, U_test, ldu, U, ldu);
        copy_matrix(datatype, "full", n_A, n_A, V_test, ldvt, V, ldvt);
        copy_realtype_vector(datatype, min_m_n, s_test, 1, s, 1);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(iwork);
        if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            free_vector(rwork);
        free_matrix(U_test);
        free_matrix(V_test);
        free_vector(s_test);
    }

    *time_min_ = time_min;

    free_matrix(A_save);
}


void invoke_gesvd(integer datatype, char* jobu, char* jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* iwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
            break;
        }
        
        case DOUBLE:
        {
            fla_lapack_dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
            break;
        }
    }
}
