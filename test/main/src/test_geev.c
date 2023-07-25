/*
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_geev_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, integer einfo, double* perf, double* t, double* residual);
void prepare_geev_run(char *jobvl, char *jobvr, integer n, void *a, integer lda, void *wr, void *wi, void *w, void *vl, integer ldvl, void *vr, integer ldvr, integer datatype, integer n_repeats, double* time_min_, integer* info);
void invoke_geev(integer datatype, char *jobvl, char *jobvr, integer *n, void *a, integer *lda, void *wr, void *wi, void *w, void *vl, integer *ldvl, void *vr, integer *ldvr, void* work, integer* lwork, void* rwork, integer* info);

void fla_test_geev(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Eigen Decomposition of non symmetric matrix";
    char* front_str = "GEEV";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_NSYM, fla_test_geev_experiment);
        tests_not_run = 0;
    }
    if(argc == 12)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[11]);
    }
    if(argc >= 11 && argc <= 12)
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
        params->eig_non_sym_paramslist[0].ldvl = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->eig_non_sym_paramslist[0].ldvr = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        g_lwork = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_geev_experiment(params, datatype,
                                          N, N,
                                          0,
                                          n_repeats, einfo,
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
        printf("\nIllegal arguments for geev\n");
        printf("./<EXE> geev <precisions - sdcz> <jobvl> <jobvr> <N> <LDA> <LDVL> <LDVR> <LWORK> <repeats>\n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
    if(g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
        g_ext_fptr = NULL;
    }
    return;

}

void fla_test_geev_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    integer einfo,
    double* perf,
    double *time_min,
    double* residual)
{
    integer n, lda, ldvl, ldvr;
    integer info = 0, vinfo = 0;
    void *A = NULL, *wr = NULL, *wi = NULL, *w = NULL, *VL = NULL, *VR = NULL;
    void *A_test = NULL;
    char jobvl, jobvr;

    /* Get input matrix dimensions.*/
    n = p_cur;
    lda = params->eig_non_sym_paramslist[pci].lda;
    ldvl = params->eig_non_sym_paramslist[pci].ldvl;
    ldvr = params->eig_non_sym_paramslist[pci].ldvr;

    *residual =  params->eig_non_sym_paramslist[pci].GenNonSymEigProblem_threshold;
    jobvl = params->eig_non_sym_paramslist[pci].jobvsl;
    jobvr = params->eig_non_sym_paramslist[pci].jobvsr;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if (config_data)
    {
        if (lda == -1)
        {
            lda = fla_max(1,n);
        }
        /* LDVL >= 1; if JOBVL = 'V', LDVL >= N */
        if (ldvl == -1)
        {
            if (jobvl == 'V')
            {
                ldvl = n;
            }
            else
            {
                ldvl = 1;
            }
        }
        /* LDVR >= 1; if JOBVR = 'V', LDVR >= N */
        if (ldvr == -1)
        {
            if (jobvr == 'V')
            {
                ldvr = n;
            }
            else
            {
                ldvr = 1;
            }
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, &A, lda, n);
    create_matrix(datatype, &VL, ldvl, n);
    create_matrix(datatype, &VR, ldvr, n);
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_vector(datatype, &w, n);
    }
    else
    {
        create_vector(datatype, &wr, n);
        create_vector(datatype, &wi, n);
    }

    init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);

    /* Make a copy of input matrix A. This is required to validate the API functionality. */
    create_matrix(datatype, &A_test, lda, n);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    prepare_geev_run(&jobvl, &jobvr, n, A_test, lda, wr, wi, w,  VL, ldvl, VR, ldvr, datatype, n_repeats, time_min, &info);

    /* performance computation
       4/3 n^3 flops if job = 'N'
       8/3 n^3 flops if job = 'V' */

    if(jobvl == 'N' && jobvr == 'N')
        *perf = (double)((4.0 / 3.0) * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    else
        *perf = (double)((8.0 / 3.0) * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if (info == 0)
        validate_geev(&jobvl, &jobvr, n, A, A_test, lda, VL, ldvl, VR, ldvr, w, wr, wi, datatype, residual, &vinfo);

    FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(VL);
    free_matrix(VR);
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

void prepare_geev_run(char *jobvl, char *jobvr, integer m_A, void *A, integer lda,
                        void *wr, void *wi, void *w,
                        void *VL, integer ldvl, void *VR, integer ldvr,
                        integer datatype, integer n_repeats, double* time_min_, integer* info)
{
    void *A_save = NULL, *rwork = NULL, *work = NULL;
    integer lwork, lrwork;
    integer i;
    double time_min = 1e9, exe_time;


    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, lda, m_A);
    copy_matrix(datatype, "full", m_A, m_A, A, lda, A_save, lda);

    /* Get rwork and iwork array size since it is not depedent on internal blocks*/
    lrwork = 2 * m_A;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_realtype_vector(datatype, &rwork, lrwork);
    }


    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.*/
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
 
        /* call to  geev API */
        invoke_geev(datatype, jobvl, jobvr, &m_A, NULL, &lda, NULL, NULL, NULL,
                    NULL, &ldvl, NULL, &ldvr, work, &lwork, rwork, info);
        if (*info == 0)
        {
            /* Get work size */
            lwork = get_work_value( datatype, work );
        }

        /* Output buffers will be freshly allocated for each iterations, free up
        the current output buffers.*/
        free_vector(work);
    }
    else
    {
         lwork = g_lwork;
    }
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        free_vector(rwork);
    }

    *info = 0;
    for (i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m_A, m_A, A_save, lda, A, lda);

        create_vector(datatype, &work, lwork);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            create_realtype_vector(datatype, &rwork, lrwork);
        }
        exe_time = fla_test_clock();

        /* call to geev API */
        invoke_geev(datatype, jobvl, jobvr, &m_A, A, &lda, wr, wi, w,
                    VL, &ldvl, VR, &ldvr, work, &lwork, rwork, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
        if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            free_vector(rwork);
        }
    }

    *time_min_ = time_min;

    free_matrix(A_save);
}


void invoke_geev(integer datatype, char *jobvl, char *jobvr, integer *n,
	       	void *a, integer *lda, void *wr, void *wi, void *w,
		void *vl, integer *ldvl, void *vr, integer *ldvr, 
		void* work, integer* lwork, void* rwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
            break;
        }
    }
}
