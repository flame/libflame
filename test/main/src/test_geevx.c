/*
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_geevx_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_geevx_run(char *balanc, char *jobvl, char *jobvr, char * sense, integer n, void *a, integer lda, void *wr, void *wi, void *w,
                       void *vl, integer ldvl, void *vr, integer ldvr, integer *ilo, integer * ihi, void *scale, void *abnrm,
		       void *rconde, void *rcondv, integer datatype, integer n_repeats, double* time_min_, integer* info);
void invoke_geevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char * sense, integer *n, void *a, integer *lda,
                  void *wr, void *wi, void *w, void *vl, integer *ldvl, void *vr, integer *ldvr, integer *ilo, integer *ihi,
                  void *scale, void *abnrm, void *rconde, void *rcondv, void* work, integer* lwork, void* rwork, integer* iwork,
		  integer* info);

/* Flag to indicate lwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
static integer g_lwork;
static FILE* g_ext_fptr = NULL;

void fla_test_geevx(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Eigen Decomposition of non symmetric matrix";
    char* front_str = "GEEVX";
    integer tests_not_run = 1, invalid_dtype = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT,  params, EIG_NSYM, fla_test_geevx_experiment);
        tests_not_run = 0;
    }
    if (argc == 14)
    {
        /* Read matrix input data from a file */
        g_ext_fptr = fopen(argv[13], "r");
        if (g_ext_fptr == NULL)
        {
            printf("\n Invalid input file argument \n");
            return;
        }
    }
    if (argc >= 13 && argc <= 14)
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
        params->eig_non_sym_paramslist[0].ldvl = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        params->eig_non_sym_paramslist[0].ldvr = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        g_lwork = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_geevx_experiment(params, datatype,
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
        printf("\nIllegal arguments for geevx\n");
        printf("./<EXE> geevx <precisions - sdcz> <balanc> <jobvl> <jobvr> <sense> <N> <LDA> <LDVL> <LDVR> <LWORK> <repeats>\n");
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

void fla_test_geevx_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    double* perf,
    double *time_min,
    double* residual)
{
    integer m, lda, ldvl, ldvr;
    integer info = 0, vinfo = 0;
    integer ilo, ihi;
    void *A = NULL, *wr = NULL, *wi = NULL, *w = NULL, *VL = NULL, *VR = NULL;
    void *scale = NULL, *abnrm = NULL, *rconde = NULL, *rcondv = NULL;
    void *A_test = NULL;
    char balanc, jobvl, jobvr, sense;

    /* Get input matrix dimensions.*/
    m = p_cur;
    lda = params->eig_non_sym_paramslist[pci].lda;
    ldvl = params->eig_non_sym_paramslist[pci].ldvl;
    ldvr = params->eig_non_sym_paramslist[pci].ldvr;

    if(lda < m || ldvl < m || ldvr < m)
    {
        *residual = DBL_MIN;
        return;
    }

    *residual =  params->eig_non_sym_paramslist[pci].GenNonSymEigProblem_threshold;
    balanc = params->eig_non_sym_paramslist[pci].balance_ggevx;
    jobvl = params->eig_non_sym_paramslist[pci].jobvsl;
    jobvr = params->eig_non_sym_paramslist[pci].jobvsr;
    sense = params->eig_non_sym_paramslist[pci].sense_ggevx;
    if(sense == 'B' || sense == 'E')
    {
        jobvl = 'V';
	    jobvr = 'V';
    }
    /* Create input matrix parameters */
    create_matrix(datatype, &A, lda, m);

    create_matrix(datatype, &VL, ldvl, m);
    create_matrix(datatype, &VR, ldvr, m);
    create_realtype_vector(datatype, &scale, m);
    create_realtype_vector(datatype, &abnrm, 1);
    create_realtype_vector(datatype, &rconde, m);
    create_realtype_vector(datatype, &rcondv, m);

    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_vector(datatype, &w, m);
    }
    else
    {
        create_vector(datatype, &wr, m);
        create_vector(datatype, &wi, m);
    }
    
    if (g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, m, m, lda, g_ext_fptr);
    }
    else
    {
        /* Initialize input matrix with random numbers */
        rand_matrix(datatype, A, m, m, lda);
    }

    /* Make a copy of input matrix A. This is required to validate the API functionality. */
    create_matrix(datatype, &A_test, lda, m);
    copy_matrix(datatype, "full", m, m, A, lda, A_test, lda);

    prepare_geevx_run(&balanc, &jobvl, &jobvr, &sense, m, A_test, lda, wr, wi, w,  VL, ldvl, VR, ldvr,
                      &ilo, &ihi, scale, abnrm, rconde , rcondv, datatype, n_repeats, time_min, &info);

    /* performance computation
       4/3 m^3 flops if job = 'N'
       8/3 m^3 + m^2 flops if job = 'V' */
    if(jobvl == 'N' && jobvr == 'N')
        *perf = (double)((4.0 / 3.0) * m * m * m) / *time_min / FLOPS_PER_UNIT_PERF;
    else
        *perf = (double)(((8.0 / 3.0) * m * m * m) + (m * m)) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if (info == 0)
        validate_geevx(&jobvl, &jobvr, &sense, &balanc, m, A, A_test, lda, VL, ldvl, VR, ldvr, w, wr, wi, scale,
                   abnrm, rconde, rcondv, datatype, residual, &vinfo);

    /* Assigning bigger value to residual as execution fails */
    if(info < 0 || vinfo < 0)
        *residual = DBL_MAX;

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(VL);
    free_matrix(VR);
    free_vector(scale);
    free_vector(abnrm);
    free_vector(rconde);
    free_vector(rcondv);
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

void prepare_geevx_run(char *balanc, char *jobvl, char *jobvr, char * sense,
                        integer m_A, void *A, integer lda,
                        void *wr, void *wi, void *w,
                        void *VL, integer ldvl, void *VR, integer ldvr,
                        integer *ilo, integer *ihi, void *scale, void *abnrm,
                        void *rconde, void *rcondv,
                        integer datatype, integer n_repeats, double* time_min_, integer* info)
{
    void *A_save = NULL, *rwork = NULL, *iwork = NULL, *work = NULL;
    integer lwork, liwork, lrwork;
    integer i;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, lda, m_A);
    copy_matrix(datatype, "full", m_A, m_A, A, lda, A_save, lda);

    /* Get rwork and iwork array size since it is not depedent on internal blocks*/
    lrwork = 2 * m_A;
    liwork = 2 * m_A - 2;

    if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_realtype_vector(datatype, &rwork, lrwork);
    }

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.*/
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
        /* call to  geevx API */
        invoke_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, NULL, &lda,
                    NULL, NULL, NULL, NULL, &ldvl, NULL, &ldvr,
                    ilo, ihi, NULL, NULL, NULL, NULL, work, &lwork, rwork, NULL, info);
        if(*info < 0)
        {
            free_matrix(A_save);
            free_vector(work);
            return;
        }

        /* Get work size */
        lwork = get_work_value( datatype, work );

        /* Output buffers will be freshly allocated for each iterations, free up
        the current output buffers.*/
        free_vector(work);
    }
    else
    {
         lwork = g_lwork;
    }
    
    if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        free_vector(rwork);
    }

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
        else
        {
            create_vector(INTEGER, &iwork, liwork);
        }

        exe_time = fla_test_clock();

        /* call to geevx API */
        invoke_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, A, &lda, wr, wi, w, VL, &ldvl, VR, &ldvr,
                     ilo, ihi, scale, abnrm, rconde, rcondv, work, &lwork, rwork, iwork, info);
 
        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
        if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            free_vector(rwork);
        }
        else
        {
            free_vector(iwork);
        }
    }
    *time_min_ = time_min;

    free_matrix(A_save);
}

void invoke_geevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char *sense,
                            integer *n, void *a, integer *lda, void *wr, void *wi, void *w,
                            void *vl, integer *ldvl, void *vr, integer *ldvr, integer *ilo, integer *ihi,
                            void *scale, void *abnrm, void *rconde, void *rcondv,
                            void* work, integer* lwork, void* rwork, integer* iwork, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgeevx(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info);
            break;
        }
        
        case DOUBLE:
        {
            fla_lapack_dgeevx(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgeevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgeevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info);
            break;
        }
    }
}
