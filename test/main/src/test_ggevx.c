/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"


/* Local prototypes */
void fla_test_ggevx_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_ggevx_run(char* balanc, char* jobvl, char* jobvr, char* sense, integer n, 
                        void* a, void* b, void* alpha, void* alphar, void* alphai, void* beta, void* vl, void* vr,
                        integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, 
                        void* rcondv, integer datatype, integer n_repeats, double* time_min_);
void invoke_ggevx(integer datatype, char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, 
                        void* b, integer* ldb, void* alpha, void* alphar, void* alphai, void* beta, void* vl, integer* ldvl, 
                        void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, 
                        void* rconde, void* rcondv, void* work, integer* lwork, void* rwork, integer* iwork, integer* bwork, integer* info);


void fla_test_ggevx(test_params_t *params)
{
    char* op_str = "Computing Eigen value and Eigen vectors with condition numbers";
    char* front_str = "GGEVX";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_NSYM, fla_test_ggevx_experiment);
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
    integer m, n, cs_A,ldvl, ldvr;
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
    m = p_cur;
    n = q_cur;
    cs_A = n;
    ldvl = n;
    ldvr = n;

    /* Create input matrix parameters */
    create_matrix(datatype, &A, n, n);
    create_matrix(datatype, &B, n, n);
    create_matrix(datatype, &VL, ldvl, n);
    create_matrix(datatype, &VR, ldvr, n);
    create_vector(datatype, &beta, n);

    if (datatype == FLOAT || datatype == DOUBLE)
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


    /* Initialize input matrix A with random numbers */
    rand_matrix(datatype, A, n, n, cs_A);
    rand_matrix(datatype, B, n, n, cs_A);

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, &A_test, n, n);
    create_matrix(datatype, &B_test, n, n);
    copy_matrix(datatype, "full", n, n, A, cs_A, A_test, cs_A);
    copy_matrix(datatype, "full", n, n, B, cs_A, B_test, cs_A);

    prepare_ggevx_run(&BALANC, &JOBVL, &JOBVR, &SENSE, n, A, B, alpha, alphar, alphai, beta, VL, VR,
                        &ilo, &ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, datatype, n_repeats, &time_min);

    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2mn^2 - (2/3)n^3 flops */
    *perf = (double)((2.0 * m * m * m) - ((2.0 / 3.0) * m * m * m)) / time_min / FLOPS_PER_UNIT_PERF;
    if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if (JOBVL == 'V' || JOBVR == 'V')
    {
        validate_ggevx(&BALANC, &JOBVL, &JOBVR, &SENSE, n, A, B, alpha, alphar, alphai, beta, VL, ldvl, VR, ldvr, datatype, residual);
    }

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(B);
    free_matrix(B_test);
    free_matrix(VL);
    free_matrix(VR);
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
    free_vector(lscale);
    free_vector(rscale);
    free_vector(rconde);
    free_vector(rcondv);
    free_vector(abnrm);
    free_vector(bbnrm);

}

void prepare_ggevx_run(char* balanc, char* jobvl, char* jobvr, char* sense, integer n_A, 
                        void* A, void* B, void* alpha, void* alphar, void* alphai, void* beta, void* VL, void* VR,
                        integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, 
                        void* rcondv, integer datatype, integer n_repeats, double* time_min_)
{
    void *A_save = NULL, *B_save = NULL , *work = NULL, *rwork = NULL, *iwork = NULL, *bwork=NULL;;
    integer i, lda, ldb, ldvl, ldvr;
    integer lwork, info = 0;
    double time_min = 1e9, exe_time;


    lda = n_A;
    ldb = n_A;
    ldvl = n_A;
    ldvr = n_A;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, n_A, n_A);
    copy_matrix(datatype, "full", n_A, n_A, A, n_A, A_save, n_A);
    create_matrix(datatype, &B_save, n_A, n_A);
    copy_matrix(datatype, "full", n_A, n_A, B, n_A, B_save, n_A);

    if(*sense != 'E')
    {
        create_vector(INTEGER, &iwork, 6 + n_A);
    }
    if (*sense != 'N')
    {
        create_vector(INTEGER, &bwork, n_A);
    }


    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size. */
    lwork = -1;
    create_vector(datatype, &work, 8*n_A);

    /* call to  ggevx API */
    invoke_ggevx(datatype, balanc, jobvl, jobvr, sense, &n_A, A, &lda, B, &ldb, alpha, alphar, alphai, beta, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork, &info);

    /* Get work size */
    lwork = get_work_value( datatype, work);

    /* Output buffers will be freshly allocated for each iterations, free up 
       the current output buffers.*/ 
    free_vector(work);

    create_realtype_vector(datatype, &rwork, 8 * n_A);


    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers for each iteration */
        copy_matrix(datatype, "full", n_A, n_A, A_save, n_A, A, n_A);
        copy_matrix(datatype, "full", n_A, n_A, B_save, n_A, B, n_A);
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

    copy_matrix(datatype, "full", n_A, n_A, A_save, n_A, A, n_A);
    copy_matrix(datatype, "full", n_A, n_A, B_save, n_A, B, n_A);

    free_matrix(A_save);
    free_matrix(B_save);
    if (*sense != 'E')
    {
        free_vector(iwork);
    }
    if (*sense != 'N')
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
