/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_geevx_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_geevx_run(char *balanc, char *jobvl, char *jobvr, char * sense, integer n, void *a, void *wr, void *wi, void *w,
                       void *vl, integer ldvl, void *vr, integer ldvr, integer *ilo, integer * ihi, void *scale, void *abnrm,
		       void *rconde, void *rcondv, integer datatype, integer n_repeats, double* time_min_);
void invoke_geevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char * sense, integer *n, void *a, integer *lda,
                  void *wr, void *wi, void *w, void *vl, integer *ldvl, void *vr, integer *ldvr, integer *ilo, integer *ihi,
                  void *scale, void *abnrm, void *rconde, void *rcondv, void* work, integer* lwork, void* rwork, integer* iwork,
		  integer* info);


void fla_test_geevx(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Eigen Decomposition of non symmetric matrix";
    char* front_str = "GEEVX";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, SQUARE_INPUT,  params, EIG_NSYM, fla_test_geevx_experiment);
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
    integer m, n, cs_A, ldvl, ldvr;
    integer ilo, ihi;
    void *A = NULL, *wr = NULL, *wi = NULL, *w = NULL, *VL = NULL, *VR = NULL;
    void *scale = NULL, *abnrm = NULL, *rconde = NULL, *rcondv = NULL;
    void *A_test = NULL;
    char balanc, jobvl, jobvr, sense;

    /* Get input matrix dimensions.*/
    m = p_cur;
    n = p_cur;
    cs_A = m;
    ldvl = m;
    ldvr = m;
    
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
    create_matrix(datatype, &A, m, n);

    create_matrix(datatype, &VL, ldvl, m);
    create_matrix(datatype, &VR, ldvr, m);
    create_realtype_vector(datatype, &scale, m);
    create_realtype_vector(datatype, &abnrm, 1);
    create_realtype_vector(datatype, &rconde, m);
    create_realtype_vector(datatype, &rcondv, n);

    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_vector(datatype, &w, m);
    }
    else
    {
        create_vector(datatype, &wr, m);
        create_vector(datatype, &wi, m);
    }
    
    /* Initialize input matrix A with random numbers */
    rand_matrix(datatype, A, m, n, cs_A);

    /* Make a copy of input matrix A. This is required to validate the API functionality. */
    create_matrix(datatype, &A_test, m, n);
    copy_matrix(datatype, "full", m, n, A, cs_A, A_test, cs_A);

    prepare_geevx_run(&balanc, &jobvl, &jobvr, &sense, m, A_test, wr, wi, w,  VL, ldvl, VR, ldvr,
                      &ilo, &ihi, scale, abnrm, rconde , rcondv, datatype, n_repeats, time_min);

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
    validate_geevx(&jobvl, &jobvr, &sense, &balanc, m, A, A_test, VL, VR, w, wr, wi, scale,
                   abnrm, rconde, rcondv, datatype, residual);

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
                        integer m_A, void *A,
                        void *wr, void *wi, void *w,
                        void *VL, integer ldvl, void *VR, integer ldvr,
                        integer *ilo, integer *ihi, void *scale, void *abnrm,
                        void *rconde, void *rcondv,
                        integer datatype, integer n_repeats, double* time_min_)
{
    integer cs_A;
    void *A_save = NULL, *rwork = NULL, *iwork = NULL, *work = NULL;
    integer lwork, liwork, lrwork;
    integer i;
    integer info = 0;
    double time_min = 1e9, exe_time;

    cs_A = m_A;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, m_A, m_A);
    copy_matrix(datatype, "full", m_A, m_A, A, m_A, A_save, m_A);

    /* Get rwork and iwork array size since it is not depedent on internal blocks*/
    lrwork = 2 * m_A;
    liwork = 2 * m_A - 2;

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.*/
    lwork = -1;
    create_vector(datatype, &work, 1);
    if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_realtype_vector(datatype, &rwork, lrwork);
    }
    /* call to  geevx API */
    invoke_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, NULL, &cs_A,
                    NULL, NULL, NULL, NULL, &ldvl, NULL, &ldvr,
                    ilo, ihi, NULL, NULL, NULL, NULL, work, &lwork, rwork, NULL, &info);

    /* Get work size */
    lwork = get_work_value( datatype, work );

    /* Output buffers will be freshly allocated for each iterations, free up
       the current output buffers.*/
    free_vector(work);
    if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        free_vector(rwork);
    }

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m_A, m_A, A_save, m_A, A, m_A);

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
        invoke_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, A, &cs_A, wr, wi, w, VL, &ldvl, VR, &ldvr,
                     ilo, ihi, scale, abnrm, rconde, rcondv, work, &lwork, rwork, iwork, &info);
 
        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

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

    free(A_save);
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
            sgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info);
            break;
        }
        
        case DOUBLE:
        {
            dgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info);
            break;
        }

        case COMPLEX:
        {
            cgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info);
            break;
        }
    }
}
