/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1


// Local prototypes.
void fla_test_geevx_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, double* perf, double* t, double* residual);
void prepare_geevx_run(char *balanc, char *jobvl, char *jobvr, char * sense, integer n, void *a, void *wr, void *wi, void *w,
							void *vl, integer ldvl, void *vr, integer ldvr, integer *ilo, integer * ihi, void *scale, void *abnrm, void *rconde,
							void *rcondv, integer datatype, integer n_repeats, double* time_min_);
void invoke_geevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char * sense, integer *n, void *a, integer *lda, void *wr, void *wi, void *w,
							void *vl, integer *ldvl, void *vr, integer *ldvr, integer *ilo, integer *ihi, void *scale, void *abnrm, void *rconde, void *rcondv,
							void* work, integer* lwork, void* rwork, integer* iwork, integer* info);

void fla_test_geevx(test_params_t *params)
{
    char* op_str = "Eigen value and Eigen vectors";
    char* front_str = "GEEVX";
    char* lapack_str = "LAPACK";
    char* pc_str[NUM_PARAM_COMBOS] = { "" };

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS,
                            params, EIG_SYM, fla_test_geevx_experiment);
}


void fla_test_geevx_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    double* perf,
    double *t,
    double* residual)
{
    integer m, n, cs_A, ldvl, ldvr;
    integer ilo, ihi;
    void *A = NULL, *wr = NULL, *wi = NULL, *w = NULL, *VL = NULL, *VR = NULL;
    void *scale = NULL, *abnrm = NULL, *rconde = NULL, *rcondv = NULL;
    void *A_test = NULL;
    double time_min = 1e9;


    // Get input matrix dimensions.
    m = p_cur;
    n = p_cur;
    cs_A = m;
    ldvl = m;
    ldvr = m;

    // Create input matrix parameters
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


    // Initialize input matrix A with random numbers
    rand_matrix(datatype, A, m, n, cs_A);

    // Make a copy of input matrix A. This is required to validate the API functionality.
    create_matrix(datatype, &A_test, m, n);
    copy_matrix(datatype, "full", m, n, A, cs_A, A_test, cs_A);

    //**prepare_gesdd_run(VECTORS_ALL, m, n, A_test, s, U, V, datatype, n_repeats, &time_min);
    prepare_geevx_run("B", "N", "N", "N", m, A_test, wr, wi, w,  VL, ldvl, VR, ldvr, &ilo, &ihi, scale, abnrm, rconde, rcondv, datatype, n_repeats, &time_min);

    // execution time
    *t = time_min;

    // Todo
    // performance computation will be added later
    *perf = 0;

    // Todo
    // output validation

    // Free up the buffers
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
    void *w_test = NULL, *wi_test = NULL, *wr_test = NULL, *scale_test = NULL, *abnrm_test = NULL;
    //void *VL_test = NULL, *VR_test = NULL, *rconde_test = NULL, *rcondv_test = NULL;
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

    // Make a workspace query the first time through. This will provide us with
    // and ideal workspace size based on an internal block size.
    lwork = -1;
    create_vector(datatype, &work, 1);

    // call to  geevx API
    invoke_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, NULL, &cs_A,
                    NULL, NULL, NULL, NULL, &ldvl, NULL, &ldvr,
                    ilo, ihi, NULL, NULL, NULL, NULL, work, &lwork, NULL, NULL, &info);

    // Get work size
    lwork = get_work_value( datatype, work );

    /* Output buffers will be freshly allocated for each iterations, free up 
       the current output buffers.*/ 
    free_vector(work);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m_A, m_A, A_save, m_A, A, m_A);

        //create_matrix(datatype, &VL_test, ldvl, m_A);
        //create_matrix(datatype, &VR_test, ldvr, m_A);
        create_realtype_vector(datatype, &scale_test, m_A);
        create_realtype_vector(datatype, &abnrm_test, 1);
        //create_realtype_vector(datatype, &rconde_test, m_A);
        //create_realtype_vector(datatype, &rcondv_test, m_A);
        create_vector(datatype, &work, lwork);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            create_vector(datatype, &w_test, m_A);
            create_realtype_vector(datatype, &rwork, lrwork);
        }
        else
        {
            create_vector(datatype, &wr_test, m_A);
            create_vector(datatype, &wi_test, m_A);
            create_vector(INTEGER, &iwork, liwork);
        }

        exe_time = fla_test_clock();

        // call to API
        /*invoke_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, A, &cs_A,
                    wr_test, wi_test, w_test, VL_test, &ldvl, VR_test, &ldvr,
                    ilo, ihi, scale_test, abnrm_test, rconde_test, rcondv_test, work, &lwork, rwork, iwork, &info);*/
        invoke_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, A, &cs_A,
                    wr_test, wi_test, w_test, NULL, &ldvl, NULL, &ldvr,
                    ilo, ihi, scale_test, abnrm_test, NULL, NULL, work, &lwork, rwork, iwork, &info);
        
        exe_time = fla_test_clock() - exe_time;

        // Get the best execution time
        time_min = min(time_min, exe_time);

        // Todo
        // Make a copy of the output buffers. This is required to validate the API functionality.


        // Free up the output buffers
        free_vector(work);
        if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            free_vector(w_test);
            free_vector(rwork);
        }
        else
        {
            free_vector(wr_test);
            free_vector(wi_test);
            free_vector(iwork);
        }
        //free_matrix(VL_test);
        //free_matrix(VR_test);
        free_vector(scale_test);
        free_vector(abnrm_test);
        //free_vector(rconde_test);
        //free_vector(rcondv_test);
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
