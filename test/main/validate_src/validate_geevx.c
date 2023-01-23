/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_geevx.c
 *  @brief Defines validate function of GEEVX() to use in test suite.
 *  */

#include "test_common.h"

/* TODO validation of balanced matrix and condition numbers for eigen values and eigen vectors */ 
void validate_geevx(char* jobvl, char* jobvr,
    char* sense,
    char* balanc,
    integer m,
    void* A,
    void* A_test,
    integer lda,
    void* VL,
    integer ldvl,
    void* VR,
    integer ldvr,
    void* w,
    void* wr,
    void* wi,
    void* scale,
    void* abnrm,
    void* rconde,
    void* rcondv,
    integer datatype,
    double* residual,
    integer* info)
{
    void *work = NULL;
    void *lambda = NULL, *Vlambda = NULL;
    *info = 0;

    create_matrix(datatype, &lambda, m, m);
    create_matrix(datatype, &Vlambda, m, m);

    reset_matrix(datatype, m, m, lambda, m);
    reset_matrix(datatype, m, m, Vlambda, m);

    if (datatype == FLOAT || datatype == DOUBLE)
    {
        create_block_diagonal_matrix(datatype, wr, wi, lambda, m, m, m);
    }

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");
            if(*jobvl == 'V' && *jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                sgemm_("N", "N", &m, &m, &m, &s_one, A, &lda, VR, &ldvr, &s_zero, Vlambda, &m);
                norm_A = fla_lapack_slange("1", &m, &m, Vlambda, &m, work);
                sgemm_("N", "N", &m, &m, &m, &s_one, VR, &ldvr, lambda, &m, &s_n_one, Vlambda, &m);
                norm = fla_lapack_slange("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (float)m);

                /* Test 2
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                sgemm_("C", "N", &m, &m, &m, &s_one, A, &lda, VL, &ldvl, &s_zero, Vlambda, &m);
                norm_A = fla_lapack_slange("1", &m, &m, Vlambda, &m, work);
                sgemm_("N", "C", &m, &m, &m, &s_one, VL, &ldvl, lambda, &m, &s_n_one, Vlambda, &m);
                norm = fla_lapack_slange("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (float)m);
                *residual = (double)fla_max(resid1, resid2);
            }
            else if(*jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                sgemm_("N", "N", &m, &m, &m, &s_one, A, &lda, VR, &ldvr, &s_zero, Vlambda, &m);
                norm_A = fla_lapack_slange("1", &m, &m, Vlambda, &m, work);
                sgemm_("N", "N", &m, &m, &m, &s_one, VR, &ldvr, lambda, &m, &s_n_one, Vlambda, &m);
                norm = fla_lapack_slange("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (float)m);
                *residual = (double)resid1;
            }
            else if(*jobvl == 'V')
            {
                /* Test 1
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                sgemm_("C", "N", &m, &m, &m, &s_one, A, &lda, VL, &ldvl, &s_zero, Vlambda, &m);
                norm_A = fla_lapack_slange("1", &m, &m, Vlambda, &m, work);
                sgemm_("N", "C", &m, &m, &m, &s_one, VL, &ldvl, lambda, &m, &s_n_one, Vlambda, &m);
                norm = fla_lapack_slange("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (float)m);
                *residual = (double)resid1;
            }
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            if(*jobvl == 'V' && *jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                dgemm_("N", "N", &m, &m, &m, &d_one, A, &lda, VR, &ldvr, &d_zero, Vlambda, &m);
                norm_A = fla_lapack_dlange("1", &m, &m, Vlambda, &m, work);
                dgemm_("N", "N", &m, &m, &m, &d_one, VR, &ldvr, lambda, &m, &d_n_one, Vlambda, &m);
                norm = fla_lapack_dlange("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (double)m);

                /* Test 2
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                dgemm_("C", "N", &m, &m, &m, &d_one, A, &lda, VL, &ldvl, &d_zero, Vlambda, &m);
                norm_A = fla_lapack_dlange("1", &m, &m, Vlambda, &m, work);
                dgemm_("N", "C", &m, &m, &m, &d_one, VL, &ldvl, lambda, &m, &d_n_one, Vlambda, &m);
                norm = fla_lapack_dlange("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (double)m);
                *residual = (double)fla_max(resid1, resid2);
            }
            else if(*jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                dgemm_("N", "N", &m, &m, &m, &d_one, A, &lda, VR, &ldvr, &d_zero, Vlambda, &m);
                norm_A = fla_lapack_dlange("1", &m, &m, Vlambda, &m, work);
                dgemm_("N", "N", &m, &m, &m, &d_one, VR, &ldvr, lambda, &m, &d_n_one, Vlambda, &m);
                norm = fla_lapack_dlange("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (double)m);
                *residual = (double)resid1;
            }
            else if(*jobvl == 'V')
            {
                /* Test 1
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                dgemm_("C", "N", &m, &m, &m, &d_one, A, &lda, VL, &ldvl, &d_zero, Vlambda, &m);
                norm_A = fla_lapack_dlange("1", &m, &m, Vlambda, &m, work);
                dgemm_("N", "C", &m, &m, &m, &d_one, VL, &ldvl, lambda, &m, &d_n_one, Vlambda, &m);
                norm = fla_lapack_dlange("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (double)m);
                *residual = (double)resid2;
            }
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2;
            integer incr;
            incr = m+1;
            eps = fla_lapack_slamch("P");
            ccopy_(&m, w, &i_one, lambda, &incr);

            if(*jobvl == 'V' && *jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                cgemm_("N", "N", &m, &m, &m, &c_one, A, &lda, VR, &ldvr, &c_zero, Vlambda, &m);
                norm_A = fla_lapack_clange("1", &m, &m, Vlambda, &m, work);
                cgemm_("N", "N", &m, &m, &m, &c_one, VR, &ldvr, lambda, &m, &c_n_one, Vlambda, &m);
                norm = fla_lapack_clange("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (float)m);

                /* Test 2
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                cgemm_("C", "N", &m, &m, &m, &c_one, A, &lda, VL, &ldvl, &c_zero, Vlambda, &m);
                norm_A = fla_lapack_clange("1", &m, &m, Vlambda, &m, work);
                cgemm_("N", "C", &m, &m, &m, &c_one, VL, &ldvl, lambda, &m, &c_n_one, Vlambda, &m);
                norm = fla_lapack_clange("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (float)m);
                *residual = (double)fla_max(resid1, resid2);
            }
            else if(*jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                cgemm_("N", "N", &m, &m, &m, &c_one, A, &lda, VR, &ldvr, &c_zero, Vlambda, &m);
                norm_A = fla_lapack_clange("1", &m, &m, Vlambda, &m, work);
                cgemm_("N", "N", &m, &m, &m, &c_one, VR, &ldvr, lambda, &m, &c_n_one, Vlambda, &m);
                norm = fla_lapack_clange("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (float)m);
                *residual = (double)resid1;
            }
            else if(*jobvl == 'V')
            {
                /* Test 1
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                cgemm_("C", "N", &m, &m, &m, &c_one, A, &lda, VL, &ldvl, &c_zero, Vlambda, &m);
                norm_A = fla_lapack_clange("1", &m, &m, Vlambda, &m, work);
                cgemm_("N", "C", &m, &m, &m, &c_one, VL, &ldvl, lambda, &m, &c_n_one, Vlambda, &m);
                norm = fla_lapack_clange("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (float)m);
                *residual = (double)resid2;
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2;
            integer incr;
            incr = m+1;
            eps = fla_lapack_dlamch("P");
            zcopy_(&m, w, &i_one, lambda, &incr);

            if(*jobvl == 'V' && *jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                zgemm_("N", "N", &m, &m, &m, &z_one, A, &lda, VR, &ldvr, &z_zero, Vlambda, &m);
                norm_A = fla_lapack_zlange("1", &m, &m, Vlambda, &m, work);
                zgemm_("N", "N", &m, &m, &m, &z_one, VR, &ldvr, lambda, &m, &z_n_one, Vlambda, &m);
                norm = fla_lapack_zlange("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (double)m);

                /* Test 2
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                zgemm_("C", "N", &m, &m, &m, &z_one, A, &lda, VL, &ldvl, &z_zero, Vlambda, &m);
                norm_A = fla_lapack_zlange("1", &m, &m, Vlambda, &m, work);
                zgemm_("N", "C", &m, &m, &m, &z_one, VL, &ldvl, lambda, &m, &z_n_one, Vlambda, &m);
                norm = fla_lapack_zlange("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (double)m);
                *residual = (double)fla_max(resid1, resid2);
            }
            else if(*jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                zgemm_("N", "N", &m, &m, &m, &z_one, A, &lda, VR, &ldvr, &z_zero, Vlambda, &m);
                norm_A = fla_lapack_zlange("1", &m, &m, Vlambda, &m, work);
                zgemm_("N", "N", &m, &m, &m, &z_one, VR, &ldvr, lambda, &m, &z_n_one, Vlambda, &m);
                norm = fla_lapack_zlange("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (double)m);
                *residual = (double)resid1;
            }
            else if(*jobvl == 'V')
            {
                /* Test 1
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                zgemm_("C", "N", &m, &m, &m, &z_one, A, &lda, VL, &ldvl, &z_zero, Vlambda, &m);
                norm_A = fla_lapack_zlange("1", &m, &m, Vlambda, &m, work);
                zgemm_("N", "C", &m, &m, &m, &z_one, VL, &ldvl, lambda, &m, &z_n_one, Vlambda, &m);
                norm = fla_lapack_zlange("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (double)m);
                *residual = (double)resid2;
            }
            break;
        }
    }
    free_matrix(lambda);
    free_matrix(Vlambda);
}
