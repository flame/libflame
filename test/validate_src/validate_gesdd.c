/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_gesdd.c
 *  @brief Defines validate function of GESDD() to use in test suite.
 *  */

#include "test_common.h"

void validate_gesdd(char *jobz, integer m, integer n, void* A, void* A_test, void* s, void* U, void* V, integer datatype, double *residual)
{
    void *sigma = NULL, *Usigma = NULL;
    void *Ibuff_U = NULL, *Ibuff_V = NULL;
    void *work = NULL;

    create_matrix(datatype, &sigma, m, n);
    create_matrix(datatype, &Usigma, m, n);
    reset_matrix(datatype, m, n, Usigma, m);

    create_matrix(datatype, &Ibuff_U, m, m);
    create_matrix(datatype, &Ibuff_V, n, n);
    reset_matrix(datatype, m, m, Ibuff_U, m);
    reset_matrix(datatype, n, n, Ibuff_V, n);

    diagonalize_vector(datatype, s, sigma, m, n, m);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2, resid3;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = slange_("1", &m, &n, A, &m, work);
            sgemm_("N", "N", &m, &n, &m, &s_one, U, &m, sigma, &m, &s_zero, Usigma, &m);
            sgemm_("N", "N", &m, &n, &n, &s_one, Usigma, &m, V, &n, &s_n_one, A, &m);
            norm = slange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            slaset_("full", &m, &m, &s_zero, &s_one, Ibuff_U, &m);
            sgemm_("N", "T", &m, &m, &m, &s_one, U, &m, U, &m, &s_n_one, Ibuff_U, &m);
            norm = slange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (float)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            slaset_("full", &n, &n, &s_zero, &s_one, Ibuff_V, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, V, &n, V, &n, &s_n_one, Ibuff_V, &n);
            norm = slange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (float)n);
            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2, resid3;
            eps = dlamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = dlange_("1", &m, &n, A, &m, work);
            dgemm_("N", "N", &m, &n, &m, &d_one, U, &m, sigma, &m, &d_zero, Usigma, &m);
            dgemm_("N", "N", &m, &n, &n, &d_one, Usigma, &m, V, &n, &d_n_one, A, &m);
            norm = dlange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (double)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            dlaset_("full", &m, &m, &d_zero, &d_one, Ibuff_U, &m);
            dgemm_("N", "T", &m, &m, &m, &d_one, U, &m, U, &m, &d_n_one, Ibuff_U, &m);
            norm = dlange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (double)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            dlaset_("full", &n, &n, &d_zero, &d_one, Ibuff_V, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, V, &n, V, &n, &d_n_one, Ibuff_V, &n);
            norm = dlange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (double)n);
            *residual = (double)max(resid1, max(resid2, resid3));
             break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2, resid3;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = clange_("1", &m, &n, A, &m, work);
            cgemm_("N", "N", &m, &n, &m, &c_one, U, &m, sigma, &m, &c_zero, Usigma, &m);
            cgemm_("N", "N", &m, &n, &n, &c_one, Usigma, &m, V, &n, &c_n_one, A, &m);
            norm = clange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            claset_("full", &m, &m, &c_zero, &c_one, Ibuff_U, &m);
            cgemm_("N", "C", &m, &m, &m, &c_one, U, &m, U, &m, &c_n_one, Ibuff_U, &m);
            norm = clange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (float)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            claset_("full", &n, &n, &c_zero, &c_one, Ibuff_V, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, V, &n, V, &n, &c_n_one, Ibuff_V, &n);
            norm = clange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (float)n);
            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2, resid3;
            eps = dlamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = zlange_("1", &m, &n, A, &m, work);
            zgemm_("N", "N", &m, &n, &m, &z_one, U, &m, sigma, &m, &z_zero, Usigma, &m);
            zgemm_("N", "N", &m, &n, &n, &z_one, Usigma, &m, V, &n, &z_n_one, A, &m);
            norm = zlange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (double)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            zlaset_("full", &m, &m, &z_zero, &z_one, Ibuff_U, &m);
            zgemm_("N", "C", &m, &m, &m, &z_one, U, &m, U, &m, &z_n_one, Ibuff_U, &m);
            norm = zlange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (double)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            zlaset_("full", &n, &n, &z_zero, &z_one, Ibuff_V, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, V, &n, V, &n, &z_n_one, Ibuff_V, &n);
            norm = zlange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (float)n);
            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }
    }
    free_matrix(sigma);
    free_matrix(Usigma);
    free_matrix(Ibuff_U);
    free_matrix(Ibuff_V);
}