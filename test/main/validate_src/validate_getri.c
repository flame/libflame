/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_getri.c
 *  @brief Defines validate function of GETRI() to use in test suite.
 *  */

#include "test_common.h"

void validate_getri(integer m_A,
    integer n_A,
    void* A,
    void* A_inv,
    integer lda,
    integer* IPIV,
    integer datatype,
    double* residual)
{
    /* System generated locals */
    void *I, * work;
    
    /* Create Identity matrix */
    create_matrix(datatype, &I, n_A, n_A);
    create_vector(datatype, &work, 2 * m_A);

    switch (datatype)
    {
        case FLOAT:
        {
            float res, norm, norm_A, norm_A_I, eps;

            eps = slamch_("Epsilon");
            norm_A = slange_("1", &m_A, &m_A, A, &lda, work);
            norm_A_I = slange_("1", &m_A, &m_A, A_inv, &lda, work);
            res = (1 / norm_A) / norm_A_I;
            /* compute I - A' * A */
            slaset_("full", &m_A, &m_A, &s_zero, &s_one, I, &m_A);
            /* compute I - A' * A */
            sgemm_("N", "N", &m_A, &n_A, &m_A, &s_n_one, A_inv, &lda, A, &lda, &s_one, I, &m_A);

            norm = slange_("1", &m_A, &m_A, I, &m_A, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            *residual = (norm * res) / eps / (double)n_A;
            break;
        }

        case DOUBLE:
        {
            double res, norm, norm_A, norm_A_I, eps;

            eps = dlamch_("Epsilon");
            norm_A = dlange_("1", &m_A, &m_A, A, &lda, work);
            norm_A_I = dlange_("1", &m_A, &m_A, A_inv, &lda, work);
            res = (1 / norm_A) / norm_A_I;
            /* compute I - A' * A */
            dlaset_("full", &m_A, &m_A, &d_zero, &d_one, I, &m_A);
            /* compute I - A' * A */
            dgemm_("N", "N", &m_A, &n_A, &m_A, &d_n_one, A_inv, &lda, A, &lda, &d_one, I, &m_A);

            norm = dlange_("1", &m_A, &m_A, I, &m_A, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            *residual = (norm * res) / eps / (double)n_A;
            break;
        }
        case COMPLEX:
        {
            float res, norm, norm_A, norm_A_I, eps;

            eps = slamch_("Epsilon");
            norm_A = clange_("1", &m_A, &m_A, A, &lda, work);
            norm_A_I = clange_("1", &m_A, &m_A, A_inv, &lda, work);
            res = (1 / norm_A) / norm_A_I;
            /* compute I - A' * A */
            claset_("full", &m_A, &m_A, &c_zero, &c_one, I, &m_A);
            /* compute I - A' * A */
            cgemm_("N", "N", &m_A, &n_A, &m_A, &c_n_one, A_inv, &lda, A, &lda, &c_one, I, &m_A);

            norm = clange_("1", &m_A, &m_A, I, &m_A, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            *residual = (norm * res) / eps / (double)n_A;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double res, norm, norm_A, norm_A_I, eps;

            eps = dlamch_("Epsilon");
            norm_A = zlange_("1", &m_A, &m_A, A, &lda, work);
            norm_A_I = zlange_("1", &m_A, &m_A, A_inv, &lda, work);
            res = (1 / norm_A) / norm_A_I;
            /* compute I - A' * A */
            zlaset_("full", &m_A, &m_A, &z_zero, &z_one, I, &m_A);
            /* compute I - A' * A */
            zgemm_("N", "N", &m_A, &n_A, &m_A, &z_n_one, A_inv, &lda, A, &lda, &z_one, I, &m_A);

            norm = zlange_("1", &m_A, &m_A, I, &m_A, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            *residual = (norm * res) / eps / (double)n_A;
            break;
        }
    }

    // Free up buffers
    free_vector(work);
    free_vector(I);
}