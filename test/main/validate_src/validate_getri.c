/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
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
    double* residual,
    integer* info)
{
    /* System generated locals */
    void *a_temp, *work;
    *info = 0;
    
    /* Create Identity matrix */
    create_matrix(datatype, &a_temp, n_A, n_A);
    create_vector(datatype, &work, 2 * m_A);

    switch (datatype)
    {
        case FLOAT:
        {
            float res, norm, norm_A, norm_A_I, eps;

            eps = fla_lapack_slamch("Epsilon");
            norm_A = fla_lapack_slange("1", &m_A, &m_A, A, &lda, work);
            norm_A_I = fla_lapack_slange("1", &m_A, &m_A, A_inv, &lda, work);
            res = (1 / norm_A) / norm_A_I;
            /* compute I - A' * A */
            fla_lapack_slaset("full", &m_A, &m_A, &s_zero, &s_one, a_temp, &m_A);
            /* compute I - A' * A */
            sgemm_("N", "N", &m_A, &n_A, &m_A, &s_n_one, A_inv, &lda, A, &lda, &s_one, a_temp, &m_A);

            norm = fla_lapack_slange("1", &m_A, &m_A, a_temp, &m_A, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            *residual = (norm * res) / eps / (double)n_A;
            break;
        }

        case DOUBLE:
        {
            double res, norm, norm_A, norm_A_I, eps;

            eps = fla_lapack_dlamch("Epsilon");
            norm_A = fla_lapack_dlange("1", &m_A, &m_A, A, &lda, work);
            norm_A_I = fla_lapack_dlange("1", &m_A, &m_A, A_inv, &lda, work);
            res = (1 / norm_A) / norm_A_I;
            /* compute I - A' * A */
            fla_lapack_dlaset("full", &m_A, &m_A, &d_zero, &d_one, a_temp, &m_A);
            /* compute I - A' * A */
            dgemm_("N", "N", &m_A, &n_A, &m_A, &d_n_one, A_inv, &lda, A, &lda, &d_one,a_temp, &m_A);

            norm = fla_lapack_dlange("1", &m_A, &m_A, a_temp, &m_A, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            *residual = (norm * res) / eps / (double)n_A;
            break;
        }
        case COMPLEX:
        {
            float res, norm, norm_A, norm_A_I, eps;

            eps = fla_lapack_slamch("Epsilon");
            norm_A = fla_lapack_clange("1", &m_A, &m_A, A, &lda, work);
            norm_A_I = fla_lapack_clange("1", &m_A, &m_A, A_inv, &lda, work);
            res = (1 / norm_A) / norm_A_I;
            /* compute I - A' * A */
            fla_lapack_claset("full", &m_A, &m_A, &c_zero, &c_one, a_temp, &m_A);
            /* compute I - A' * A */
            cgemm_("N", "N", &m_A, &n_A, &m_A, &c_n_one, A_inv, &lda, A, &lda, &c_one,a_temp, &m_A);

            norm = fla_lapack_clange("1", &m_A, &m_A, a_temp, &m_A, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            *residual = (norm * res) / eps / (double)n_A;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double res, norm, norm_A, norm_A_I, eps;

            eps = fla_lapack_dlamch("Epsilon");
            norm_A = fla_lapack_zlange("1", &m_A, &m_A, A, &lda, work);
            norm_A_I = fla_lapack_zlange("1", &m_A, &m_A, A_inv, &lda, work);
            res = (1 / norm_A) / norm_A_I;
            /* compute I - A' * A */
            fla_lapack_zlaset("full", &m_A, &m_A, &z_zero, &z_one, a_temp, &m_A);
            /* compute I - A' * A */
            zgemm_("N", "N", &m_A, &n_A, &m_A, &z_n_one, A_inv, &lda, A, &lda, &z_one,a_temp, &m_A);

            norm = fla_lapack_zlange("1", &m_A, &m_A, a_temp, &m_A, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            *residual = (norm * res) / eps / (double)n_A;
            break;
        }
    }

    // Free up buffers
    free_vector(work);
    free_vector(a_temp);
}
