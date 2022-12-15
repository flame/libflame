/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_getrf.c
 *  @brief Defines validate function of GETRF() to use in test suite.
 *  */

#include "test_common.h"

void validate_getrf(integer m_A,
    integer n_A,
    void* A,
    void* A_test,/*AFACT*/
    integer lda,
    integer* IPIV,
    integer datatype,
    double* residual,
    integer* info)
{
    /* System generated locals */
    integer m_n_vector, min_A;
    integer m_L, n_L, m_U, n_U, k;
    void *L, *U, *T, *work, *B, *B_test, *A_save;
    integer nrhs=1;
    *info = 0;

    m_n_vector = m_A * n_A;
    min_A = fla_min(m_A, n_A);
    create_vector(datatype, &B, m_A);
    create_vector(datatype, &B_test, m_A);
    if (m_A > n_A)
    {
        m_L = m_A;
        n_L = m_U = n_U = n_A;
        k = n_A;
    }
    else
    {
        m_L = n_L = m_U = m_A;
        n_U = n_A;
        k = m_A;
    }
    create_matrix(datatype, &L, m_L, n_L);
    create_matrix(datatype, &U, m_U, n_U);
    reset_matrix(datatype, m_U, n_U, U, m_U);
    create_matrix(datatype, &T, m_A, n_A);
    create_vector(datatype, &work, 2 * m_A);

    rand_vector(datatype, B, m_A, 1);
    copy_vector(datatype, m_A, B, 1, B_test, 1);

    /* Lower triangular matrix should be sqare matrix m x m */
    /* For m==i OR  m < n OR m > n -->  A(mxn) = L(mxm) * U(mxn) */
    copy_matrix(datatype, "Lower", m_L, n_L, A_test, lda, L, m_L);
    copy_matrix(datatype, "Upper", m_U, n_U, A_test, lda, U, m_U);
    
    // Create matrix of A(MxN) for subtracting in Test 2 (T-A).
    create_matrix(datatype, &A_save, m_A, n_A);
    copy_matrix(datatype, "Full", m_A, n_A, A, lda, A_save, m_A);

    switch (datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, norm_B, eps, resid1, resid2;
            eps = slamch_("Epsilon");
            /* Test 1 */
            if (m_A == n_A)
            {
                norm_B = snrm2_(&m_A, B, &i_one);
                /* Compute X by passing A and B */
                sgetrs_("N", &m_A, &nrhs, A_test, &lda, IPIV, B_test, &m_A, info);
                if(*info < 0)
                    break;
                /* Compute AX-B */
                sgemv_("N", &m_A, &m_A, &s_one, A, &lda, B_test, &i_one, &s_n_one, B, &i_one);
                norm = snrm2_(&m_A, B, &i_one);
                resid1 = norm / (float)m_A / norm_B / eps;
            }
            else
            {
                resid1 = 0.0;
            }
            
            /* Test 2 */
            /* Unity diagonal elements to Lower triangular matrix */
            slaset_("U", &m_L, &n_L, &s_zero, &s_one, L, &m_L);
            norm_A = slange_("1", &m_A, &n_A, A, &lda, work);
            /* T = L * U  */
            sgemm_("N", "N", &m_A, &n_A, &k, &s_one, L, &m_L, U, &m_U, &s_zero, T, &m_A);
            /*  Row interchanges based on IPIV values */
            slaswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
            /* T - A --> L*U - A */
            saxpy_(&m_n_vector, &s_n_one, A_save, &i_one, T, &i_one);
            /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
            norm = slange_("1", &m_A, &n_A, T, &m_A, work);

            resid2 = norm / (float)n_A / norm_A / eps;
            *residual = (float)fla_max(resid1, resid2);
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, norm_B, eps, resid1, resid2;

            eps = slamch_("Epsilon");
            /* Test 1 */
            if (m_A == n_A)
            {
                norm_B = dnrm2_(&m_A, B, &i_one);
                /* Compute X by passing A and B */
                dgetrs_("N", &m_A, &nrhs, A_test, &lda, IPIV, B_test, &m_A, info);
                if(*info < 0)
                    break;
                /* Compute AX-B */
                dgemv_("N", &m_A, &m_A, &d_one, A, &lda, B_test, &i_one, &d_n_one, B, &i_one);
                norm = dnrm2_(&m_A, B, &i_one);
                resid1 = norm / (double)m_A / norm_B / eps;
            }
            else
            {
                resid1 = 0.0;
            }
            
            /* Test 2 */
            /* Unity diagonal elements to Lower triangular matrix */
            dlaset_("U",&m_L, &n_L, &d_zero, &d_one, L, &m_L);
            norm_A = dlange_("1", &m_A, &n_A, A, &lda, work);
            /* T = L * U  */
            dgemm_("N", "N", &m_A, &n_A, &k, &d_one, L, &m_L, U, &m_U, &d_zero, T, &m_A);
            /*  Row interchanges based on IPIV values*/
            dlaswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
            /* T - A --> L*U - A */
            daxpy_(&m_n_vector, &d_n_one, A_save, &i_one, T, &i_one);
            /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
            norm = dlange_("1", &m_A, &n_A, T, &m_A, work);

            resid2 = norm / (double)n_A / norm_A / eps;
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, norm_B, eps, resid1, resid2;

            eps = slamch_("Epsilon");
            /* Test 1 */
            if (m_A == n_A)
            {
                norm_B = snrm2_(&m_A, B, &i_one);
                /* Compute X by passing A and B */
                cgetrs_("N", &m_A, &nrhs, A_test, &lda, IPIV, B_test, &m_A, info);
                if(*info < 0)
                    break;
                /* Compute AX-B */
                cgemv_("N", &m_A, &m_A, &c_one, A, &lda, B_test, &i_one, &c_n_one, B, &i_one);
                norm = snrm2_(&m_A, B, &i_one);
                resid1 = norm / (float)m_A / norm_B / eps;
            }
            else
            {
                resid1 = 0.0;
            }
            
            /* Test 2 */
            /* Unity diagonal elements to Lower triangular matrix */
            claset_("U", &m_L, &n_L, &c_zero, &c_one, L, &m_L);

            norm_A = clange_("1", &m_A, &n_A, A, &lda, work);
            /* T = L * U  */
            cgemm_("N", "N", &m_A, &n_A, &k, &c_one, L, &m_L, U, &m_U, &c_zero, T, &m_A);
            /*  Row interchanges based on IPIV values*/
            claswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
            /* T - A --> L*U - A */
            caxpy_(&m_n_vector, &c_n_one, A_save, &i_one, T, &i_one);
            /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
            norm = clange_("1", &m_A, &n_A, T, &m_A, work);

            resid2 = norm / (float)n_A / norm_A / eps;
            *residual = (float)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, norm_B, eps, resid1, resid2;

            eps = slamch_("Epsilon");
            /* Test 1 */
            if (m_A == n_A)
            {
                norm_B = dnrm2_(&m_A, B, &i_one);
                /* Compute X by passing A and B */
                zgetrs_("N", &m_A, &nrhs, A_test, &lda, IPIV, B_test, &m_A, info);
                if(*info < 0)
                    break;
                /* Compute AX-B */
                zgemv_("N", &m_A, &m_A, &z_one, A, &lda, B_test, &i_one, &z_n_one, B, &i_one);
                norm = dnrm2_(&m_A, B, &i_one);
                resid1 = norm / (double)m_A / norm_B / eps;
            }
            else
            {
                resid1 = 0.0;
            }
            
            /* Test 2 */
            /* Unity diagonal elements to Lower triangular matrix */
            zlaset_("U", &m_L, &n_L, &z_zero, &z_one, L, &m_L);

            norm_A = zlange_("1", &m_A, &n_A, A, &lda, work);
            /* T = L * U  */
            zgemm_("N", "N", &m_A, &n_A, &k, &z_one, L, &m_L, U, &m_U, &z_zero, T, &m_A);
            /*  Row interchanges based on IPIV values*/
            zlaswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
            /* T - A --> L*U - A */
            zaxpy_(&m_n_vector, &z_n_one, A_save, &i_one, T, &i_one);
            /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
            norm = zlange_("1", &m_A, &n_A, T, &m_A, work);

            resid2 = norm / (double)n_A / norm_A / eps;
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }

    // Free up buffers
    free_matrix(L);
    free_matrix(U);
    free_matrix(T);
    free_vector(work);
    free_vector(B);
    free_vector(B_test);
    free_vector(A_save);
}