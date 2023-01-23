/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_gesv.c
 *  @brief Defines validate function of GESV() to use in test suite.
 *  */

#include "test_common.h"

void validate_gesv(integer n,
    integer nrhs,
    void* A,
    integer lda,
    void* B,
    integer ldb,
    void* X,
    integer datatype,
    double* residual,
    integer* info)
{
    void* work = NULL;
    integer ldx;
    *info = 0;
    ldx = ldb;

    switch (datatype)
    {
        case FLOAT:
        {
             float norm_a, norm_b, norm_x, norm, eps, resid;
        
            /* Test 1 */
             norm_a = fla_lapack_slange("1", &n, &n, A, &lda, work);
             norm_b = fla_lapack_slange("1", &n, &nrhs, B, &ldb, work);
             norm_x = fla_lapack_slange("1", &n, &nrhs, X, &ldx, work);
             eps = fla_lapack_slamch("E");

            /* Compute AX-B */ 
            sgemm_("N", "N", &n, &nrhs, &n, &s_one, A, &lda, X, &ldx, &s_n_one, B, &ldb);
            norm = fla_lapack_slange("1", &n, &nrhs, B, &ldb, work);
        
             resid = norm / (((norm_a * norm_x + norm_b) * (float)n) * eps);
        
            *residual = (double)resid;
            break;
        }
        case DOUBLE:
        {
             double norm_a, norm_b, norm_x, norm, eps, resid;
        
            /* Test 1 */
             norm_a = fla_lapack_dlange("1", &n, &n, A, &lda, work);
             norm_b = fla_lapack_dlange("1", &n, &nrhs, B, &ldb, work);
             norm_x = fla_lapack_dlange("1", &n, &nrhs, X, &ldx, work);
             eps = fla_lapack_dlamch("E");
        
            /* Compute AX-B */ 
            dgemm_("N", "N", &n, &nrhs, &n, &d_one, A, &lda, X, &ldx, &d_n_one, B, &ldb);
            norm = fla_lapack_dlange("1", &n, &nrhs, B, &ldb, work);
        
             resid = norm / (((norm_a * norm_x + norm_b) * (double)n) * eps);
        
            *residual = (double)resid;
            break;
        }
        case COMPLEX:
        {
             float norm_a, norm_b, norm_x, norm, eps, resid;
        
            /* Test 1 */
             norm_a = fla_lapack_clange("1", &n, &n, A, &lda, work);
             norm_b = fla_lapack_clange("1", &n, &nrhs, B, &ldb, work);
             norm_x = fla_lapack_clange("1", &n, &nrhs, X, &ldx, work);
             eps = fla_lapack_slamch("E");

            /* Compute AX-B */ 
            cgemm_("N", "N", &n, &nrhs, &n, &c_one, A, &lda, X, &ldx, &c_n_one, B, &ldb);
            norm = fla_lapack_clange("1", &n, &nrhs, B, &ldb, work);
        
             resid = norm / (((norm_a * norm_x + norm_b) * (float)n) * eps);
        
            *residual = (double)resid;
            break;
        }
        case DOUBLE_COMPLEX:
        {
             double norm_a, norm_b, norm_x, norm, eps, resid;
        
            /* Test 1 */
             norm_a = fla_lapack_zlange("1", &n, &n, A, &lda, work);
             norm_b = fla_lapack_zlange("1", &n, &nrhs, B, &ldb, work);
             norm_x = fla_lapack_zlange("1", &n, &nrhs, X, &ldx, work);
             eps = fla_lapack_dlamch("E");

            /* Compute AX-B */ 
            zgemm_("N", "N", &n, &nrhs, &n, &z_one, A, &lda, X, &ldx, &z_n_one, B, &ldb);
            norm = fla_lapack_zlange("1", &n, &nrhs, B, &ldb, work);
        
             resid = norm / (((norm_a * norm_x + norm_b) * (double)n) * eps);
        
            *residual = (double)resid;
            break;
        }
    }
}

