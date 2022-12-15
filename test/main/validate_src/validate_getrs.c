/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_getrs.c
 *  @brief Defines validate function of GETRS() to use in test suite.
 *  */

#include "test_common.h"

void validate_getrs(char *trans,
    integer n,
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
    ldx = n;

    switch (datatype)
    {
         case FLOAT:
         {
             float norm_a, norm_b, norm_x, norm, eps, resid;
         
             /* Test 1 */
             norm_a = slange_("1", &n, &n, A, &lda, work);
             norm_b = slange_("1", &n, &nrhs, B, &ldb, work);
             norm_x = slange_("1", &n, &nrhs, X, &ldx, work);
             eps = slamch_("E");
         
             /* Compute AX-B */
             sgemm_(trans, "N", &n, &nrhs, &n, &s_one, A, &lda, X, &ldx, &s_n_one, B, &ldb);
             norm = slange_("1", &n, &nrhs, B, &ldb, work);
         
             resid = norm / (((norm_a * norm_x + norm_b) * (float)n) * eps);
         
             *residual = (double)resid;
             break;
         }
         case DOUBLE:
         {
             double norm_a, norm_b, norm_x, norm, eps, resid;
         
             /* Test 1 */
             norm_a = dlange_("1", &n, &n, A, &lda, work);
             norm_b = dlange_("1", &n, &nrhs, B, &ldb, work);
             norm_x = dlange_("1", &n, &nrhs, X, &ldx, work);
             eps = dlamch_("E");
         
             /* Compute AX-B */
             dgemm_(trans, "N", &n, &nrhs, &n, &d_one, A, &lda, X, &ldx, &d_n_one, B, &ldb);
             norm = dlange_("1", &n, &nrhs, B, &ldb, work);
         
             resid = norm / (((norm_a * norm_x + norm_b) * (double)n) * eps);
         
             *residual = (double)resid;
             break;
         }
         case COMPLEX:
         {
             float norm_a, norm_b, norm_x, norm, eps, resid;
        
            /* Test 1 */
             norm_a = clange_("1", &n, &n, A, &lda, work);
             norm_b = clange_("1", &n, &nrhs, B, &ldb, work);
             norm_x = clange_("1", &n, &nrhs, X, &ldx, work);
             eps = slamch_("E");
         
             /* Compute AX-B */
             cgemm_(trans, "N", &n, &nrhs, &n, &c_one, A, &lda, X, &ldx, &c_n_one, B, &ldb);
             norm = clange_("1", &n, &nrhs, B, &ldb, work);
         
             resid = norm / (((norm_a * norm_x + norm_b) * (float)n) * eps);
         
             *residual = (double)resid;
             break;
         }
         case DOUBLE_COMPLEX:
         {
             double norm_a, norm_b, norm_x, norm, eps, resid;
         
             /* Test 1 */
             norm_a = zlange_("1", &n, &n, A, &lda, work);
             norm_b = zlange_("1", &n, &nrhs, B, &ldb, work);
             norm_x = zlange_("1", &n, &nrhs, X, &ldx, work);
             eps = dlamch_("E");
         
             /* Compute AX-B */
             zgemm_(trans, "N", &n, &nrhs, &n, &z_one, A, &lda, X, &ldx, &z_n_one, B, &ldb);
             norm = zlange_("1", &n, &nrhs, B, &ldb, work);
         
             resid = norm / (((norm_a * norm_x + norm_b) * (double)n) * eps);
         
             *residual = (double)resid;
             break;
         }
    }
}