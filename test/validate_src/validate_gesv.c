/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_gesv.c
 *  @brief Defines validate function of GESV() to use in test suite.
 *  */

#include "test_common.h"

void validate_gesv(integer n,
    integer nrhs,
    void* A,
    void* B,
    void* X,
    integer datatype,
    double* residual)
{
    void* work = NULL;
    integer ldx, ldb;

    ldx = n;
    ldb = n;

    switch (datatype)
    {
        case FLOAT:
        {
            float norm_b, norm, eps, resid;
        
            /* Test 1 */
            norm_b = slange_("1", &n, &nrhs, B, &ldb, work);
            eps = slamch_("P");

            /* Compute AX-B */ 
            sgemm_("N", "N", &n, &nrhs, &n, &s_one, A, &n, X, &ldx, &s_n_one, B, &ldb);
            norm = slange_("1", &n, &nrhs, B, &nrhs, work);
        
            resid = norm / (eps * norm_b * (float)n);
        
            *residual = (double)resid;
            break;
        }
        case DOUBLE:
        {
            double norm_b, norm, eps, resid;
        
            /* Test 1 */
            norm_b = dlange_("1", &n, &nrhs, B, &ldb, work);
            eps = dlamch_("P");
        
            /* Compute AX-B */ 
            dgemm_("N", "N", &n, &nrhs, &n, &d_one, A, &n, X, &ldx, &d_n_one, B, &ldb);
            norm = dlange_("1", &n, &nrhs, B, &nrhs, work);
        
            resid = norm / (eps * norm_b * (double)n);
        
            *residual = (double)resid;
            break;
        }
        case COMPLEX:
        {
            float norm_b, norm, eps, resid;
        
            /* Test 1 */
            norm_b = clange_("1", &n, &nrhs, B, &ldb, work);
            eps = slamch_("P");

            /* Compute AX-B */ 
            cgemm_("N", "N", &n, &nrhs, &n, &c_one, A, &n, X, &ldx, &c_n_one, B, &ldb);
            norm = clange_("1", &n, &nrhs, B, &nrhs, work);
        
            resid = norm / (eps * norm_b * (float)n);
        
            *residual = (double)resid;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_b, norm, eps, resid;
        
            /* Test 1 */
            norm_b = zlange_("1", &n, &nrhs, B, &ldb, work);
            eps = dlamch_("P");

            /* Compute AX-B */ 
            zgemm_("N", "N", &n, &nrhs, &n, &z_one, A, &n, X, &ldx, &z_n_one, B, &ldb);
            norm = zlange_("1", &n, &nrhs, B, &nrhs, work);
        
            resid = norm / (eps * norm_b * (double)n);
        
            *residual = (double)resid;
            break;
        }
    }
}

