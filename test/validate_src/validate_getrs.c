/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_getrs.c
 *  @brief Defines validate function of GETRS() to use in test suite.
 *  */

#include "test_common.h"

void validate_getrs(char *trans,
    integer m,
    integer nrhs,
    void* A,
    void* B,
    void* X,
    integer datatype,
    double* residual)
{
    void* work = NULL;
    integer ldx, ldb;

    ldx = m;
    ldb = m;
    /* Creating work buffer */
    create_vector(datatype, &work, 2 * m);

    switch (datatype)
    {
         case FLOAT:
         {
             float norm_b, norm, eps, resid;
         
             /* Test 1 */
             norm_b = slange_("1", &m, &nrhs, B, &nrhs, work);
             eps = slamch_("P");
         
             /* Compute AX-B */
             sgemm_(trans, "N", &m, &nrhs, &m, &s_one, A, &m, X, &ldx, &s_n_one, B, &ldb);
             norm = slange_("1", &m, &nrhs, B, &nrhs, work);
         
             resid = norm / (eps * norm_b * (float)m);
         
             *residual = (double)resid;
             break;
         }
         case DOUBLE:
         {
             double norm_b, norm, eps, resid;
         
             /* Test 1 */
             norm_b = dlange_("1", &m, &nrhs, B, &nrhs, work);
             eps = dlamch_("P");
         
             /* Compute AX-B */
             dgemm_(trans, "N", &m, &nrhs, &m, &d_one, A, &m, X, &ldx, &d_n_one, B, &ldb);
             norm = dlange_("1", &m, &nrhs, B, &nrhs, work);
         
             resid = norm / (eps * norm_b * (double)m);
         
             *residual = (double)resid;
             break;
         }
         case COMPLEX:
         {
             float norm_b, norm, eps, resid;
         
             /* Test 1 */
             norm_b = clange_("1", &m, &nrhs, B, &nrhs, work);
             eps = slamch_("P");
         
             /* Compute AX-B */
             cgemm_(trans, "N", &m, &nrhs, &m, &c_one, A, &m, X, &ldx, &c_n_one, B, &ldb);
             norm = clange_("1", &m, &nrhs, B, &nrhs, work);
         
             resid = norm / (eps * norm_b * (float)m);
         
             *residual = (double)resid;
             break;
         }
         case DOUBLE_COMPLEX:
         {
             double norm_b, norm, eps, resid;
         
             /* Test 1 */
             norm_b = zlange_("1", &m, &nrhs, B, &nrhs, work);
             eps = dlamch_("P");
         
             /* Compute AX-B */
             zgemm_(trans, "N", &m, &nrhs, &m, &z_one, A, &m, X, &ldx, &z_n_one, B, &ldb);
             norm = zlange_("1", &m, &nrhs, B, &nrhs, work);
         
             resid = norm / (eps * norm_b * (double)m);
         
             *residual = (double)resid;
             break;
         }
    }
    free_vector(work);
}