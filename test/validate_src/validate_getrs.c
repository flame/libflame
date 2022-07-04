/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_getrs.c
 *  @brief Defines validate function of GETRS() to use in test suite.
 *  */

#include "test_common.h"

void validate_getrs(char *trans,
    integer m,
    integer n,
    void* A,
    void* B,
    void* X,
    integer datatype,
    double* residual)
{
    switch (datatype)
    {
        case FLOAT:
        {
            float norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = snrm2_(&m, B, &i_one);
            eps = slamch_("P");

            /* Compute Ax-b */
            sgemv_(trans, &m, &m, &s_one, A, &m, X, &i_one, &s_n_one, B, &i_one);
            norm = snrm2_(&m, B, &i_one);

            resid = norm / (eps * norm_b * (float)m);

            *residual = (double)resid;
            break;
        }
        case DOUBLE:
        {
            double norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = dnrm2_(&m, B, &i_one);
            eps = dlamch_("P");

            /* Compute Ax-b */
            dgemv_(trans, &m, &m, &d_one, A, &m, X, &i_one, &d_n_one, B, &i_one);
            norm = dnrm2_(&m, B, &i_one);

            resid = norm / (eps * norm_b * (double)m);

            *residual = (double)resid;
            break;
        }
        case COMPLEX:
        {
            float norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = scnrm2_(&m, B, &i_one);
            eps = slamch_("P");

            /* Compute Ax-b */
            cgemv_(trans, &m, &m, &c_one, A, &m, X, &i_one, &c_n_one, B, &i_one);
            norm = scnrm2_(&m, B, &i_one);

            resid = norm / (eps * norm_b * (float)m);

            *residual = (double)resid;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = dznrm2_(&m, B, &i_one);
            eps = dlamch_("P");

            /* Compute Ax-b */
            zgemv_(trans, &m, &m, &z_one, A, &m, X, &i_one, &z_n_one, B, &i_one);
            norm = dznrm2_(&m, B, &i_one);

            resid = norm / (eps * norm_b * (double)m);

            *residual = (double)resid;
            break;
        }
    }
}