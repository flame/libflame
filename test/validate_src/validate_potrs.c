/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_potrs.c
 *  @brief Defines validate function of POTRS() to use in test suite.
 *  */

#include "test_common.h"

void validate_potrs(char *uplo, integer m,
    void *A,
    void *A_test,
    integer datatype,
    void *x,
    void *b,
    double* residual)
{
    integer incx = 1, incy = 1;

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = snrm2_(&m, b, &incx);
            eps = slamch_("P");

            /* Compute Ax-b */
            sgemv_("N", &m, &m, &s_one, A, &m, x, &incx, &s_n_one, b, &incy);
            norm = snrm2_(&m, b, &incx);

            resid = norm/(eps * norm_b * (float)m);

            *residual = (double)resid;
            break;
        }
        case DOUBLE:
        {
            double norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = dnrm2_(&m, b, &incx);
            eps = dlamch_("P");

            /* Compute Ax-b */
            dgemv_("N", &m, &m, &d_one, A, &m, x, &incx, &d_n_one, b, &incy);
            norm = dnrm2_(&m, b, &incx);

            resid = norm/(eps * norm_b * (double)m);

            *residual = (double)resid;
            break;
        }
        case COMPLEX:
        {
            float norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = scnrm2_(&m, b, &incx);
            eps = slamch_("P");

            /* Compute Ax-b */
            cgemv_("N", &m, &m, &c_one, A, &m, x, &incx, &c_n_one, b, &incy);
            norm = scnrm2_(&m, b, &incx);

            resid = norm/(eps * norm_b * (float)m);

            *residual = (double)resid;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = dznrm2_(&m, b, &incx);
            eps = dlamch_("P");

            /* Compute Ax-b */
            zgemv_("N", &m, &m, &z_one, A, &m, x, &incx, &z_n_one, b, &incy);
            norm = dznrm2_(&m, b, &incx);

            resid = norm/(eps * norm_b * (double)m);

            *residual = (double)resid;
            break;
        }
    }
}