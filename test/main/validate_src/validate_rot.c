/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_rot.c
 *  @brief Defines validate function of ROT() to use in test suite.
 *  */

#include "test_common.h"
void validate_rot(integer datatype, integer n, void *cx, void *cx_test, integer incx, void *cy, void *cy_test, integer incy, void *c, void *s, double* residual)
{
    switch(datatype)
    {
        case FLOAT:
        {
            float resid1 = 0.0;
            float resid2 = 0.0;
            float eps = fla_lapack_slamch("P");
            float norm_resid1, norm_resid2, norm_cx, norm_cy;
            
            norm_cx = snrm2_(&n, cx, &incx);
            norm_cy = snrm2_(&n, cy, &incy);
            /*Hermitian Transpose of original rotation vector*/
            ((float *)s)[0] = -((float *)s)[0];

            /*Rotating output vectors*/
            fla_lapack_srot(&n, cx, &incx, cy, &incy,((float *)c), ((float *)s));

            /*cx-cx_test*/
            saxpy_(&n, &s_n_one, cx_test, &incx, cx, &incx);
            saxpy_(&n, &s_n_one, cy_test, &incy, cy, &incy);

            norm_resid1 = fla_max(resid1, snrm2_(&n, cx, &incx));
            resid1 = (norm_resid1)/ eps/ norm_cx/ ((float)n);
            norm_resid2 = fla_max(resid2, snrm2_(&n, cy, &incy));
            resid2 = (norm_resid2)/ eps/ norm_cy/ ((float)n);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }

        case DOUBLE:
        {
            double resid1 = 0.f;
            double resid2 = 0.f;
            double eps = fla_lapack_slamch("P");
            double norm_resid1, norm_resid2, norm_cx, norm_cy;

            norm_cx = dnrm2_(&n, cx, &incx);
            norm_cy = dnrm2_(&n, cy, &incy);
            /*Hermitian Transpose of original rotation vector*/
            ((double *)s)[0] = -((double *)s)[0];

            /*Rotating output vectors*/
            fla_lapack_drot(&n, cx, &incx, cy, &incy,((double *)c), ((double *)s));

            /*cx-cx_test*/
            daxpy_(&n, &d_n_one, cx_test, &incx, cx, &incx);
            daxpy_(&n, &d_n_one, cy_test, &incy, cy, &incy);

            norm_resid1 = fla_max(resid1, dnrm2_(&n, cx, &incx));
            resid1 = (norm_resid1)/ eps/ norm_cx/ ((double)n);
            norm_resid2 = fla_max(resid2, dnrm2_(&n, cy, &incy));
            resid2 = (norm_resid2)/ eps/ norm_cy/ ((double)n);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            float resid1 = 0.0;
            float resid2 = 0.0;
            float eps = fla_lapack_slamch("P");
            float norm_resid1, norm_resid2, norm_cx, norm_cy;

            norm_cx = scnrm2_(&n, cx, &incx);
            norm_cy = scnrm2_(&n, cy, &incy);

            /*Hermitian Transpose of original rotation vector*/
            ((scomplex *)s)[0].real = -((scomplex *)s)[0].real;
            ((scomplex *)s)[0].imag = -((scomplex *)s)[0].imag;

            /*Rotating output vectors*/
            fla_lapack_crot(&n, cx, &incx, cy, &incy,((float *)c), ((scomplex *)s));

            /*cx-cx_test*/
            caxpy_(&n, &c_n_one, cx_test, &incx, cx, &incx);
            caxpy_(&n, &c_n_one, cy_test, &incy, cy, &incy);

            norm_resid1 = fla_max(resid1, scnrm2_(&n, cx, &incx));
            resid1 = (norm_resid1)/ eps/ norm_cx/ ((float)n);
            norm_resid2 = fla_max(resid2, scnrm2_(&n, cy, &incy));
            resid2 = (norm_resid2)/ eps/ norm_cy/ ((float)n);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double resid1 = 0.f;
            double resid2 = 0.f;
            double eps = fla_lapack_slamch("P");
            double norm_resid1, norm_resid2, norm_cx, norm_cy;

            norm_cx = dnrm2_(&n, cx, &incx);
            norm_cy = dnrm2_(&n, cy, &incy);

            /*Hermitian Transpose of original rotation vector*/
            ((dcomplex *)s)[0].real = -((dcomplex *)s)[0].real;
            ((dcomplex *)s)[0].imag = -((dcomplex *)s)[0].imag;

            /*Rotating output vectors*/
            fla_lapack_zrot(&n, cx, &incx, cy, &incy,((double *)c), ((dcomplex *)s));

            /*cx-cx_test*/
            zaxpy_(&n, &z_n_one, cx_test, &incx, cx, &incx);
            zaxpy_(&n, &z_n_one, cy_test, &incy, cy, &incy);

            norm_resid1 = fla_max(resid1, dznrm2_(&n, cx, &incx));
            resid1 = (norm_resid1)/ eps/ norm_cx/ ((double)n);
            norm_resid2 = fla_max(resid2, dznrm2_(&n, cy, &incy));
            resid2 = (norm_resid2)/ eps/ norm_cy/ ((double)n);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }
}