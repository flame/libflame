/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_rot.c
 *  @brief Defines validate function of ROT() to use in test suite.
 *  */

#include "test_common.h"
void validate_lartg(integer datatype, void *f, void *g, void *r, void *c, void *s, double* residual)
{
    void *out_zero = NULL;
    create_vector(datatype, &out_zero, 1);
    reset_vector(datatype, out_zero, 1, 1);

    switch(datatype)
    {
        case FLOAT:
        {
            float resid1 = 0.0;
            float resid2 = 0.0;
            float eps, norm_r, norm_f, norm_res, norm_1;
            float res = 0.0;

            eps = fla_lapack_slamch("P");

            /*Test 1 validating C and S */
            norm_1 = snrm2_(&i_one, &s_one, &i_one);
            /*square root of ((c*c) + (s*s)) = 1*/
            res = sqrt((((float *)c)[0]*((float *)c)[0]) + (((float *)s)[0]*((float *)s)[0]));
            /*res->res-1*/
            saxpy_(&i_one, &s_n_one, &s_one, &i_one, &res, &i_one);
            norm_res = snrm2_(&i_one, &res, &i_one);
            resid1 = (norm_res/ norm_1/ eps);

            /*Test 2 Validating R*/
            norm_r = snrm2_(&i_one, r, &i_one);
            /*Hermitian Transpose of original rotation vector*/
            ((float *)s)[0] = -((float *)s)[0];
            /*Rotating output vectors*/
            fla_lapack_srot(&i_one, r, &i_one, out_zero, &i_one,((float *)c), ((float *)s));

            /*r-f*/
            saxpy_(&i_one, &s_n_one, f, &i_one, r, &i_one);
            saxpy_(&i_one, &s_n_one, g, &i_one, out_zero, &i_one);

            norm_f = snrm2_(&i_one, r, &i_one);
            resid2 = (norm_f/ norm_r/ eps);
    
            *residual = (double)fla_max(resid1, resid2);
            break;
        }

        case DOUBLE:
        {
            double resid1 = 0.0;
            double resid2 = 0.0;
            double eps, norm_r, norm_f, norm_res, norm_1;
            double res = 0.0;

            eps = fla_lapack_slamch("P");

            /*Test 1 validating C and S */
            norm_1 = dnrm2_(&i_one, &d_one, &i_one);
            /*square root of ((c*c) + (s*s)) = 1*/
            res = sqrt((((double *)c)[0]*((double *)c)[0]) + (((double *)s)[0]*((double *)s)[0]));
            /*res->res -1*/
            daxpy_(&i_one, &d_n_one, &d_one, &i_one, &res, &i_one);
            norm_res = dnrm2_(&i_one, &res, &i_one);
            resid1 = (norm_res/ norm_1/ eps);

            /*Test 2 Validating R*/
            norm_r = dnrm2_(&i_one, r, &i_one);
            /*Hermitian Transpose of original rotation vector*/
            ((double *)s)[0] = -((double *)s)[0];
            /*Rotating output vectors*/
            fla_lapack_drot(&i_one, r, &i_one, out_zero, &i_one,((double *)c), ((double *)s));

            /*r-f*/
            daxpy_(&i_one, &d_n_one, f, &i_one, r, &i_one);
            daxpy_(&i_one, &d_n_one, g, &i_one, out_zero, &i_one);
            norm_f = dnrm2_(&i_one, r, &i_one);
            resid2 = (norm_f/ norm_r/ eps);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }

        case COMPLEX:
        {
            float resid1 = 0.0;
            float resid2 = 0.0;
            float eps, norm_r, norm_f, norm_res, norm_1;
            float res = 0.0;;

            eps = fla_lapack_slamch("P");

            /*Test 1 validating C and S */
            norm_1 = snrm2_(&i_one, &c_one, &i_one);
            /*square root of ((c*c) + (s*s)) = 1*/
            res = sqrt((((float *)c)[0]*((float *)c)[0]) + (((scomplex *)s)[0].real*((scomplex *)s)[0].real) + 
                       (((scomplex *)s)[0].imag*((scomplex *)s)[0].imag));
            /*res->res - 1*/
            saxpy_(&i_one, &s_n_one, &s_one, &i_one, &res, &i_one);
            norm_res = snrm2_(&i_one, &res, &i_one);
            resid1 = (norm_res/ norm_1/ eps);

            /*Test 2 Validating R*/
            norm_r = scnrm2_(&i_one, r, &i_one);
            /*Hermitian Transpose of original rotation vector*/
            ((scomplex *)s)[0].real = -((scomplex *)s)[0].real;
            ((scomplex *)s)[0].imag = -((scomplex *)s)[0].imag;
            /*Rotating output vectors*/
            fla_lapack_crot(&i_one, r, &i_one, out_zero, &i_one,((float *)c), ((scomplex *)s));

            /*r-f*/
            caxpy_(&i_one, &c_n_one, f, &i_one, r, &i_one);
            caxpy_(&i_one, &c_n_one, g, &i_one, out_zero, &i_one);
            norm_f = scnrm2_(&i_one, r, &i_one);
            resid2 = (norm_f/ norm_r/ eps);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double resid1 = 0.0;
            double resid2 = 0.0;
            double eps, norm_r, norm_f, norm_res, norm_1;
            double res = 0.0;

            eps = fla_lapack_slamch("P");

            /*Test 1 validating C and S */
            norm_1 = snrm2_(&i_one, &c_one, &i_one);
            /*square root of ((c*c) + (s*s)) = 1*/
            res = sqrt((((double *)c)[0]*((double *)c)[0]) + (((dcomplex *)s)[0].real*((dcomplex *)s)[0].real) + 
                       (((dcomplex *)s)[0].imag*((dcomplex *)s)[0].imag));
            /*res-> res - 1*/
            daxpy_(&i_one, &d_n_one, &d_one, &i_one, &res, &i_one);
            norm_res = dnrm2_(&i_one, &res, &i_one);
            resid1 = (norm_res/ norm_1/ eps);

            /*Test 2 Validating R*/
            norm_r = dznrm2_(&i_one, r, &i_one);
            /*Hermitian Transpose of original rotation vector*/
            ((dcomplex *)s)[0].real = -((dcomplex *)s)[0].real;
            ((dcomplex *)s)[0].imag = -((dcomplex *)s)[0].imag;
            /*Rotating output vectors*/
            fla_lapack_zrot(&i_one, r, &i_one, out_zero, &i_one,((double *)c), ((dcomplex *)s));

            /*r-f*/
            zaxpy_(&i_one, &z_n_one, f, &i_one, r, &i_one);
            zaxpy_(&i_one, &z_n_one, g, &i_one, out_zero, &i_one);
            resid2 = (norm_f/ norm_r/ eps);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }
    free_vector(out_zero);
}