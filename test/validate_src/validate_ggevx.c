/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_ggevx.c
 *  @brief Defines validate function of GGEVX() to use in test suite.
 *  */

#include "test_common.h"

void validate_ggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer n, void* A, integer lda, void* B, integer ldb, void* alpha, void* alphar, void* alphai, void* beta, void* VL, integer ldvl, void* VR, integer ldvr, integer datatype, double *residual)
{
    integer i, j;
    void *work = NULL;

    switch (datatype)
    {
        case FLOAT:
        {
            float eps, norm_A, alphar_t, max_val = 0.0;
            void *YC = NULL, * Y = NULL, * VRTemp = NULL, * VLTemp = NULL, *VRC = NULL, * VLC = NULL;
            scomplex alphac;

            create_vector(datatype, &Y, n);
            create_vector(COMPLEX, &YC, n);
            eps = slamch_("P");
            norm_A = slange_("1", &n, &n, A, &lda, work);
            /* Test 1 */
            /* Validation for 'V' 'V' combination */
            if (*jobvr == 'V')
            {
                create_vector(COMPLEX, &VRC, n);
                create_vector(datatype, &VRTemp, n);

                for (i = 0; i < n; i++)
                {
                    if (((float*)alphai)[i] != 0.0)
                    {
                        reset_matrix(COMPLEX, n, 1, YC, lda);
                        for (j = 0; j < n; j++)
                        {
                            ((scomplex*)VRC)[j].real = ((float*)VR)[i * ldvr + j];
                            ((scomplex*)VRC)[j].imag = ((float*)VR)[(i + 1) * ldvr + j];
                        }
                        alphac.real = ((float*)beta)[i];
                        alphac.imag = 0;
                        /* beta * A * VR */
                        scgemv('N',0, n, n, &alphac, A, lda, VRC, i_one, s_zero, YC, i_one);
                        alphac.real = ((float*)alphar)[i];
                        alphac.imag = ((float*)alphai)[i];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        scgemv('N',0, n, n, &alphac, B, ldb, VRC, i_one, s_n_one, YC, i_one);
                        max_val = max(max_val, scnrm2_(&n, YC, &i_one));

                        reset_matrix(COMPLEX, n, 1, YC, lda);
                        alphac.real = ((float*)beta)[i + 1];
                        alphac.imag = 0;
                        for (j = 0; j < n; j++)
                        {
                            ((scomplex*)VRC)[j].real = ((float*)VR)[i * ldvr + j];
                            ((scomplex*)VRC)[j].imag = -((float*)VR)[(i + 1) * ldvr + j];
                        }
                        /* beta * A * VR */
                        scgemv('N', 0, n, n, &alphac, A, lda, VRC, i_one, s_zero, YC, i_one);
                        alphac.real = ((float*)alphar)[i + 1];
                        alphac.imag = ((float*)alphai)[i + 1];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        scgemv('N', 0, n, n, &alphac, B, ldb, VRC, i_one, s_n_one, YC, i_one);
                        i++;
                        max_val = max(max_val, scnrm2_(&n, YC, &i_one));
                    }
                    else
                    {
                        reset_matrix(FLOAT, n, 1, Y, lda);
                        for (j = 0; j < n; j++)
                        {
                            ((float*)VRTemp)[j] = ((float*)VR)[i * ldvr + j];
                        }
                        /* beta*A*vr(i) */
                        alphar_t = ((float*)beta)[i];
                        sgemv_("N", &n, &n, &alphar_t, A, &lda, VRTemp, &i_one, &s_zero, Y, &i_one);
                        /* alpha*B*vr(i) - beta*A*vr(i) */
                        alphar_t = ((float*)alphar)[i];
                        sgemv_("N", &n, &n, &alphar_t, B, &ldb, VRTemp, &i_one, &s_n_one, Y, &i_one);
                        max_val = max(max_val, snrm2_(&n, Y, &i_one));
                    }
                }
                free_vector(VRC);
                free_vector(VRTemp);

            }
            if (*jobvl == 'V')
            {
                create_vector(COMPLEX, &VLC, n);
                create_vector(FLOAT, &VLTemp, n);

                for (i = 0; i < n; i++)
                {
                    if (((float*)alphai)[i] != 0.0)
                    {
                        alphac.real = ((float*)beta)[i];
                        alphac.imag = 0;
                        for (j = 0; j < n; j++)
                        {
                            ((scomplex*)VLC)[j].real = ((float*)VL)[i * ldvl + j];
                            ((scomplex*)VLC)[j].imag = ((float*)VL)[(i + 1) * ldvl + j];
                        }
                        reset_matrix(COMPLEX, n, 1, YC, lda);
                        /* beta * A * VL */
                        scgemv('T',0, n, n, &alphac, A, lda, VLC, i_one, s_zero, YC, i_one);
                        alphac.real = ((float*)alphar)[i];
                        alphac.imag = -((float*)alphai)[i];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        scgemv('T',0, n, n, &alphac, B, ldb, VLC, i_one, s_n_one, YC, i_one);
                        max_val = max(max_val, scnrm2_(&n, YC, &i_one));

                        reset_matrix(COMPLEX, n, 1, YC, lda);
                        alphac.real = ((float*)beta)[i + 1];
                        alphac.imag = 0;
                        for (j = 0; j < n; j++)
                        {
                            ((scomplex*)VLC)[j].real = ((float*)VL)[i * ldvl + j];
                            ((scomplex*)VLC)[j].imag = -((float*)VL)[(i + 1) * ldvl + j];
                        }
                        /* beta * A * VR */
                        scgemv('T', 0, n, n, &alphac, A, lda, VLC, i_one, s_zero, YC, i_one);
                        alphac.real = ((float*)alphar)[i + 1];
                        alphac.imag = -((float*)alphai)[i + 1];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        scgemv('T', 0, n, n, &alphac, B, ldb, VLC, i_one, s_n_one, YC, i_one);
                        i++;
                        max_val = max(max_val, scnrm2_(&n, YC, &i_one));
                    }
                    else
                    {
                        reset_matrix(FLOAT, n, 1, Y, lda);
                        for (j = 0; j < n; j++)
                        {
                            ((float*)VLTemp)[j] = ((float*)VL)[i * ldvl + j];
                        }
                        /* beta*A*vl(i) */
                        alphar_t = ((float*)beta)[i];
                        sgemv_("C", &n, &n, &alphar_t, A, &lda, VLTemp, &i_one, &s_zero, Y, &i_one);
                        /* alpha*B*vl**H(i) - beta*A*vl(i) */
                        alphar_t = ((float*)alphar)[i];
                        sgemv_("C", &n, &n, &alphar_t, B, &ldb, VLTemp, &i_one, &s_n_one, Y, &i_one);
                        max_val = max(max_val, snrm2_(&n, Y, &i_one));
                    }
                }
                free_vector(VLC);
                free_vector(VLTemp);
            }

            *residual = (double)max_val/(eps * norm_A * (float)n);
            free_vector(Y);
            free_vector(YC);
            break;
        }
        case DOUBLE:
        {
            double norm_A, eps, alphar_t, max_val = 0.0;
            void *YC = NULL, *Y=NULL, *VRTemp = NULL, *VLTemp = NULL, *VRC = NULL, *VLC = NULL;
            dcomplex alphac;

            create_vector(datatype, &Y, n);
            create_vector(DOUBLE_COMPLEX, &YC, n);
            eps = dlamch_("P");
            norm_A = dlange_("1", &n, &n, A, &lda, work);
            /* Test 1 */
            /* Validation for 'V' 'V' combination */
            if (*jobvr == 'V')
            {
                create_vector(DOUBLE_COMPLEX, &VRC, n);
                create_vector(datatype, &VRTemp, n);

                for (i = 0; i < n; i++)
                {
                    if (((double*)alphai)[i] != 0.0)
                    {
                        alphac.real = ((double*)beta)[i];
                        alphac.imag = 0;
                        for (j = 0; j < n; j++)
                        {
                            ((dcomplex*)VRC)[j].real = ((double*)VR)[i * ldvr + j];
                            ((dcomplex*)VRC)[j].imag = ((double*)VR)[(i + 1) * ldvr + j];
                        }
                        reset_matrix(DOUBLE_COMPLEX, n, 1, YC, lda); 
                        /* beta * A * VR */ 
                        dzgemm_("N", "N", &n, &i_one, &n, &alphac, A, &lda, VRC, &ldvr, &z_zero, YC, &n);
                        alphac.real = ((double*)alphar)[i];
                        alphac.imag = ((double*)alphai)[i];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        dzgemm_("N", "N", &n, &i_one, &n, &alphac, B, &ldb, VRC, &ldvr, &z_n_one, YC, &n);
                        max_val = max(max_val, dznrm2_(&n, YC, &i_one));

                        reset_matrix(DOUBLE_COMPLEX, n, 1, YC, lda);
                        alphac.real = ((double*)beta)[i + 1];
                        alphac.imag = 0;
                        for (j = 0; j < n; j++)
                        {
                            ((dcomplex*)VRC)[j].real = ((double*)VR)[i * ldvr + j];
                            ((dcomplex*)VRC)[j].imag = -((double*)VR)[(i + 1) * ldvr + j];
                        }
                        /* beta * A * VR */
                        dzgemm_("N", "N", &n, &i_one, &n, &alphac, A, &lda, VRC, &ldvr, &z_zero, YC, &n);
                        alphac.real = ((double*)alphar)[i + 1];
                        alphac.imag = ((double*)alphai)[i + 1];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        dzgemm_("N", "N", &n, &i_one, &n, &alphac, B, &ldb, VRC, &ldvr, &z_n_one, YC, &n);
                        i++;
                        max_val = max(max_val, dznrm2_(&n, YC, &i_one));
                    }
                    else
                    {
                        reset_matrix(datatype, n, 1, Y, lda);
                        for (j = 0; j < n; j++)
                        {
                            ((double*)VRTemp)[j]= ((double*)VR)[i * ldvr + j];
                        }
                        /* beta*A*vr(i) */
                        alphar_t = ((double*)beta)[i];
                        dgemv_("N", &n, &n, &alphar_t, A, &lda, VRTemp, &i_one, &d_zero, Y, &i_one);
                        /* alpha*B*vr(i) - beta*A*vr(i) */
                        alphar_t = ((double*)alphar)[i];
                        dgemv_("N", &n, &n, &alphar_t, B, &ldb, VRTemp, &i_one, &d_n_one, Y, &i_one);
                        max_val = max(max_val, dnrm2_(&n, Y, &i_one));
                    }
                }
                free_vector(VRC);
                free_vector(VRTemp);

            }
            if (*jobvl == 'V')
            {
                create_vector(DOUBLE_COMPLEX, &VLC, n);
                create_vector(datatype, &VLTemp, n);

                for (i = 0; i < n; i++)
                {
                    if (((double*)alphai)[i] != 0.0)
                    {
                        alphac.real = ((double*)beta)[i];
                        alphac.imag = 0;
                        for (j = 0; j < n; j++)
                        {
                            ((dcomplex*)VLC)[j].real = ((double*)VL)[i * ldvl + j];
                            ((dcomplex*)VLC)[j].imag = ((double*)VL)[(i + 1) * ldvl + j];
                        }
                        reset_matrix(DOUBLE_COMPLEX, n, 1, YC, lda);
                        /* beta * VL**H * A --> beta A**H * VL */

                        dzgemm_("C", "N", &n, &i_one, &n, &alphac, A, &lda, VLC, &ldvl, &z_zero, YC, &n);
                        alphac.real = ((double*)alphar)[i];
                        alphac.imag = -((double*)alphai)[i];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        dzgemm_("C", "N", &n, &i_one, &n, &alphac, B, &ldb, VLC, &ldvl, &z_n_one, YC, &n);
                        max_val = max(max_val, dznrm2_(&n, YC, &i_one));

                        reset_matrix(DOUBLE_COMPLEX, n, 1, YC, lda);
                        alphac.real = ((double*)beta)[i + 1];
                        alphac.imag = 0;
                        for (j = 0; j < n; j++)
                        {
                            ((dcomplex*)VLC)[j].real = ((double*)VL)[i * ldvl + j];
                            ((dcomplex*)VLC)[j].imag = -((double*)VL)[(i + 1) * ldvl + j];
                        }
                        /* beta * A * VR */
                        dzgemm_("C", "N", &n, &i_one, &n, &alphac, A, &lda, VLC, &ldvl, &z_zero, YC, &n);
                        alphac.real = ((double*)alphar)[i + 1];
                        alphac.imag = -((double*)alphai)[i + 1];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        dzgemm_("C", "N", &n, &i_one, &n, &alphac, B, &ldb, VLC, &ldvl, &z_n_one, YC, &n);
                        i++;
                        max_val = max(max_val, dznrm2_(&n, YC, &i_one));
                    }
                    else
                    {
                        reset_matrix(datatype, n, 1, Y, lda);
                        for (j = 0; j < n; j++)
                        {
                            ((double*)VLTemp)[j] = ((double*)VL)[i * ldvl + j];
                        }
                        /* beta*A*vl(i) */ 
                        alphar_t = ((double*)beta)[i];
                        dgemv_("C", &n, &n, &alphar_t, A, &lda, VLTemp, &i_one, &d_zero, Y, &i_one);
                        /* alpha*B*vl(i) - beta*A*vl(i) */
                        alphar_t = ((double*)alphar)[i];
                        dgemv_("C", &n, &n, &alphar_t, B, &ldb, VLTemp, &i_one, &d_n_one, Y, &i_one);
                        max_val = max(max_val, dnrm2_(&n, Y, &i_one));
                    }
                }
                free_vector(VLC);
                free_vector(VLTemp);
            }

            *residual = (double)max_val/(eps * norm_A * (double)n);
            free_vector(Y);
            free_vector(YC);
            break;
        }
        case COMPLEX:
        {
            float eps, norm_A, max_val=0.0;
            void* VRTemp = NULL, * VLTemp = NULL;
            void *Y = NULL;
            scomplex alphar_t;
            eps = slamch_("P");
            norm_A = clange_("1", &n, &n, A, &n, work);
            /* Test 1 */
            /* Validation for 'V' 'V' combination */
            create_vector(datatype, &VRTemp, n);
            create_vector(datatype, &VLTemp, n);
            create_vector(datatype, &Y, n);

            if (*jobvr == 'V')
            {
                for (i = 0; i < n; i++)
                {
                    reset_matrix(datatype, n, 1, Y, lda);
                    for (j = 0; j < n; j++)
                    {
                        ((scomplex*)VRTemp)[j] = ((scomplex*)VR)[i * ldvr + j];
                    }
                    /* beta * A * VR */
                    alphar_t = ((scomplex*)beta)[i];
                    cgemv_("N", &n, &n, &alphar_t, A, &lda, VRTemp, &i_one, &c_zero, Y, &i_one);
                    /* alpha * B * VR - beta * A * VR */
                    alphar_t = ((scomplex*)alpha)[i];
                    cgemv_("N", &n, &n, &alphar_t, B, &ldb, VRTemp, &i_one, &c_n_one, Y, &i_one);
                    max_val = max(max_val, scnrm2_(&n, Y, &i_one));
                }
            }

            if (*jobvl == 'V')
            {

                for (i = 0; i < n; i++)
                {
                    reset_matrix(datatype, n, 1, Y, lda);
                    for (j = 0; j < n; j++)
                    {
                        ((scomplex*)VLTemp)[j] = ((scomplex*)VL)[i * ldvl + j];
                    }
                    /* betaC * a**H * vl */
                    ((scomplex*)beta)[i].imag = -((scomplex*)beta)[i].imag;
                    alphar_t = ((scomplex*)beta)[i];
                    cgemv_("C", &n, &n, &alphar_t, A, &lda, VLTemp, &i_one, &c_zero, Y, &i_one);
                    /* alphaC * b**H * vl - beta * a**H * vr */
                    ((scomplex*)alpha)[i].imag = -((scomplex*)alpha)[i].imag;
                    alphar_t = ((scomplex*)alpha)[i];
                    cgemv_("C", &n, &n, &alphar_t, B, &ldb, VLTemp, &i_one, &c_n_one, Y, &i_one);
                    max_val = max(max_val, scnrm2_(&n, Y, &i_one));
                }
            }

            *residual = (double)max_val/(eps * norm_A * (float)n);
            free_vector(VLTemp);
            free_vector(Y);
            free_vector(VRTemp);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps, norm_A, max_val=0.0;
            void* VRTemp = NULL, * VLTemp = NULL;
            void *Y = NULL;
            dcomplex alphar_t;
            eps = dlamch_("P");
            norm_A = zlange_("1", &n, &n, A, &n, work);
            /* Test 1 */
            /* Validation for 'V' 'V' combination */
            create_vector(datatype, &VRTemp, n);
            create_vector(datatype, &VLTemp, n);
            create_vector(datatype, &Y, n);
 
            if (*jobvr == 'V')
            {
                for (i = 0; i < n; i++)
                {

                    reset_matrix(datatype, n, 1, Y, lda);
                    for (j = 0; j < n; j++)
                    {
                        ((dcomplex*)VRTemp)[j] = ((dcomplex*)VR)[i * ldvr + j];
                    }
                    /* beta * A * VR */
                    alphar_t = ((dcomplex*)beta)[i];
                    zgemv_("NoTran", &n, &n, &alphar_t, A, &lda, VRTemp, &i_one, &z_zero, Y, &i_one);
                    /* alpha * B * VR - beta * A * VR */
                    alphar_t = ((dcomplex*)alpha)[i];
                    zgemv_("NoTran", &n, &n, &alphar_t, B, &ldb, VRTemp, &i_one, &z_n_one, Y, &i_one);
                    max_val = max(max_val, dznrm2_(&n, Y, &i_one));
                }
            }

            if (*jobvl == 'V')
            {

                for (i = 0; i < n; i++)
                {
                    reset_matrix(datatype, n, 1, Y, lda);
                    for (j = 0; j < n; j++)
                    {
                        ((dcomplex*)VLTemp)[j] = ((dcomplex*)VL)[i * ldvl + j];
                    }
                    /* betaC * a**H * vl */
                    ((dcomplex*)beta)[i].imag = -((dcomplex*)beta)[i].imag;
                    alphar_t = ((dcomplex*)beta)[i];
                    zgemv_("C", &n, &n, &alphar_t, A, &lda, VLTemp, &i_one, &z_zero, Y, &i_one);
                    /* alphaC * b**H * vl - beta * a**H * vr */
                    ((dcomplex*)alpha)[i].imag = -((dcomplex*)alpha)[i].imag;
                    alphar_t = ((dcomplex*)alpha)[i];
                    zgemv_("C", &n, &n, &alphar_t, B, &ldb, VLTemp, &i_one, &z_n_one, Y, &i_one);
                    max_val = max(max_val, dznrm2_(&n, Y, &i_one));
                }
            }
 
            *residual = (double)max_val/(eps * norm_A * (double)n);
            free_vector(VLTemp);
            free_vector(Y);
            free_vector(VRTemp);
            break;
        }
    }
}
