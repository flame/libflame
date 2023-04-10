/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_spffrtx.c
 *  @brief Defines validate function of SPFFRTX() to use in test suite.
 *  */

#include "test_common.h"

void validate_spffrtx(integer n,
    integer ncolm,
    void* A,
    void* AP,
    integer datatype,
    double* residual)
{
    if(n == 0 || ncolm == 0)
        return;
    integer i, j, di;
    void* work = NULL;
    void* L = NULL, *D = NULL, *T = NULL;


    create_matrix(datatype, &L, n, n);
    create_matrix(datatype, &D, n, n);
    create_matrix(datatype, &T, n, n);
    /* Creating work buffer */
    create_vector(datatype, &work, 2*n);
    set_identity_matrix(datatype, n, n, L, n);
    set_identity_matrix(datatype, n, n, D, n);
    reset_matrix(datatype, n, n, T, n);

    switch (datatype)
    {
        case FLOAT:
        {
            float norma, norm, eps, resid;
            /* Test 1 */
            /* Unpack L, D and L' */
            di = 0;
            for (i = 0; i < ncolm; i++)
            {
                /* Form L */
                for (j = i+1; j < n; j++)
                {
                    ((float*)L)[i * n + j] = ((float*)AP)[j - i + di];
                }

                ((float*)D)[i * n + i] = ((float*)AP)[di];

                di += n - i;
            }

            for (; i < n; i++)
            {
                for (j = i; j < n; j++)
                {
                    ((float*)D)[j + i * n] = ((float*)AP)[di + j - i];
                    ((float*)D)[i + j * n] = ((float*)AP)[di + j - i];
                }

                di += n - i;
            }

            norma = slange_("1", &n, &n, A, &n, work);
            eps = slamch_("P");
            /* Compute L * D */ 
            sgemm_("N", "N", &n, &n, &n, &s_one, L, &n, D, &n, &s_zero, T, &n);
            /* Compute (L * D * L') - A */ 
            sgemm_("N", "T", &n, &n, &n, &s_n_one, T, &n, L, &n, &s_one, A, &n);
            norm = slange_("1", &n, &n, A, &n, work);
            resid = norm / (eps * norma * (float)n);        
            *residual = (double)resid;
            break;
        }
        case DOUBLE:
        {
            double norma, norm, eps, resid;
            /* Test 1 */
            /* Unpack L, D and L' */
            di = 0;
            for (i = 0; i < ncolm; i++)
            {
                /* Form L */
                for (j = i+1; j < n; j++)
                {
                    ((double*)L)[i * n + j] = ((double*)AP)[j - i + di];
                }

                ((double*)D)[i * n + i] = ((double*)AP)[di];

                di += n - i;
            }

            for (; i < n; i++)
            {
                for (j = i; j < n; j++)
                {
                    ((double*)D)[j + i * n] = ((double*)AP)[di + j - i];
                    ((double*)D)[i + j * n] = ((double*)AP)[di + j - i];
                }

                di += n - i;
            }
            norma = dlange_("1", &n, &n, A, &n, work);
            eps = dlamch_("P");
            /* Compute L * D */ 
            dgemm_("N", "N", &n, &n, &n, &d_one, L, &n, D, &n, &d_zero, T, &n);
            /* Compute (L * D * L') - A */ 
            dgemm_("N", "T", &n, &n, &n, &d_one, T, &n, L, &n, &d_n_one, A, &n);
            norm = dlange_("1", &n, &n, A, &n, work);
            resid = norm / (eps * norma * (double)n);
            *residual = (double)resid;
            break;
        }
        case COMPLEX:
        {
            float norma, norm, eps, resid;
            /* Test 1 */
            /* Unpack L, D and L' */
            di = 0;
            for (i = 0; i < ncolm; i++)
            {
                /* Form L */
                for (j = i+1; j < n; j++)
                {
                    ((scomplex*)L)[i * n + j].real = ((scomplex*)AP)[j - i + di].real;
                    ((scomplex*)L)[i * n + j].imag = ((scomplex*)AP)[j - i + di].imag;
                }

                ((scomplex*)D)[i * n + i].real = ((scomplex*)AP)[di].real;
                ((scomplex*)D)[i * n + i].imag = ((scomplex*)AP)[di].imag;

                di += n - i;
            }

            for (; i < n; i++)
            {
                for (j = i; j < n; j++)
                {
                    ((scomplex*)D)[j + i * n].real = ((scomplex*)AP)[di + j - i].real;
                    ((scomplex*)D)[j + i * n].imag = ((scomplex*)AP)[di + j - i].imag;
                    ((scomplex*)D)[i + j * n].real = ((scomplex*)AP)[di + j - i].real;
                    ((scomplex*)D)[i + j * n].imag = ((scomplex*)AP)[di + j - i].imag; 
                }

                di += n - i;
            }
            norma = clange_("1", &n, &n, A, &n, work);
            eps = slamch_("P");
            /* Compute L * D */ 
            cgemm_("N", "N", &n, &n, &n, &c_one, L, &n, D, &n, &c_zero, T, &n);
            /* Compute (L * D * L') - A */ 
            cgemm_("N", "T", &n, &n, &n, &c_one, T, &n, L, &n, &c_n_one, A, &n);
            norm = clange_("1", &n, &n, A, &n, work);
            resid = norm / (eps * norma * (float)n);
            *residual = (double)resid;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norma, norm, eps, resid;
            /* Test 1 */
            /* Unpack L, D and L' */
            di = 0;
            for (i = 0; i < ncolm; i++)
            {
                /* Form L */
                for (j = i+1; j < n; j++)
                {
                    ((dcomplex*)L)[i * n + j].real = ((dcomplex*)AP)[j - i + di].real;
                    ((dcomplex*)L)[i * n + j].imag = ((dcomplex*)AP)[j - i + di].imag;
                }

                ((dcomplex*)D)[i * n + i].real = ((dcomplex*)AP)[di].real;
                ((dcomplex*)D)[i * n + i].imag = ((dcomplex*)AP)[di].imag;

                di += n - i;
            }

            for (; i < n; i++)
            {
                for (j = i; j < n; j++)
                {
                    ((dcomplex*)D)[j + i * n].real = ((dcomplex*)AP)[di + j - i].real;
                    ((dcomplex*)D)[j + i * n].imag = ((dcomplex*)AP)[di + j - i].imag;
                    ((dcomplex*)D)[i + j * n].real = ((dcomplex*)AP)[di + j - i].real;
                    ((dcomplex*)D)[i + j * n].imag = ((dcomplex*)AP)[di + j - i].imag;
                }

                di += n - i;
            }
            norma = zlange_("1", &n, &n, A, &n, work);
            eps = dlamch_("P");
            /* Compute L * D */ 
            zgemm_("N", "N", &n, &n, &n, &z_one, L, &n, D, &n, &z_zero, T, &n);
            /* Compute (L * D * L') -A */ 
            zgemm_("N", "T", &n, &n, &n, &z_one, T, &n, L, &n, &z_n_one, A, &n);
            norm = zlange_("1", &n, &n, A, &n, work);
            resid = norm / (eps * norma * (double)n);
            *residual = (double)resid;
            break;
        }
    }

    *residual = *residual / 10.0;

    free_vector(work);
    free_matrix(L);
    free_matrix(D);
    free_matrix(T);
}

