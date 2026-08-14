/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

static inline __m256 cgetrf_avx2_load4sc(const scomplex *v, aocl_int64_t off)
{
    __m128 lo = _mm_loadu_ps((const float *)&v[off].real);
    __m128 hi = _mm_loadu_ps((const float *)&v[off + 2].real);
    return _mm256_insertf128_ps(_mm256_castps128_ps256(lo), hi, 1);
}

static inline void cgetrf_avx2_store4sc(scomplex *v, aocl_int64_t off, __m256 val)
{
    _mm_storeu_ps((float *)&v[off].real, _mm256_castps256_ps128(val));
    _mm_storeu_ps((float *)&v[off + 2].real, _mm256_extractf128_ps(val, 1));
}

static inline __m512 cgetrf_avx512_load4sc(const scomplex *v, aocl_int64_t off)
{
    return _mm512_zextps256_ps512(cgetrf_avx2_load4sc(v, off));
}

static inline void cgetrf_avx512_store4sc(scomplex *v, aocl_int64_t off, __m512 val)
{
    cgetrf_avx2_store4sc(v, off, _mm512_castps512_ps256(val));
}

/*
 * LU with partial pivoting for tiny matrices
 *
 * All the computations are done inline without using
 * corresponding BLAS APIs to reduce function overheads.
 */
int fla_cgetrf_small_avx512(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                            aocl_int_t *ipiv, aocl_int64_t *info)
{
    aocl_int64_t mi, ni;
    aocl_int64_t i, j, i_1, i_2, i_3, i_4, i_5, i_6, i_7;
    real max_val, t_val, z_val;
    scomplex *acur, *apiv, *asrc;
    scomplex z__1;
    aocl_int64_t p_idx;
    aocl_int64_t min_m_n = fla_min(*m, *n);
    __m512 alpha_real, alpha_img, x_real[2], x_img[2];
    __m512 bv[2], bv_p[2], xv0[2], xv1[2], yv0[2], yv1[2];

    *info = 0;

    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }

    if(*info != 0)
    {
        return 0;
    }

    for(i = 0; i < min_m_n; i++)
    {
        mi = *m - i;
        ni = *n - i;

        acur = &a[i + *lda * i];

        // Find the pivot element
        max_val = 0;
        p_idx = i;
        for(i_1 = 0; i_1 < mi; i_1++)
        {
            t_val = f2c_abs(acur[i_1].real) + f2c_abs(acur[i_1].imag);
            if(t_val > max_val)
            {
                max_val = t_val;
                p_idx = i + i_1;
            }
        }

        apiv = a + p_idx;
        asrc = a + i;
        ipiv[i] = (aocl_int_t)(p_idx + 1);

        // Swap rows
        if(apiv[*lda * i].real != 0. || apiv[*lda * i].imag != 0.)
        {
            if(p_idx != i)
            {
                for(i_1 = 0; i_1 < *n; i_1++)
                {
                    i_2 = i_1 * *lda;
                    t_val = apiv[i_2].real;
                    z_val = apiv[i_2].imag;
                    apiv[i_2].real = asrc[i_2].real;
                    apiv[i_2].imag = asrc[i_2].imag;
                    asrc[i_2].real = t_val;
                    asrc[i_2].imag = z_val;
                }
            }

            /*----------------unblocked LU algorithm-------------------------

            A00 | a01    A02
            ----|-----------
            a10 | a11    a12
            A20 | a21    A22

            alpha = 1 / a11
            a21 := a21 * alpha
            A22 := A22 - a21 * a12
            ------------------------------------------------------------------*/

            // Calculate scalefactors (a21) & update trailing matrix
            {
                real piv_r = acur->real;
                real piv_i = acur->imag;
                real inv = 1.0f / (piv_r * piv_r + piv_i * piv_i);

                z__1.real = piv_r * inv;
                z__1.imag = -piv_i * inv;
            }

            // Load alpha from memory
            alpha_real = _mm512_set1_ps(z__1.real);
            alpha_img = _mm512_set1_ps(z__1.imag);

            // Updates 8 rows of trailing matrix per iteration
            for(i_1 = 1; i_1 < mi - 7; i_1 += 8)
            {
                /*-----------Trailing matrix update for LU factorisation-----------

                A00 | a01    A02
                ----|-----------
                a10 | alpha  x
                A20 | b      Y

                b := alpha * b
                Y := Y - b * x


                SIMD algorithm:

                alpha_real = aR1  aR1  aR1  aR1
                alpha_img  = aI1  aI1  aI1  aI1
                bv         = bR1  bI1  bR2  bI2
                bv_p       = bI1  bR1  bI2  bR2
                x_real     = xR1  xR1  xR1  xR1
                x_img      = xI1  xI1  xI1  xI1
                yv         = yR1  yI1  yR2  yI2

                step 1 => b := alpha * b
                    bv_p = alpha_img * bv_p
                    bv = ( alpha_real * bv ) - bv_p

                step 2 => Y := Y - b * x
                    xv = x_img * bv_p
                    xv = ( x_real * bv ) - xv
                    yv = yv - xv
                ------------------------------------------------------------------*/

                i_2 = i_1 + 4;

                // Load alpha from memory
                bv[0] = cgetrf_avx512_load4sc(acur, i_1);
                bv[1] = cgetrf_avx512_load4sc(acur, i_2);
                bv_p[0] = _mm512_permute_ps(bv[0], 0xB1);
                bv_p[1] = _mm512_permute_ps(bv[1], 0xB1);

                // b := alpha * b
                bv_p[0] = _mm512_mul_ps(alpha_img, bv_p[0]);
                bv_p[1] = _mm512_mul_ps(alpha_img, bv_p[1]);
                bv[0] = _mm512_fmaddsub_ps(alpha_real, bv[0], bv_p[0]);
                bv[1] = _mm512_fmaddsub_ps(alpha_real, bv[1], bv_p[1]);

                cgetrf_avx512_store4sc(acur, i_1, bv[0]);
                cgetrf_avx512_store4sc(acur, i_2, bv[1]);

                bv_p[0] = _mm512_permute_ps(bv[0], 0xB1);
                bv_p[1] = _mm512_permute_ps(bv[1], 0xB1);

                for(j = 1; j < ni - 1; j = j + 2)
                {
                    i_3 = j * *lda;
                    i_2 = i_1 + i_3;

                    i_4 = (j + 1) * *lda;
                    i_5 = i_1 + i_4;

                    i_6 = i_2 + 4;
                    i_7 = i_5 + 4;

                    // Load x from memory
                    x_real[0] = _mm512_set1_ps(acur[i_3].real);
                    x_img[0] = _mm512_set1_ps(acur[i_3].imag);
                    x_real[1] = _mm512_set1_ps(acur[i_4].real);
                    x_img[1] = _mm512_set1_ps(acur[i_4].imag);

                    // Y := Y - b * x
                    xv0[0] = _mm512_mul_ps(x_img[0], bv_p[0]);
                    xv0[1] = _mm512_mul_ps(x_img[0], bv_p[1]);
                    xv0[0] = _mm512_fmaddsub_ps(x_real[0], bv[0], xv0[0]);
                    xv0[1] = _mm512_fmaddsub_ps(x_real[0], bv[1], xv0[1]);

                    xv1[0] = _mm512_mul_ps(x_img[1], bv_p[0]);
                    xv1[1] = _mm512_mul_ps(x_img[1], bv_p[1]);
                    xv1[0] = _mm512_fmaddsub_ps(x_real[1], bv[0], xv1[0]);
                    xv1[1] = _mm512_fmaddsub_ps(x_real[1], bv[1], xv1[1]);

                    yv0[0] = cgetrf_avx512_load4sc(acur, i_2);
                    yv0[1] = cgetrf_avx512_load4sc(acur, i_6);
                    yv1[0] = cgetrf_avx512_load4sc(acur, i_5);
                    yv1[1] = cgetrf_avx512_load4sc(acur, i_7);

                    yv0[0] = _mm512_sub_ps(yv0[0], xv0[0]);
                    yv0[1] = _mm512_sub_ps(yv0[1], xv0[1]);
                    yv1[0] = _mm512_sub_ps(yv1[0], xv1[0]);
                    yv1[1] = _mm512_sub_ps(yv1[1], xv1[1]);

                    cgetrf_avx512_store4sc(acur, i_2, yv0[0]);
                    cgetrf_avx512_store4sc(acur, i_6, yv0[1]);
                    cgetrf_avx512_store4sc(acur, i_5, yv1[0]);
                    cgetrf_avx512_store4sc(acur, i_7, yv1[1]);
                }
                if(ni - j > 0)
                {
                    i_3 = j * *lda;
                    i_2 = i_1 + i_3;

                    i_4 = i_2 + 4;

                    // Load x from memory
                    x_real[0] = _mm512_set1_ps(acur[i_3].real);
                    x_img[0] = _mm512_set1_ps(acur[i_3].imag);

                    // Y := Y - b * x
                    xv0[0] = _mm512_mul_ps(x_img[0], bv_p[0]);
                    xv0[1] = _mm512_mul_ps(x_img[0], bv_p[1]);
                    xv0[0] = _mm512_fmaddsub_ps(x_real[0], bv[0], xv0[0]);
                    xv0[1] = _mm512_fmaddsub_ps(x_real[0], bv[1], xv0[1]);

                    yv0[0] = cgetrf_avx512_load4sc(acur, i_2);
                    yv0[1] = cgetrf_avx512_load4sc(acur, i_4);

                    yv0[0] = _mm512_sub_ps(yv0[0], xv0[0]);
                    yv0[1] = _mm512_sub_ps(yv0[1], xv0[1]);

                    cgetrf_avx512_store4sc(acur, i_2, yv0[0]);
                    cgetrf_avx512_store4sc(acur, i_4, yv0[1]);
                }
            }

            // Updates 4 rows of trailing matrix per iteration
            for(; i_1 < mi - 3; i_1 += 4)
            {
                /*-----------Trailing matrix update for LU factorisation-----------

                A00 | a01    A02
                ----|-----------
                a10 | alpha  x
                A20 | b      Y

                b := alpha * b
                Y := Y - b * x


                SIMD algorithm:

                alpha = aR1  aI1  aR1  aI1
                bv    = bR1  bI1  bR2  bI2
                bv_p  = bI1  bR1  bI2  bR2
                xv    = xR1  xI1  xR1  xI1
                yv    = yR1  yI1  yR2  yI2

                step 1 => b := alpha * b
                    bv = alpha * bv
                    bv_p = alpha * (-bv_p)
                    bv = bv - bv_p
                    bv_p = shuffle(bv)

                step 2 => Y := Y - b * x
                    temp = xv * bv
                    temp1 = xv * (-bv_p)
                    xv = temp0 - temp1
                    yv = yv - xv
                ------------------------------------------------------------------*/

                // Load alpha from memory
                bv[0] = cgetrf_avx512_load4sc(acur, i_1);
                bv_p[0] = _mm512_permute_ps(bv[0], 0xB1);

                // b := alpha * b
                bv_p[0] = _mm512_mul_ps(alpha_img, bv_p[0]);
                bv[0] = _mm512_fmaddsub_ps(alpha_real, bv[0], bv_p[0]);

                cgetrf_avx512_store4sc(acur, i_1, bv[0]);

                bv_p[0] = _mm512_permute_ps(bv[0], 0xB1);

                for(j = 1; j < ni - 1; j = j + 2)
                {
                    i_3 = j * *lda;
                    i_2 = i_1 + i_3;
                    i_4 = (j + 1) * *lda;
                    i_5 = i_1 + i_4;

                    // Load x from memory
                    x_real[0] = _mm512_set1_ps(acur[i_3].real);
                    x_img[0] = _mm512_set1_ps(acur[i_3].imag);
                    x_real[1] = _mm512_set1_ps(acur[i_4].real);
                    x_img[1] = _mm512_set1_ps(acur[i_4].imag);

                    // Y := Y - b * x
                    xv0[0] = _mm512_mul_ps(x_img[0], bv_p[0]);
                    xv0[0] = _mm512_fmaddsub_ps(x_real[0], bv[0], xv0[0]);

                    xv1[0] = _mm512_mul_ps(x_img[1], bv_p[0]);
                    xv1[0] = _mm512_fmaddsub_ps(x_real[1], bv[0], xv1[0]);

                    yv0[0] = cgetrf_avx512_load4sc(acur, i_2);
                    yv1[0] = cgetrf_avx512_load4sc(acur, i_5);

                    yv0[0] = _mm512_sub_ps(yv0[0], xv0[0]);
                    yv1[0] = _mm512_sub_ps(yv1[0], xv1[0]);

                    cgetrf_avx512_store4sc(acur, i_2, yv0[0]);
                    cgetrf_avx512_store4sc(acur, i_5, yv1[0]);
                }
                if(ni - j > 0)
                {
                    i_3 = j * *lda;
                    i_2 = i_1 + i_3;

                    // Load x from memory
                    x_real[0] = _mm512_set1_ps(acur[i_3].real);
                    x_img[0] = _mm512_set1_ps(acur[i_3].imag);

                    // Y := Y - b * x
                    xv0[0] = _mm512_mul_ps(x_img[0], bv_p[0]);
                    xv0[0] = _mm512_fmaddsub_ps(x_real[0], bv[0], xv0[0]);

                    yv0[0] = cgetrf_avx512_load4sc(acur, i_2);

                    yv0[0] = _mm512_sub_ps(yv0[0], xv0[0]);

                    cgetrf_avx512_store4sc(acur, i_2, yv0[0]);
                }
            }

            // Updates 1 row of trailing matrix per iteration
            for(; i_1 < mi; i_1++)
            {
                t_val = acur[i_1].real;
                acur[i_1].real = (t_val * z__1.real - acur[i_1].imag * z__1.imag);
                acur[i_1].imag = (t_val * z__1.imag + acur[i_1].imag * z__1.real);

                t_val = acur[i_1].real;
                z_val = acur[i_1].imag;

                for(j = 1; j < ni - 1; j = j + 2)
                {
                    i_3 = j * *lda;
                    i_2 = i_1 + i_3;

                    i_4 = (j + 1) * *lda;
                    i_5 = i_1 + i_4;

                    acur[i_2].real
                        = acur[i_2].real - t_val * acur[i_3].real + z_val * acur[i_3].imag;
                    acur[i_2].imag
                        = acur[i_2].imag - t_val * acur[i_3].imag - z_val * acur[i_3].real;

                    acur[i_5].real
                        = acur[i_5].real - t_val * acur[i_4].real + z_val * acur[i_4].imag;
                    acur[i_5].imag
                        = acur[i_5].imag - t_val * acur[i_4].imag - z_val * acur[i_4].real;
                }

                if(ni - j > 0)
                {
                    i_3 = j * *lda;
                    i_2 = i_1 + i_3;

                    acur[i_2].real
                        = acur[i_2].real - t_val * acur[i_3].real + z_val * acur[i_3].imag;
                    acur[i_2].imag
                        = acur[i_2].imag - t_val * acur[i_3].imag - z_val * acur[i_3].real;
                }
            }
        }
        else
        {
            *info = (*info == 0) ? p_idx + 1 : *info;
        }
    }
    return *info;
}

#endif
