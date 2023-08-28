/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_dgeqrf_small_avx2.h
 *  @brief QR Kernels for small sizes.
 *  */

#ifdef FLA_ENABLE_AMD_OPT

/* Declaration of local variables for QR Small */
#define FLA_GEQRF_INIT_DSMALL()                                           \
    integer i, j, k;                                                      \
    integer min_m_n, kcnt, slen;                                          \
    integer acols, arows;                                                 \
    doublereal xnorm, vnorm, dtmp;                                        \
    doublereal fnorm, scale;                                              \
    doublereal med_sum, sml_sum, big_sum;                                 \
    doublereal alpha, beta, ntau;                                         \
    doublereal *v, *A, *ac;                                               \
                                                                          \
    static TLS_CLASS_SPEC int r_once = 1;                                 \
    static doublereal safmin, rsafmin;                                    \
    /* Constants chosen to minimize roundoff, */                          \
    /* according to Blue's algorithm          */                          \
    static doublereal thres_sml = 1.491668e-154;                          \
    static doublereal thres_big = 1.997919e+146;                          \
    static doublereal scale_sml = 4.498914e+161;                          \
    static doublereal scale_big = 1.111379e-162;                          \
                                                                          \
    __m128d vd2_inp, vd2_abs_inp;                                         \
    __m128d vd2_sth, vd2_bth, vd2_sscl, vd2_bscl;                         \
    __m128d vd2_sinp, vd2_binp, vd2_minp;                                 \
    __m128d vd2_smsk, vd2_bmsk, vd2_mmsk;                                 \
    __m128d vd2_ssum, vd2_bsum, vd2_msum;                                 \
    __m128d vd2_norm, vd2_vj1, vd2_dtmp;                                  \
    __m128d vd2_ntau, vd2_ltmp, vd2_htmp;                                 \
    __m128d vd2_zero = _mm_set1_pd(-0.0f);                                \
                                                                          \
    __m256d vd4_inp, vd4_abs_inp;                                         \
    __m256d vd4_sth, vd4_bth, vd4_sscl, vd4_bscl;                         \
    __m256d vd4_sinp, vd4_binp, vd4_minp;                                 \
    __m256d vd4_smsk, vd4_bmsk, vd4_mmsk;                                 \
    __m256d vd4_ssum, vd4_bsum, vd4_msum;                                 \
    __m256d vd4_norm, vd4_vj, vd4_dtmp;                                   \
    __m256d vd4_zero = _mm256_set1_pd(-0.0f);                             \
                                                                          \
    if (r_once)                                                           \
    {                                                                     \
        safmin = dlamch_("S") / dlamch_("E");                             \
        rsafmin = 1. / safmin;                                            \
        r_once = 0;                                                       \
    }                                                                     \
                                                                          \
    vd2_sth = _mm_set1_pd(thres_sml);                                     \
    vd2_bth = _mm_set1_pd(thres_big);                                     \
    vd2_sscl = _mm_set1_pd(scale_sml);                                    \
    vd2_bscl = _mm_set1_pd(scale_big);

/* Combining three accumulators of norm to get final norm */
#define FLA_GEQRF_SMALL_GET_NORM()                                        \
    fnorm = 0.0;                                                          \
    if (big_sum > 0.0)                                                    \
    {                                                                     \
        fnorm = big_sum;                                                  \
        if (med_sum > 0.0)                                                \
        {                                                                 \
            fnorm = fnorm + (med_sum * scale_big) * scale_big;            \
        }                                                                 \
        scale = scale_big;                                                \
    }                                                                     \
    else /* small sum must be non-zero */                                 \
    {                                                                     \
        doublereal ymin, ymax;                                            \
        if (med_sum > 0.0)                                                \
        {                                                                 \
            med_sum = sqrt(med_sum);                                      \
            sml_sum = sqrt(sml_sum) / scale_sml;                          \
                                                                          \
            ymin = fla_min(med_sum, sml_sum);                             \
            ymax = fla_max(med_sum, sml_sum);                             \
                                                                          \
            scale = 1.0;                                                  \
            fnorm = ymax * ymax * (1.0 + (ymin / ymax) * (ymin / ymax));  \
        }                                                                 \
        else                                                              \
        {                                                                 \
            scale = scale_sml;                                            \
            fnorm = sml_sum;                                              \
        }                                                                 \
    }                                                                     \
    xnorm = sqrt(fnorm) / scale;

/* NORM computation using 256-bit AVX2 intrinsics */
#define FLA_GEQRF_SMALL_CALC_NORM4(idx)                                   \
    /* load input and get its absolute values */                          \
    vd4_inp = _mm256_loadu_pd(&iptr[idx]);                                \
    vd4_abs_inp = _mm256_andnot_pd(vd4_zero, vd4_inp);                    \
                                                                          \
    /* segregate input values into small, medium and big */               \
    vd4_smsk = _mm256_cmp_pd(vd4_abs_inp, vd4_sth, _CMP_LT_OQ);           \
    vd4_bmsk = _mm256_cmp_pd(vd4_abs_inp, vd4_bth, _CMP_GT_OQ);           \
    vd4_mmsk = _mm256_or_pd(vd4_smsk, vd4_bmsk);                          \
                                                                          \
    /* if all inputs are in  medium range */                              \
    if (_mm256_testz_pd(vd4_mmsk, vd4_mmsk))                              \
    {                                                                     \
        vd4_msum = _mm256_fmadd_pd(vd4_inp, vd4_inp, vd4_msum);           \
    }                                                                     \
    else /* for small and large inputs */                                 \
    {                                                                     \
        has_outliers = 1;                                                 \
        vd4_sinp = _mm256_blendv_pd(vd4_zero, vd4_abs_inp, vd4_smsk);     \
        vd4_binp = _mm256_blendv_pd(vd4_zero, vd4_abs_inp, vd4_bmsk);     \
        vd4_minp = _mm256_blendv_pd(vd4_abs_inp, vd4_zero, vd4_mmsk);     \
                                                                          \
        /* scale, square and add as applicable */                         \
        vd4_sinp = _mm256_mul_pd(vd4_sinp, vd4_sscl);                     \
        vd4_binp = _mm256_mul_pd(vd4_binp, vd4_bscl);                     \
                                                                          \
        vd4_msum = _mm256_fmadd_pd(vd4_minp, vd4_minp, vd4_msum);         \
        vd4_ssum = _mm256_fmadd_pd(vd4_sinp, vd4_sinp, vd4_ssum);         \
        vd4_bsum = _mm256_fmadd_pd(vd4_binp, vd4_binp, vd4_bsum);         \
    }

/* NORM computation using 128-bit AVX intrinsics */
#define FLA_GEQRF_SMALL_CALC_NORM2()                                      \
    /* get absolute value of the vector input */                          \
    vd2_abs_inp = _mm_andnot_pd(vd2_zero, vd2_inp);                       \
                                                                          \
    /* compute flags to detect out-of-range values */                     \
    vd2_smsk = _mm_cmp_pd(vd2_abs_inp, vd2_sth, _CMP_LT_OQ);              \
    vd2_bmsk = _mm_cmp_pd(vd2_abs_inp, vd2_bth, _CMP_GT_OQ);              \
    vd2_mmsk = _mm_or_pd(vd2_smsk, vd2_bmsk);                             \
                                                                          \
    /* if all inputs are in  medium range */                              \
    if (_mm_testz_pd(vd2_mmsk, vd2_mmsk))                                 \
    {                                                                     \
        vd2_msum = _mm_fmadd_pd(vd2_inp, vd2_inp, vd2_msum);              \
    }                                                                     \
    else /* for small and large inputs */                                 \
    {                                                                     \
        has_outliers = 1;                                                 \
        vd2_sinp = _mm_blendv_pd(vd2_zero, vd2_abs_inp, vd2_smsk);        \
        vd2_binp = _mm_blendv_pd(vd2_zero, vd2_abs_inp, vd2_bmsk);        \
        vd2_minp = _mm_blendv_pd(vd2_abs_inp, vd2_zero, vd2_mmsk);        \
                                                                          \
        /* scale, square and add as applicable */                         \
        vd2_sinp = _mm_mul_pd(vd2_sinp, vd2_sscl);                        \
        vd2_binp = _mm_mul_pd(vd2_binp, vd2_bscl);                        \
        vd2_ssum = _mm_fmadd_pd(vd2_sinp, vd2_sinp, vd2_ssum);            \
        vd2_bsum = _mm_fmadd_pd(vd2_binp, vd2_binp, vd2_bsum);            \
        vd2_msum = _mm_fmadd_pd(vd2_minp, vd2_minp, vd2_msum);            \
    }

#define FLA_ELEM_REFLECTOR_GENERATE_DSMALL(i, m, n, tau)                     \
    /* calculate norm of sub-diagonal elements in current column */          \
    med_sum = sml_sum = big_sum = 0.;                                        \
    vd2_msum = _mm_setzero_pd();                                             \
    vd2_ssum = _mm_setzero_pd();                                             \
    vd2_bsum = _mm_setzero_pd();                                             \
                                                                             \
    /* process two inputs per iteration */                                   \
    for (j = 1; j <= (slen - 1); j += 2)                                     \
    {                                                                        \
        vd2_inp = _mm_loadu_pd(&iptr[j]);                                    \
        FLA_GEQRF_SMALL_CALC_NORM2();                                        \
    }                                                                        \
                                                                             \
    if (j == slen)                                                           \
    {                                                                        \
        /* load input and get its absolute values */                         \
        vd2_inp = _mm_load_sd(&iptr[j]);                                     \
        FLA_GEQRF_SMALL_CALC_NORM2();                                        \
    }                                                                        \
                                                                             \
    /* Get all the three sums */                                             \
    med_sum = vd2_msum[0] + vd2_msum[1];                                     \
    /* Combining outlier accumulators if non-zero */                         \
    if (has_outliers)                                                        \
    {                                                                        \
        sml_sum = vd2_ssum[0] + vd2_ssum[1];                                 \
        big_sum = vd2_bsum[0] + vd2_bsum[1];                                 \
        FLA_GEQRF_SMALL_GET_NORM();                                          \
    }                                                                        \
    else                                                                     \
    {                                                                        \
        xnorm = sqrt(med_sum);                                               \
    }                                                                        \
                                                                             \
    /* Compute Householder Reflector parameters */                           \
    if (xnorm == 0.) /* Sub-diagonal elements are already zero */            \
    {                                                                        \
        tau[i] = 0.;                                                         \
        beta = 0.;                                                           \
    }                                                                        \
    else /* Non-zero sub-diagonal elements */                                \
    {                                                                        \
        /* Part 1: Compute Householder vector 'v' and tau */                 \
                                                                             \
        v = &a[i + i * *lda - 1];                                            \
        A = &a[i + (i + 1) * *lda];                                          \
        alpha = v[1];                                                        \
        /* check for not a number */                                         \
        if (alpha != alpha || xnorm != xnorm)                                \
        {                                                                    \
            beta = alpha + xnorm;                                            \
        }                                                                    \
        else                                                                 \
        {                                                                    \
            doublereal px0, w, z;                                            \
                                                                             \
            px0 = f2c_abs(alpha);                                            \
            w = fla_max(px0, xnorm);                                         \
            z = fla_min(px0, xnorm);                                         \
                                                                             \
            z = z / w;                                                       \
            beta = w * sqrt(z * z + 1);                                      \
        }                                                                    \
        beta = (alpha >= 0.) ? -beta : beta;                                 \
                                                                             \
        /* Scale-up the inputs for small norm */                             \
        for (kcnt = 0; (f2c_abs(beta) < safmin && kcnt <= 20); kcnt++)       \
        {                                                                    \
            dscal_(&slen, &rsafmin, &v[2], &c__1);                           \
            beta = beta * rsafmin;                                           \
            alpha = alpha * rsafmin;                                         \
        }                                                                    \
                                                                             \
        /* Calculate tau and v */                                            \
        tau[i] = (beta - alpha) / beta;                                      \
        vnorm = 1. / (alpha - beta);                                         \
        /* Scale current column by norm to get v */                          \
        vd2_norm = _mm_set1_pd(vnorm);                                       \
                                                                             \
        /* Normalize using SIMD */                                           \
        for (j = 1; j <= (slen - 1); j += 2)                                 \
        {                                                                    \
            vd2_vj1 = _mm_loadu_pd((doublereal const *) &v[j + 1]);          \
            vd2_vj1 = _mm_mul_pd(vd2_vj1, vd2_norm);                         \
            _mm_storeu_pd((doublereal *) &v[j + 1], vd2_vj1);                \
        }                                                                    \
        if (j == slen)                                                       \
        {                                                                    \
            vd2_vj1 = _mm_loaddup_pd((doublereal const *) &v[j + 1]);        \
            vd2_vj1 = _mm_mul_pd(vd2_vj1, vd2_norm);                         \
            _mm_storel_pd((doublereal *) &v[j + 1], vd2_vj1);                \
        }                                                                    \
        /* Scale-down beta */                                                \
        for ( ; kcnt >= 1; kcnt--)                                           \
        {                                                                    \
            beta = beta * safmin;                                            \
        }                                                                    \
    }

#define FLA_ELEM_REFLECTOR_APPLY_DSMALL(i, m, n, tau)                        \
    if (xnorm != 0.) /* Sub-diagonal elements are already zero */            \
    {                                                                        \
        /* Part 2: Apply the Householder rotation              */            \
        /* on the rest of the matrix                           */            \
        /*    A = A - tau * v * v**T * A                       */            \
        /*      = A - v * tau * (A**T * v)**T                  */            \
                                                                             \
        arows = *m - i + 1;                                                  \
        acols = *n - i;                                                      \
        v[1] = 1.;                                                           \
        ntau = -tau[i];                                                      \
        vd2_ntau = _mm_set1_pd(ntau);                                        \
                                                                             \
        /* work = A**T * v */                                                \
        for (j = 1; j <= acols; j++) /* for every column c_A of A */         \
        {                                                                    \
            ac = &A[(j - 1) * *lda - 1];                                     \
                                                                             \
            /* Compute tmp = c_A**T . v */                                   \
            dtmp = 0;                                                        \
            for (k = 1; k <= (arows - 1); k += 2)                            \
            {                                                                \
                dtmp = dtmp + ac[k] * v[k];                                  \
                dtmp = dtmp + ac[k + 1] * v[k + 1];                          \
            }                                                                \
            if (k == arows)                                                  \
            {                                                                \
                dtmp = dtmp + ac[k] * v[k];                                  \
            }                                                                \
            vd2_dtmp = _mm_set1_pd(dtmp);                                    \
                                                                             \
            /* Compute tmp = -tau * tmp */                                   \
            vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);                       \
                                                                             \
            /* Compute c_A + tmp * v */                                      \
            for (k = 1; k <= (arows - 1); k += 2)                            \
            {                                                                \
                /* load column elements of c_A and v */                      \
                vd2_inp = _mm_loadu_pd((const doublereal *) &ac[k]);         \
                vd2_vj1 = _mm_loadu_pd((const doublereal *) &v[k]);          \
                                                                             \
                /* mul by dtmp, add and store */                             \
                vd2_inp = _mm_fmadd_pd(vd2_vj1, vd2_dtmp, vd2_inp);          \
                _mm_storeu_pd((doublereal *) &ac[k], vd2_inp);               \
            }                                                                \
            if (k == arows)                                                  \
            {                                                                \
                /* load single remaining element from c_A and v */           \
                vd2_inp = _mm_load_sd((const doublereal *) &ac[k]);          \
                vd2_vj1 = _mm_load_sd((const doublereal *) &v[k]);           \
                                                                             \
                /* multiply with tau and store */                            \
                vd2_inp = _mm_fmadd_pd(vd2_vj1, vd2_dtmp, vd2_inp);          \
                _mm_store_sd((doublereal *) &ac[k], vd2_inp);                \
            }                                                                \
        }                                                                    \
        v[1] = beta;                                                         \
    }

#define FLA_ELEM_REFLECTOR_GENERATE_DLARGE(i, m, n, tau)                     \
    /* calculate norm of sub-diagonal elements in current column */          \
    med_sum = sml_sum = big_sum = 0.;                                        \
    vd4_sth = _mm256_set1_pd(thres_sml);                                     \
    vd4_bth = _mm256_set1_pd(thres_big);                                     \
    vd4_sscl = _mm256_set1_pd(scale_sml);                                    \
    vd4_bscl = _mm256_set1_pd(scale_big);                                    \
                                                                             \
    vd4_msum = _mm256_setzero_pd();                                          \
    vd4_ssum = _mm256_setzero_pd();                                          \
    vd4_bsum = _mm256_setzero_pd();                                          \
                                                                             \
    /* process four inputs per iteration */                                  \
    for (j = 1; j <= (slen - 3); j += 4)                                     \
    {                                                                        \
        FLA_GEQRF_SMALL_CALC_NORM4(j);                                       \
    }                                                                        \
                                                                             \
    if (j <= slen)                                                           \
    {   /* process remaining iterations */                                   \
        vd2_msum = _mm_setzero_pd();                                         \
        vd2_ssum = _mm_setzero_pd();                                         \
        vd2_bsum = _mm_setzero_pd();                                         \
                                                                             \
        /* process two inputs per iteration */                               \
        for ( ; j <= slen; j++)                                              \
        {                                                                    \
            vd2_inp = _mm_loaddup_pd(&iptr[j]);                              \
            FLA_GEQRF_SMALL_CALC_NORM2();                                    \
        }                                                                    \
        /* Get all the three sums */                                         \
        med_sum = vd4_msum[0] + vd4_msum[1] + vd4_msum[2]                    \
                + vd4_msum[3] + vd2_msum[0];                                 \
        sml_sum = vd4_ssum[0] + vd4_ssum[1] + vd4_ssum[2]                    \
                + vd4_ssum[3] + vd2_ssum[0];                                 \
        big_sum = vd4_bsum[0] + vd4_bsum[1] + vd4_bsum[2]                    \
                + vd4_bsum[3] + vd2_bsum[0];                                 \
    }                                                                        \
    else                                                                     \
    {                                                                        \
        /* Get all the three sums in case of no remaining iterations */      \
        med_sum = vd4_msum[0] + vd4_msum[1] + vd4_msum[2] + vd4_msum[3];     \
        sml_sum = vd4_ssum[0] + vd4_ssum[1] + vd4_ssum[2] + vd4_ssum[3];     \
        big_sum = vd4_bsum[0] + vd4_bsum[1] + vd4_bsum[2] + vd4_bsum[3];     \
    }                                                                        \
                                                                             \
    /* Combining outlier accumulators if non-zero */                         \
    if (has_outliers)                                                        \
    {                                                                        \
        FLA_GEQRF_SMALL_GET_NORM();                                          \
    }                                                                        \
    else                                                                     \
    {                                                                        \
        xnorm = sqrt(med_sum);                                               \
    }                                                                        \
                                                                             \
    /* Compute Householder Reflector parameters */                           \
    if (xnorm == 0.) /* Sub-diagonal elements are already zero */            \
    {                                                                        \
        tau[i] = 0.;                                                         \
        beta = 0.;                                                           \
    }                                                                        \
    else /* Non-zero sub-diagonal elements */                                \
    {                                                                        \
        /* Part 1: Compute Householder vector 'v' and tau */                 \
                                                                             \
        v = &a[i + i * *lda - 1];                                            \
        A = &a[i + (i + 1) * *lda];                                          \
        alpha = v[1];                                                        \
                                                                             \
        /* Compute Householder rotated vector */                             \
        if (alpha != alpha || xnorm != xnorm) /* check for not a number */   \
        {                                                                    \
            beta = alpha + xnorm;                                            \
        }                                                                    \
        else                                                                 \
        {                                                                    \
            doublereal px0, w, z;                                            \
                                                                             \
            px0 = f2c_abs(alpha);                                            \
            w = fla_max(px0, xnorm);                                         \
            z = fla_min(px0, xnorm);                                         \
                                                                             \
            z = z / w;                                                       \
            beta = w * sqrt(z * z + 1);                                      \
        }                                                                    \
        beta = (alpha >= 0.) ? -beta : beta;                                 \
                                                                             \
        /* Scale-up the inputs for small norm */                             \
        for (kcnt = 0; (f2c_abs(beta) < safmin && kcnt <= 20); kcnt++)       \
        {                                                                    \
            dscal_(&slen, &rsafmin, &v[2], &c__1);                           \
            beta = beta * rsafmin;                                           \
            alpha = alpha * rsafmin;                                         \
        }                                                                    \
                                                                             \
        /* Calculate tau and v */                                            \
        tau[i] = (beta - alpha) / beta;                                      \
        vnorm = 1. / (alpha - beta);                                         \
        /* Scale current column by norm to get v */                          \
        vd4_norm = _mm256_set1_pd(vnorm);                                    \
                                                                             \
        /* Normalize using SIMD */                                           \
        for (j = 1; j <= (slen - 3); j += 4)                                 \
        {                                                                    \
            vd4_vj = _mm256_loadu_pd((doublereal const *) &v[j + 1]);        \
            vd4_vj = _mm256_mul_pd(vd4_vj, vd4_norm);                        \
            _mm256_storeu_pd((doublereal *) &v[j + 1], vd4_vj);              \
        }                                                                    \
        /* Remaining iterations through 128-bit SIMD */                      \
        if (j <= slen)                                                       \
        {                                                                    \
            vd2_norm = _mm_set1_pd(vnorm);                                   \
            for ( ; j <= (slen - 1); j += 2)                                 \
            {                                                                \
                vd2_vj1 = _mm_loadu_pd((doublereal const *) &v[j + 1]);      \
                vd2_vj1 = _mm_mul_pd(vd2_vj1, vd2_norm);                     \
                _mm_storeu_pd((doublereal *) &v[j + 1], vd2_vj1);            \
            }                                                                \
            if (j == slen)                                                   \
            {                                                                \
                vd2_vj1 = _mm_loaddup_pd((doublereal const *) &v[j + 1]);    \
                vd2_vj1 = _mm_mul_pd(vd2_vj1, vd2_norm);                     \
                _mm_storel_pd((doublereal *) &v[j + 1], vd2_vj1);            \
            }                                                                \
        }                                                                    \
        /* Scale-down beta */                                                \
        for ( ; kcnt >= 1; kcnt--)                                           \
        {                                                                    \
            beta = beta * safmin;                                            \
        }                                                                    \
    }

#define FLA_ELEM_REFLECTOR_APPLY_DLARGE(i, m, n, tau)                        \
    if (xnorm != 0.) /* Sub-diagonal elements are already zero */            \
    {                                                                        \
        /* Part 2: Apply the Householder rotation              */            \
        /* on the rest of the matrix                           */            \
        /*    A = A - tau * v * v**T * A                       */            \
        /*      = A - v * tau * (A**T * v)**T                  */            \
                                                                             \
        arows = *m - i + 1;                                                  \
        acols = *n - i;                                                      \
        v[1] = 1.;                                                           \
        ntau = -tau[i];                                                      \
        vd2_ntau = _mm_set1_pd(ntau);                                        \
                                                                             \
        /* work = A**T * v */                                                \
        for (j = 1; j <= acols; j++) /* for every column c_A of A */         \
        {                                                                    \
            ac = &A[(j - 1) * *lda - 1];                                     \
            vd2_dtmp = _mm_setzero_pd();                                     \
            vd4_dtmp = _mm256_setzero_pd();                                  \
                                                                             \
            /* Compute tmp = c_A**T . v */                                   \
            for (k = 1; k <= (arows - 3); k += 4)                            \
            {                                                                \
                /* load column elements of A and v */                        \
                vd4_inp = _mm256_loadu_pd((const doublereal *) &ac[k]);      \
                vd4_vj  = _mm256_loadu_pd((const doublereal *) &v[k]);       \
                                                                             \
                /* take dot product */                                       \
                vd4_dtmp = _mm256_fmadd_pd(vd4_inp, vd4_vj, vd4_dtmp);       \
            }                                                                \
            if (k < arows)                                                   \
            {                                                                \
                /* load column elements of A and v */                        \
                vd2_inp = _mm_loadu_pd((const doublereal *) &ac[k]);         \
                vd2_vj1 = _mm_loadu_pd((const doublereal *) &v[k]);          \
                                                                             \
                /* take dot product */                                       \
                vd2_dtmp = _mm_fmadd_pd(vd2_inp, vd2_vj1, vd2_dtmp);         \
                k += 2;                                                      \
            }                                                                \
            if (k == arows)                                                  \
            {                                                                \
                /* load single remaining element from c_A and v */           \
                vd2_inp = _mm_load_sd((const doublereal *) &ac[k]);          \
                vd2_vj1 = _mm_load_sd((const doublereal *) &v[k]);           \
                                                                             \
                /* take dot product */                                       \
                vd2_dtmp = _mm_fmadd_pd(vd2_inp, vd2_vj1, vd2_dtmp);         \
            }                                                                \
            /* Horizontal add of dtmp */                                     \
            vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);                     \
            vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);                 \
                                                                             \
            vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);                       \
            vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_htmp);                       \
            vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);                      \
                                                                             \
            /* Compute tmp = - tau * tmp */                                  \
            vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);                       \
            vd4_dtmp = _mm256_castpd128_pd256(vd2_dtmp);                     \
            vd4_dtmp = _mm256_insertf128_pd(vd4_dtmp, vd2_dtmp, 0x1);        \
                                                                             \
            /* alternate for above 2 instructions which do not  */           \
            /* compile for older gcc versions (7 and below).    */           \
            /* Both will be same in terms of latency though     */           \
            /* vd4_dtmp = _mm256_set_m128d(vd2_dtmp, vd2_dtmp); */           \
                                                                             \
            /* Compute c_A + tmp * v */                                      \
            for (k = 1; k <= (arows - 3); k += 4)                            \
            {                                                                \
                /* load column elements of c_A and v */                      \
                vd4_inp = _mm256_loadu_pd((const doublereal *) &ac[k]);      \
                vd4_vj  = _mm256_loadu_pd((const doublereal *) &v[k]);       \
                                                                             \
                /* mul by dtmp, add and store */                             \
                vd4_inp = _mm256_fmadd_pd(vd4_dtmp, vd4_vj, vd4_inp);        \
                _mm256_storeu_pd((doublereal *) &ac[k], vd4_inp);            \
            }                                                                \
            if (k < arows)                                                   \
            {                                                                \
                /* load column elements of c_A and v */                      \
                vd2_inp = _mm_loadu_pd((const doublereal *) &ac[k]);         \
                vd2_vj1 = _mm_loadu_pd((const doublereal *) &v[k]);          \
                                                                             \
                /* mul by dtmp, add and store */                             \
                vd2_inp = _mm_fmadd_pd(vd2_dtmp, vd2_vj1, vd2_inp);          \
                _mm_storeu_pd((doublereal *) &ac[k], vd2_inp);               \
                k += 2;                                                      \
            }                                                                \
            if (k == arows)                                                  \
            {                                                                \
                /* load single remaining element from c_A and v */           \
                vd2_inp = _mm_load_sd((const doublereal *) &ac[k]);          \
                vd2_vj1 = _mm_load_sd((const doublereal *) &v[k]);           \
                                                                             \
                /* mul by dtmp, add and store */                             \
                vd2_inp = _mm_fmadd_pd(vd2_dtmp, vd2_vj1, vd2_inp);          \
                _mm_storel_pd((doublereal *) &ac[k], vd2_inp);               \
            }                                                                \
        }                                                                    \
        v[1] = beta;                                                         \
    }
#endif
