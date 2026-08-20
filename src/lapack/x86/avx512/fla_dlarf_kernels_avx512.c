/******************************************************************************
 * * Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

/**
 * Applies the Householder reflector for DLARF (incv == 1, plain variant).
 * Required: m >= 1
 */
void fla_dlarf_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                       aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                       doublereal *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp, vd2_ntau, vd2_ltmp, vd2_htmp;
    __m128d vd2_dtmp, vd2_dtmp2, vd2_vj1;
    __m256d vd4_inp, vd4_dtmp, vd4_vj, vd4_dtmp2;
    __m256d vd4_ltmp, vd4_htmp;
    __m512d vd8_dtmp, vd8_inp, vd8_vj, vd8_dtmp2;

    /* Apply the Householder rotation                      */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */

    arows = m;
    acols = n;
    vd2_ntau = _mm_set1_pd(ntau);
    --v;
    a_buff -= ldr + 1;
    --work;

    /* Compute A**T * v */
    for(j = 1; j <= acols; j++) /* for every column c_A of A */
    {
        vd2_dtmp = _mm_setzero_pd();
        vd4_dtmp = _mm256_setzero_pd();
        vd8_dtmp = _mm512_setzero_pd();

        /* Compute tmp = c_A**T . v */
        for(k = 1; k <= (arows - 7); k += 8)
        {
            /* load column elements of A and v */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd8_dtmp2 = _mm512_mul_pd(vd8_inp, vd8_vj);
            vd8_dtmp = _mm512_add_pd(vd8_dtmp, vd8_dtmp2);
        }
        if(k <= (arows - 3))
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp2 = _mm256_mul_pd(vd4_inp, vd4_vj);
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_dtmp2);

            k += 4;
        }
        if(k < arows)
        {
            /* load column elements of A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp2 = _mm_mul_pd(vd2_inp, vd2_vj1);
            vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_dtmp2);
            k += 2;
        }
        if(k == arows)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp2 = _mm_mul_pd(vd2_inp, vd2_vj1);
            vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_dtmp2);
        }

        /* Reduce add the values in vd8_dtmp, vd4_dtmp and vd2_dtmp*/

        /* Etract Upper and lower 256 bits of vd8_dtmp*/
        vd4_ltmp = _mm512_castpd512_pd256(vd8_dtmp);
        vd4_htmp = _mm512_extractf64x4_pd(vd8_dtmp, 0x1);

        /* Add the lower and upper 256 bits with vd4_dtmp */
        vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_ltmp);
        vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_htmp);

        /* Horizontal add of dtmp */
        vd4_dtmp = _mm256_hadd_pd(vd4_dtmp, vd4_dtmp);

        /* Etract Upper and lower 128 bits of vd4_dtmp*/
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);
        /* Add the lower and upper 128 bits and store in vd2_ltmp */
        vd2_ltmp = _mm_add_pd(vd2_htmp, vd2_ltmp);

        /* Horizontal add of vd2_ltmp and vd2_dtmp */
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);

        /* Store the result in work */
        _mm_storel_pd((doublereal *)&work[j], vd2_dtmp);

        /* Compute tmp = - tau * tmp */
        vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);
        vd4_dtmp = _mm256_castpd128_pd256(vd2_dtmp);
        vd4_dtmp = _mm256_insertf128_pd(vd4_dtmp, vd2_dtmp, 0x1);

        vd8_dtmp = _mm512_castpd256_pd512(vd4_dtmp);
        vd8_dtmp = _mm512_insertf64x4(vd8_dtmp, vd4_dtmp, 0x1);

        /* Compute c_A + tmp * v */
        for(k = 1; k <= (arows - 7); k += 8)
        {
            /* load column elements of c_A and v */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd8_dtmp2 = _mm512_mul_pd(vd8_dtmp, vd8_vj);
            vd8_inp = _mm512_add_pd(vd8_dtmp2, vd8_inp);
            _mm512_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd8_inp);
        }
        if(k <= (arows - 3))
        {
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_dtmp2 = _mm256_mul_pd(vd4_dtmp, vd4_vj);
            vd4_inp = _mm256_add_pd(vd4_dtmp2, vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
            k += 4;
        }
        if(k < arows)
        {
            /* load column elements of c_A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_dtmp2 = _mm_mul_pd(vd2_dtmp, vd2_vj1);
            vd2_inp = _mm_add_pd(vd2_dtmp2, vd2_inp);
            _mm_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
            k += 2;
        }
        if(k == arows)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_dtmp2 = _mm_mul_pd(vd2_dtmp, vd2_vj1);
            vd2_inp = _mm_add_pd(vd2_dtmp2, vd2_inp);
            _mm_storel_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
        }
    }
}

/**
 * Applies the Householder reflector for DLARF1F (incv == 1, v(1) implicitly 1).
 * Required: m >= 8
 */
void fla_dlarf1f_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                         aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                         doublereal *work)
{
    aocl_int64_t k, j, ktail, nrem;
    __mmask8 tmask;
    __m512d vd8_vhead, vd8_dtmp, vd8_inp, vd8_vj;
    __m256d vd4_dtmp, vd4_inp, vd4_vj;
    __m128d vd2_dtmp, vd2_inp, vd2_vj;

    /* Apply the Householder rotation                      */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */
    /* v(1) is not stored explicitly and is assumed to be 1 */

    --v;
    a_buff -= ldr + 1;
    --work;

    /* The implicit v(1) is fed in as lane 0 of the leading vector block, so the
     * accesses to A keep the same layout as the plain DLARF kernel */
    vd8_vhead = _mm512_mask_loadu_pd(_mm512_set1_pd(1.), 0xfe, &v[1]);

    nrem = m & 7; /* rows left over by the 8-row blocks */
    ktail = m - nrem + 1;
    tmask = (__mmask8)((1 << nrem) - 1);

    for(j = 1; j <= n; j++) /* for every column c_A of A */
    {
        doublereal *acol = &a_buff[j * ldr];
        doublereal dtmp;

        /* Compute tmp = c_A**T . v */
        vd8_inp = _mm512_loadu_pd(&acol[1]);
        vd8_dtmp = _mm512_mul_pd(vd8_inp, vd8_vhead);

        for(k = 9; k < ktail; k += 8)
        {
            vd8_inp = _mm512_loadu_pd(&acol[k]);
            vd8_vj = _mm512_loadu_pd(&v[k]);
            vd8_dtmp = _mm512_add_pd(vd8_dtmp, _mm512_mul_pd(vd8_inp, vd8_vj));
        }
        if(nrem)
        {
            vd8_inp = _mm512_maskz_loadu_pd(tmask, &acol[ktail]);
            vd8_vj = _mm512_maskz_loadu_pd(tmask, &v[ktail]);
            vd8_dtmp = _mm512_add_pd(vd8_dtmp, _mm512_mul_pd(vd8_inp, vd8_vj));
        }

        /* Store the result in work */
        dtmp = _mm512_reduce_add_pd(vd8_dtmp);
        work[j] = dtmp;

        /* Compute c_A + (- tau * tmp) * v. Partially masked stores are slow, so
         * the left over rows are written with narrower full width stores */
        dtmp *= ntau;
        vd8_dtmp = _mm512_set1_pd(dtmp);
        vd4_dtmp = _mm512_castpd512_pd256(vd8_dtmp);
        vd2_dtmp = _mm256_castpd256_pd128(vd4_dtmp);

        vd8_inp = _mm512_loadu_pd(&acol[1]);
        vd8_inp = _mm512_add_pd(_mm512_mul_pd(vd8_dtmp, vd8_vhead), vd8_inp);
        _mm512_storeu_pd(&acol[1], vd8_inp);

        for(k = 9; k < ktail; k += 8)
        {
            vd8_inp = _mm512_loadu_pd(&acol[k]);
            vd8_vj = _mm512_loadu_pd(&v[k]);
            vd8_inp = _mm512_add_pd(_mm512_mul_pd(vd8_dtmp, vd8_vj), vd8_inp);
            _mm512_storeu_pd(&acol[k], vd8_inp);
        }
        k = ktail;
        if(k <= m - 3)
        {
            vd4_inp = _mm256_loadu_pd(&acol[k]);
            vd4_vj = _mm256_loadu_pd(&v[k]);
            vd4_inp = _mm256_add_pd(_mm256_mul_pd(vd4_dtmp, vd4_vj), vd4_inp);
            _mm256_storeu_pd(&acol[k], vd4_inp);
            k += 4;
        }
        if(k <= m - 1)
        {
            vd2_inp = _mm_loadu_pd(&acol[k]);
            vd2_vj = _mm_loadu_pd(&v[k]);
            vd2_inp = _mm_add_pd(_mm_mul_pd(vd2_dtmp, vd2_vj), vd2_inp);
            _mm_storeu_pd(&acol[k], vd2_inp);
            k += 2;
        }
        if(k == m)
        {
            acol[k] += dtmp * v[k];
        }
    }
}

/**
 * Applies the Householder reflector for DLARF1L (incv == 1, v(m) implicitly 1).
 * Required: m >= 1
 */
void fla_dlarf1l_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                         aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                         doublereal *work)
{
    aocl_int64_t k, j, klast, nlast;
    __mmask8 lmask;
    __m512d vd8_vlast, vd8_dtmp, vd8_inp, vd8_vj;
    __m256d vd4_dtmp, vd4_inp, vd4_vj;
    __m128d vd2_dtmp, vd2_inp, vd2_vj;

    /* Apply the Householder rotation                      */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */
    /* v(m) is not stored explicitly and is assumed to be 1 */

    --v;
    a_buff -= ldr + 1;
    --work;

    nlast = ((m - 1) & 7) + 1; /* rows in the trailing vector block */
    klast = m - nlast + 1;
    lmask = (__mmask8)((1 << nlast) - 1);

    /* The implicit v(m) is fed in as the last active lane of that block */
    vd8_vlast = _mm512_mask_loadu_pd(
        _mm512_maskz_mov_pd((__mmask8)(lmask ^ (lmask >> 1)), _mm512_set1_pd(1.)),
        (__mmask8)(lmask >> 1), &v[klast]);

    for(j = 1; j <= n; j++) /* for every column c_A of A */
    {
        doublereal *acol = &a_buff[j * ldr];
        doublereal dtmp;

        /* Compute tmp = c_A**T . v */
        vd8_dtmp = _mm512_setzero_pd();
        for(k = 1; k < klast; k += 8)
        {
            vd8_inp = _mm512_loadu_pd(&acol[k]);
            vd8_vj = _mm512_loadu_pd(&v[k]);
            vd8_dtmp = _mm512_add_pd(vd8_dtmp, _mm512_mul_pd(vd8_inp, vd8_vj));
        }
        vd8_inp = _mm512_maskz_loadu_pd(lmask, &acol[klast]);
        vd8_dtmp = _mm512_add_pd(vd8_dtmp, _mm512_mul_pd(vd8_inp, vd8_vlast));

        /* Store the result in work */
        dtmp = _mm512_reduce_add_pd(vd8_dtmp);
        work[j] = dtmp;

        /* Compute c_A + (- tau * tmp) * v */
        dtmp *= ntau;
        vd8_dtmp = _mm512_set1_pd(dtmp);

        for(k = 1; k < klast; k += 8)
        {
            vd8_inp = _mm512_loadu_pd(&acol[k]);
            vd8_vj = _mm512_loadu_pd(&v[k]);
            vd8_inp = _mm512_add_pd(_mm512_mul_pd(vd8_dtmp, vd8_vj), vd8_inp);
            _mm512_storeu_pd(&acol[k], vd8_inp);
        }
        if(nlast == 8)
        {
            vd8_inp = _mm512_loadu_pd(&acol[klast]);
            vd8_inp = _mm512_add_pd(_mm512_mul_pd(vd8_dtmp, vd8_vlast), vd8_inp);
            _mm512_storeu_pd(&acol[klast], vd8_inp);
        }
        else
        {
            /* Partially masked stores are slow, so rows klast..m-1 are written
             * with narrower full width stores and v(m) = 1 is applied directly */
            vd4_dtmp = _mm512_castpd512_pd256(vd8_dtmp);
            vd2_dtmp = _mm256_castpd256_pd128(vd4_dtmp);
            k = klast;
            if(k <= m - 4)
            {
                vd4_inp = _mm256_loadu_pd(&acol[k]);
                vd4_vj = _mm256_loadu_pd(&v[k]);
                vd4_inp = _mm256_add_pd(_mm256_mul_pd(vd4_dtmp, vd4_vj), vd4_inp);
                _mm256_storeu_pd(&acol[k], vd4_inp);
                k += 4;
            }
            if(k <= m - 2)
            {
                vd2_inp = _mm_loadu_pd(&acol[k]);
                vd2_vj = _mm_loadu_pd(&v[k]);
                vd2_inp = _mm_add_pd(_mm_mul_pd(vd2_dtmp, vd2_vj), vd2_inp);
                _mm_storeu_pd(&acol[k], vd2_inp);
                k += 2;
            }
            if(k == m - 1)
            {
                acol[k] += dtmp * v[k];
            }
            acol[m] += dtmp;
        }
    }
}

/**
 * Applies the Householder reflector from the right (incv == 1), for DLARF1F and
 * DLARF1L alike.
 *    A = A - tau * A * v * v**T
 * DGEMV and DGER operations are combined.
 * The column of A that pairs with the implicit v element 1 is not part of
 * a_buff, it is passed separately in c1. The remaining n columns start at
 * a_buff and pair with the n elements of v.
 * Columns of A are contiguous, so both passes reduce to axpys and the rows are
 * processed in register blocks to keep w in flight.
 *
 * @note: n is the number of columns of A excluding c1 (i.e. the number of columns of A minus 1)
 * Required: m >= 1
 */
void fla_dlarf1_right_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                         aocl_int64_t ldr, doublereal *v, doublereal *c1,
                                         doublereal ntau, doublereal *work)
{
    aocl_int64_t i, j;
    __m512d vd8_w, vd8_inp, vd8_vj, vd8_ntau;
    __m256d vd4_w, vd4_inp, vd4_vj, vd4_ntau;
    __m128d vd2_w, vd2_inp, vd2_vj, vd2_ntau;

    vd8_ntau = _mm512_set1_pd(ntau);
    vd4_ntau = _mm512_castpd512_pd256(vd8_ntau);
    vd2_ntau = _mm256_castpd256_pd128(vd4_ntau);

    for(i = 0; i <= m - 8; i += 8)
    {
        /* w = A(i:i+7,:) * v, with the implicit v element 1 taken from c1 */
        vd8_w = _mm512_loadu_pd(&c1[i]);
        for(j = 0; j < n; j++)
        {
            vd8_inp = _mm512_loadu_pd(&a_buff[i + j * ldr]);
            vd8_vj = _mm512_set1_pd(v[j]);
            vd8_w = _mm512_fmadd_pd(vd8_inp, vd8_vj, vd8_w);
        }
        _mm512_storeu_pd(&work[i], vd8_w);

        /* A(i:i+7,:) -= tau * w * v**T */
        vd8_inp = _mm512_loadu_pd(&c1[i]);
        _mm512_storeu_pd(&c1[i], _mm512_fmadd_pd(vd8_ntau, vd8_w, vd8_inp));
        for(j = 0; j < n; j++)
        {
            vd8_inp = _mm512_loadu_pd(&a_buff[i + j * ldr]);
            vd8_vj = _mm512_set1_pd(ntau * v[j]);
            vd8_inp = _mm512_fmadd_pd(vd8_w, vd8_vj, vd8_inp);
            _mm512_storeu_pd(&a_buff[i + j * ldr], vd8_inp);
        }
    }
    /* Partially masked stores are slow, so the left over rows are processed
     * with narrower full width accesses */
    if(i <= m - 4)
    {
        vd4_w = _mm256_loadu_pd(&c1[i]);
        for(j = 0; j < n; j++)
        {
            vd4_inp = _mm256_loadu_pd(&a_buff[i + j * ldr]);
            vd4_vj = _mm256_set1_pd(v[j]);
            vd4_w = _mm256_fmadd_pd(vd4_inp, vd4_vj, vd4_w);
        }
        _mm256_storeu_pd(&work[i], vd4_w);

        vd4_inp = _mm256_loadu_pd(&c1[i]);
        _mm256_storeu_pd(&c1[i], _mm256_fmadd_pd(vd4_ntau, vd4_w, vd4_inp));
        for(j = 0; j < n; j++)
        {
            vd4_inp = _mm256_loadu_pd(&a_buff[i + j * ldr]);
            vd4_vj = _mm256_set1_pd(ntau * v[j]);
            vd4_inp = _mm256_fmadd_pd(vd4_w, vd4_vj, vd4_inp);
            _mm256_storeu_pd(&a_buff[i + j * ldr], vd4_inp);
        }
        i += 4;
    }
    if(i <= m - 2)
    {
        vd2_w = _mm_loadu_pd(&c1[i]);
        for(j = 0; j < n; j++)
        {
            vd2_inp = _mm_loadu_pd(&a_buff[i + j * ldr]);
            vd2_vj = _mm_set1_pd(v[j]);
            vd2_w = _mm_fmadd_pd(vd2_inp, vd2_vj, vd2_w);
        }
        _mm_storeu_pd(&work[i], vd2_w);

        vd2_inp = _mm_loadu_pd(&c1[i]);
        _mm_storeu_pd(&c1[i], _mm_fmadd_pd(vd2_ntau, vd2_w, vd2_inp));
        for(j = 0; j < n; j++)
        {
            vd2_inp = _mm_loadu_pd(&a_buff[i + j * ldr]);
            vd2_vj = _mm_set1_pd(ntau * v[j]);
            vd2_inp = _mm_fmadd_pd(vd2_w, vd2_vj, vd2_inp);
            _mm_storeu_pd(&a_buff[i + j * ldr], vd2_inp);
        }
        i += 2;
    }
    if(i < m)
    {
        vd2_w = _mm_load_sd(&c1[i]);
        for(j = 0; j < n; j++)
        {
            vd2_inp = _mm_load_sd(&a_buff[i + j * ldr]);
            vd2_vj = _mm_load_sd(&v[j]);
            vd2_w = _mm_fmadd_sd(vd2_inp, vd2_vj, vd2_w);
        }
        _mm_store_sd(&work[i], vd2_w);

        vd2_inp = _mm_load_sd(&c1[i]);
        _mm_store_sd(&c1[i], _mm_fmadd_sd(vd2_ntau, vd2_w, vd2_inp));
        for(j = 0; j < n; j++)
        {
            vd2_inp = _mm_load_sd(&a_buff[i + j * ldr]);
            vd2_vj = _mm_set_sd(ntau * v[j]);
            vd2_inp = _mm_fmadd_sd(vd2_w, vd2_vj, vd2_inp);
            _mm_store_sd(&a_buff[i + j * ldr], vd2_inp);
        }
    }
}

/**
 * Folds the column of A that pairs with the implicit v element 1 into the GEMV
 * result and updates it, for DLARF1F and DLARF1L applied from the right.
 *    work = work + c1
 *    c1   = c1 + ntau * work
 * Both operands are contiguous, so the two passes fuse into one streaming loop.
 * Required: m >= 0
 */
void fla_dlarf1_right_update_c1_avx512(aocl_int64_t m, doublereal *restrict c1, doublereal ntau,
                                       doublereal *restrict work)
{
    aocl_int64_t i;
    __m512d vd8_ntau, vd8_c, vd8_w;
    __m256d vd4_ntau, vd4_c, vd4_w;
    __m128d vd2_ntau, vd2_c, vd2_w;

    vd8_ntau = _mm512_set1_pd(ntau);
    vd4_ntau = _mm512_castpd512_pd256(vd8_ntau);
    vd2_ntau = _mm256_castpd256_pd128(vd4_ntau);

    for(i = 0; i <= m - 8; i += 8)
    {
        vd8_c = _mm512_loadu_pd(&c1[i]);
        vd8_w = _mm512_add_pd(_mm512_loadu_pd(&work[i]), vd8_c);
        _mm512_storeu_pd(&work[i], vd8_w);
        _mm512_storeu_pd(&c1[i], _mm512_fmadd_pd(vd8_ntau, vd8_w, vd8_c));
    }
    /* Partially masked stores are slow, so the left over rows are written with
     * narrower full width stores */
    if(i <= m - 4)
    {
        vd4_c = _mm256_loadu_pd(&c1[i]);
        vd4_w = _mm256_add_pd(_mm256_loadu_pd(&work[i]), vd4_c);
        _mm256_storeu_pd(&work[i], vd4_w);
        _mm256_storeu_pd(&c1[i], _mm256_fmadd_pd(vd4_ntau, vd4_w, vd4_c));
        i += 4;
    }
    if(i <= m - 2)
    {
        vd2_c = _mm_loadu_pd(&c1[i]);
        vd2_w = _mm_add_pd(_mm_loadu_pd(&work[i]), vd2_c);
        _mm_storeu_pd(&work[i], vd2_w);
        _mm_storeu_pd(&c1[i], _mm_fmadd_pd(vd2_ntau, vd2_w, vd2_c));
        i += 2;
    }
    if(i < m)
    {
        vd2_c = _mm_load_sd(&c1[i]);
        vd2_w = _mm_add_sd(_mm_load_sd(&work[i]), vd2_c);
        _mm_store_sd(&work[i], vd2_w);
        _mm_store_sd(&c1[i], _mm_fmadd_sd(vd2_ntau, vd2_w, vd2_c));
    }
}
#endif