/******************************************************************************
 * * Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT
__attribute__((aligned(512))) void fla_dlarf_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n,
                                                                   doublereal *a_buff,
                                                                   aocl_int64_t ldr, doublereal *v,
                                                                   doublereal ntau,
                                                                   doublereal *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp, vd2_vj1, vd2_dtmp2;
    __m256d vd4_dtmp, vd4_inp, vd4_vj, vd4_dtmp2;
    __m128d vd2_ltmp, vd2_htmp;

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

        /* Compute tmp = c_A**T . v */
        for(k = 1; k <= (arows - 3); k += 4)
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp2 = _mm256_mul_pd(vd4_inp, vd4_vj);
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_dtmp2);
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
        /* Horizontal add of dtmp */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);
        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_htmp);
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);

        /* Store the result in work */
        _mm_storel_pd((doublereal *)&work[j], vd2_dtmp);

        /* Compute tmp = - tau * tmp */
        vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);
        vd4_dtmp = _mm256_castpd128_pd256(vd2_dtmp);
        vd4_dtmp = _mm256_insertf128_pd(vd4_dtmp, vd2_dtmp, 0x1);

        /* alternate for above 2 instructions which do not  */
        /* compile for older gcc versions (7 and below).    */
        /* Both will be same in terms of latency though     */
        /* vd4_dtmp = _mm256_set_m128d(vd2_dtmp, vd2_dtmp); */

        /* Compute c_A + tmp * v */
        for(k = 1; k <= (arows - 3); k += 4)
        {
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_dtmp2 = _mm256_mul_pd(vd4_dtmp, vd4_vj);
            vd4_inp = _mm256_add_pd(vd4_dtmp2, vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
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
 * Required: m >= 4
 */
__attribute__((aligned(512))) void fla_dlarf1f_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n,
                                                                     doublereal *a_buff,
                                                                     aocl_int64_t ldr,
                                                                     doublereal *v, doublereal ntau,
                                                                     doublereal *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp, vd2_vj1;
    __m256d vd4_dtmp, vd4_inp, vd4_vj, vd4_vhead;
    __m128d vd2_ltmp, vd2_htmp;

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

    /* v(1) is not stored and is assumed to be 1. Substitute it in the leading */
    /* vector block of v so that no scalar prologue/epilogue is needed for the */
    /* vector sized blocks.                                                    */
    /* The API does not reference v(1), so it is masked out of the load instead */
    /* of being loaded and discarded.                                          */
    vd4_vhead = _mm256_blend_pd(
        _mm256_maskload_pd((const doublereal *)&v[1], _mm256_set_epi64x(-1, -1, -1, 0)),
        _mm256_set1_pd(1.), 0x1);

    /* Compute A**T * v */
    for(j = 1; j <= acols; j++) /* for every column c_A of A */
    {
        vd2_dtmp = _mm_setzero_pd();

        /* Compute tmp = c_A**T . v, the first block carries v(1) = 1 */
        vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[1 + j * ldr]);
        vd4_dtmp = _mm256_mul_pd(vd4_inp, vd4_vhead);

        for(k = 5; k <= (arows - 3); k += 4)
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, _mm256_mul_pd(vd4_inp, vd4_vj));
        }
        if(k < arows)
        {
            /* load column elements of A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp = _mm_add_pd(vd2_dtmp, _mm_mul_pd(vd2_inp, vd2_vj1));
            k += 2;
        }
        if(k == arows)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp = _mm_add_pd(vd2_dtmp, _mm_mul_pd(vd2_inp, vd2_vj1));
        }
        /* Horizontal add of dtmp */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);
        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_htmp);
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);

        /* Store the result in work */
        _mm_storel_pd((doublereal *)&work[j], vd2_dtmp);

        /* Compute tmp = - tau * tmp */
        vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);
        vd4_dtmp = _mm256_castpd128_pd256(vd2_dtmp);
        vd4_dtmp = _mm256_insertf128_pd(vd4_dtmp, vd2_dtmp, 0x1);

        /* alternate for above 2 instructions which do not  */
        /* compile for older gcc versions (7 and below).    */
        /* Both will be same in terms of latency though     */
        /* vd4_dtmp = _mm256_set_m128d(vd2_dtmp, vd2_dtmp); */

        /* Compute c_A + tmp * v, the first block carries v(1) = 1 */
        vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[1 + j * ldr]);
        vd4_inp = _mm256_add_pd(_mm256_mul_pd(vd4_dtmp, vd4_vhead), vd4_inp);
        _mm256_storeu_pd((doublereal *)&a_buff[1 + j * ldr], vd4_inp);

        for(k = 5; k <= (arows - 3); k += 4)
        {
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_inp = _mm256_add_pd(_mm256_mul_pd(vd4_dtmp, vd4_vj), vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
        }
        if(k < arows)
        {
            /* load column elements of c_A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_inp = _mm_add_pd(_mm_mul_pd(vd2_dtmp, vd2_vj1), vd2_inp);
            _mm_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
            k += 2;
        }
        if(k == arows)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_inp = _mm_add_pd(_mm_mul_pd(vd2_dtmp, vd2_vj1), vd2_inp);
            _mm_storel_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
        }
    }
}

/**
 * Applies the Householder reflector for DLARF1L (incv == 1, v(m) implicitly 1).
 * Required: m >= 1
 */
__attribute__((aligned(512))) void fla_dlarf1l_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n,
                                                                     doublereal *a_buff,
                                                                     aocl_int64_t ldr,
                                                                     doublereal *v, doublereal ntau,
                                                                     doublereal *work)
{
    aocl_int64_t acols, arows, nlast, klast;
    aocl_int64_t k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp, vd2_vlast;
    __m256d vd4_dtmp, vd4_inp, vd4_vj, vd4_vlast;
    __m128d vd2_ltmp, vd2_htmp;

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

    /* v(arows) is not stored and is assumed to be 1. Substitute it in the      */
    /* trailing vector block of v so that no scalar prologue/epilogue is needed */
    /* for the vector sized blocks.                                             */
    /* The API does not reference v(arows), so it is masked out of the loads     */
    /* instead of being loaded and discarded.                                    */
    nlast = ((arows - 1) & 3) + 1;
    klast = arows - nlast + 1;
    vd4_vlast = _mm256_setzero_pd();
    vd2_vlast = _mm_setzero_pd();
    if(nlast == 4)
    {
        vd4_vlast = _mm256_blend_pd(
            _mm256_maskload_pd((const doublereal *)&v[klast], _mm256_set_epi64x(0, -1, -1, -1)),
            _mm256_set1_pd(1.), 0x8);
    }
    else if(nlast == 3)
    {
        vd2_vlast = _mm_loadu_pd((const doublereal *)&v[klast]);
    }
    else if(nlast == 2)
    {
        vd2_vlast
            = _mm_blend_pd(_mm_maskload_pd((const doublereal *)&v[klast], _mm_set_epi64x(0, -1)),
                           _mm_set1_pd(1.), 0x2);
    }

    /* Compute A**T * v */
    for(j = 1; j <= acols; j++) /* for every column c_A of A */
    {
        vd2_dtmp = _mm_setzero_pd();
        vd4_dtmp = _mm256_setzero_pd();

        /* Compute tmp = c_A**T . v */
        for(k = 1; k < klast; k += 4)
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, _mm256_mul_pd(vd4_inp, vd4_vj));
        }
        /* trailing block, it carries v(arows) = 1 */
        if(nlast == 4)
        {
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[klast + j * ldr]);
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, _mm256_mul_pd(vd4_inp, vd4_vlast));
        }
        else
        {
            if(nlast > 1)
            {
                vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[klast + j * ldr]);
                vd2_dtmp = _mm_add_pd(vd2_dtmp, _mm_mul_pd(vd2_inp, vd2_vlast));
            }
            if(nlast != 2)
            {
                /* v(arows) = 1, add the contribution of v(arows) to tmp */
                vd2_dtmp = _mm_add_sd(vd2_dtmp,
                                      _mm_load_sd((const doublereal *)&a_buff[arows + j * ldr]));
            }
        }
        /* Horizontal add of dtmp */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);
        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_htmp);
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);

        /* Store the result in work */
        _mm_storel_pd((doublereal *)&work[j], vd2_dtmp);

        /* Compute tmp = - tau * tmp */
        vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);
        vd4_dtmp = _mm256_castpd128_pd256(vd2_dtmp);
        vd4_dtmp = _mm256_insertf128_pd(vd4_dtmp, vd2_dtmp, 0x1);

        /* alternate for above 2 instructions which do not  */
        /* compile for older gcc versions (7 and below).    */
        /* Both will be same in terms of latency though     */
        /* vd4_dtmp = _mm256_set_m128d(vd2_dtmp, vd2_dtmp); */

        /* Compute c_A + tmp * v */
        for(k = 1; k < klast; k += 4)
        {
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_inp = _mm256_add_pd(_mm256_mul_pd(vd4_dtmp, vd4_vj), vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
        }
        /* trailing block, it carries v(arows) = 1 */
        if(nlast == 4)
        {
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[klast + j * ldr]);
            vd4_inp = _mm256_add_pd(_mm256_mul_pd(vd4_dtmp, vd4_vlast), vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[klast + j * ldr], vd4_inp);
        }
        else
        {
            if(nlast > 1)
            {
                vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[klast + j * ldr]);
                vd2_inp = _mm_add_pd(_mm_mul_pd(vd2_dtmp, vd2_vlast), vd2_inp);
                _mm_storeu_pd((doublereal *)&a_buff[klast + j * ldr], vd2_inp);
            }
            if(nlast != 2)
            {
                /* v(arows) = 1, add the contribution of v(arows) to c_A */
                vd2_inp = _mm_load_sd(&a_buff[arows + j * ldr]);
                vd2_inp = _mm_add_sd(vd2_inp, vd2_dtmp);
                _mm_storel_pd(&a_buff[arows + j * ldr], vd2_inp);
            }
        }
    }
}

/* Apply the Householder rotation from the right                   */
/*    A = A - tau * A * v * v**T                                   */
/* DGEMV and DGER operations are combined.                         */
/* The column of A that pairs with the implicit v element 1 is not */
/* part of a_buff, it is passed separately in c1. The remaining n  */
/* columns start at a_buff and pair with the n elements of v.      */
/* Columns of A are contiguous, so both passes reduce to axpys and */
/* the rows are processed in register blocks to keep w in flight.  */
__attribute__((aligned(512))) void
    fla_dlarf1_right_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                      aocl_int64_t ldr, doublereal *v, doublereal *c1,
                                      doublereal ntau, doublereal *work)
{
    aocl_int64_t i, j;
    __m256d vd4_w, vd4_inp, vd4_vj, vd4_ntau;
    __m128d vd2_w, vd2_inp, vd2_vj, vd2_ntau;

    vd4_ntau = _mm256_set1_pd(ntau);
    vd2_ntau = _mm_set1_pd(ntau);

    for(i = 0; i <= m - 4; i += 4)
    {
        /* w = A(i:i+3,:) * v, with the implicit v element 1 taken from c1 */
        vd4_w = _mm256_loadu_pd((const doublereal *)&c1[i]);
        for(j = 0; j < n; j++)
        {
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[i + j * ldr]);
            vd4_vj = _mm256_broadcast_sd((const doublereal *)&v[j]);
            vd4_w = _mm256_fmadd_pd(vd4_inp, vd4_vj, vd4_w);
        }
        _mm256_storeu_pd((doublereal *)&work[i], vd4_w);

        /* A(i:i+3,:) -= tau * w * v**T */
        vd4_inp = _mm256_loadu_pd((const doublereal *)&c1[i]);
        _mm256_storeu_pd((doublereal *)&c1[i], _mm256_fmadd_pd(vd4_ntau, vd4_w, vd4_inp));
        for(j = 0; j < n; j++)
        {
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[i + j * ldr]);
            vd4_vj = _mm256_set1_pd(ntau * v[j]);
            vd4_inp = _mm256_fmadd_pd(vd4_w, vd4_vj, vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[i + j * ldr], vd4_inp);
        }
    }
    if(i <= m - 2)
    {
        vd2_w = _mm_loadu_pd((const doublereal *)&c1[i]);
        for(j = 0; j < n; j++)
        {
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[i + j * ldr]);
            vd2_vj = _mm_load1_pd((const doublereal *)&v[j]);
            vd2_w = _mm_fmadd_pd(vd2_inp, vd2_vj, vd2_w);
        }
        _mm_storeu_pd((doublereal *)&work[i], vd2_w);

        vd2_inp = _mm_loadu_pd((const doublereal *)&c1[i]);
        _mm_storeu_pd((doublereal *)&c1[i], _mm_fmadd_pd(vd2_ntau, vd2_w, vd2_inp));
        for(j = 0; j < n; j++)
        {
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[i + j * ldr]);
            vd2_vj = _mm_set1_pd(ntau * v[j]);
            vd2_inp = _mm_fmadd_pd(vd2_w, vd2_vj, vd2_inp);
            _mm_storeu_pd((doublereal *)&a_buff[i + j * ldr], vd2_inp);
        }
        i += 2;
    }
    if(i < m)
    {
        vd2_w = _mm_load_sd((const doublereal *)&c1[i]);
        for(j = 0; j < n; j++)
        {
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[i + j * ldr]);
            vd2_vj = _mm_load_sd((const doublereal *)&v[j]);
            vd2_w = _mm_fmadd_sd(vd2_inp, vd2_vj, vd2_w);
        }
        _mm_storel_pd((doublereal *)&work[i], vd2_w);

        vd2_inp = _mm_load_sd((const doublereal *)&c1[i]);
        _mm_storel_pd((doublereal *)&c1[i], _mm_fmadd_sd(vd2_ntau, vd2_w, vd2_inp));
        for(j = 0; j < n; j++)
        {
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[i + j * ldr]);
            vd2_vj = _mm_set_sd(ntau * v[j]);
            vd2_inp = _mm_fmadd_sd(vd2_w, vd2_vj, vd2_inp);
            _mm_storel_pd((doublereal *)&a_buff[i + j * ldr], vd2_inp);
        }
    }
}

/* Folds the column of A that pairs with the implicit v element 1 into the  */
/* GEMV result and updates it, for DLARF1F and DLARF1L from the right.      */
/*    work = work + c1                                                      */
/*    c1   = c1 + ntau * work                                               */
/* Both operands are contiguous, so the two passes fuse into one loop.      */
__attribute__((aligned(512))) void fla_dlarf1_right_update_c1_avx2(aocl_int64_t m,
                                                                   doublereal *restrict c1,
                                                                   doublereal ntau,
                                                                   doublereal *restrict work)
{
    aocl_int64_t i;
    __m256d vd4_ntau, vd4_c, vd4_w;
    __m128d vd2_ntau, vd2_c, vd2_w;

    vd4_ntau = _mm256_set1_pd(ntau);
    vd2_ntau = _mm256_castpd256_pd128(vd4_ntau);

    for(i = 0; i <= m - 4; i += 4)
    {
        vd4_c = _mm256_loadu_pd((const doublereal *)&c1[i]);
        vd4_w = _mm256_add_pd(_mm256_loadu_pd((const doublereal *)&work[i]), vd4_c);
        _mm256_storeu_pd((doublereal *)&work[i], vd4_w);
        _mm256_storeu_pd((doublereal *)&c1[i], _mm256_fmadd_pd(vd4_ntau, vd4_w, vd4_c));
    }
    if(i <= m - 2)
    {
        vd2_c = _mm_loadu_pd((const doublereal *)&c1[i]);
        vd2_w = _mm_add_pd(_mm_loadu_pd((const doublereal *)&work[i]), vd2_c);
        _mm_storeu_pd((doublereal *)&work[i], vd2_w);
        _mm_storeu_pd((doublereal *)&c1[i], _mm_fmadd_pd(vd2_ntau, vd2_w, vd2_c));
        i += 2;
    }
    if(i < m)
    {
        vd2_c = _mm_load_sd((const doublereal *)&c1[i]);
        vd2_w = _mm_add_sd(_mm_load_sd((const doublereal *)&work[i]), vd2_c);
        _mm_storel_pd((doublereal *)&work[i], vd2_w);
        _mm_storel_pd((doublereal *)&c1[i], _mm_fmadd_sd(vd2_ntau, vd2_w, vd2_c));
    }
}

#endif
