/******************************************************************************
 * * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* 512-bit kernels for side='L', pivot='V'. See fla_dlasr_kernels_avx2.c for
   why the nest is inverted and the row slices are assembled by hand.

   8 columns are one zmm, gathered with 8-byte loads and 3 inserts. The two rows
   that a fused pair retires are adjacent within a column, so they leave as a single
   16-byte store per column. Trailing columns go to the 256-bit path. */

/* Gather one row across the 8 columns p, p+st, ... p+7*st, in column order. */
static inline __m512d fla_dlasr_row8_avx512(const doublereal *p, aocl_int64_t st)
{
    __m128d x0 = _mm_loadh_pd(_mm_load_sd(p), p + st);
    __m128d x1 = _mm_loadh_pd(_mm_load_sd(p + 2 * st), p + 3 * st);
    __m128d x2 = _mm_loadh_pd(_mm_load_sd(p + 4 * st), p + 5 * st);
    __m128d x3 = _mm_loadh_pd(_mm_load_sd(p + 6 * st), p + 7 * st);
    __m512d v = _mm512_castpd128_pd512(x0);

    v = _mm512_insertf64x2(v, x1, 1);
    v = _mm512_insertf64x2(v, x2, 2);
    v = _mm512_insertf64x2(v, x3, 3);
    return v;
}

/* Scatter two adjacent rows, r0 at p and r1 at p+1, across the same 8 columns.
   The pair is contiguous within a column, so one 16-byte store per column. */
static inline void fla_dlasr_store8x2_avx512(doublereal *p, aocl_int64_t st, __m512d r0, __m512d r1)
{
    __m512d lo = _mm512_unpacklo_pd(r0, r1);
    __m512d hi = _mm512_unpackhi_pd(r0, r1);

    _mm_storeu_pd(p, _mm512_castpd512_pd128(lo));
    _mm_storeu_pd(p + st, _mm512_castpd512_pd128(hi));
    _mm_storeu_pd(p + 2 * st, _mm512_extractf64x2_pd(lo, 1));
    _mm_storeu_pd(p + 3 * st, _mm512_extractf64x2_pd(hi, 1));
    _mm_storeu_pd(p + 4 * st, _mm512_extractf64x2_pd(lo, 2));
    _mm_storeu_pd(p + 5 * st, _mm512_extractf64x2_pd(hi, 2));
    _mm_storeu_pd(p + 6 * st, _mm512_extractf64x2_pd(lo, 3));
    _mm_storeu_pd(p + 7 * st, _mm512_extractf64x2_pd(hi, 3));
}

/* Scatter a single row across the same 8 columns. */
static inline void fla_dlasr_store8_avx512(doublereal *p, aocl_int64_t st, __m512d v)
{
    __m128d z0 = _mm512_castpd512_pd128(v);
    __m128d z1 = _mm512_extractf64x2_pd(v, 1);
    __m128d z2 = _mm512_extractf64x2_pd(v, 2);
    __m128d z3 = _mm512_extractf64x2_pd(v, 3);

    _mm_storel_pd(p, z0);
    _mm_storeh_pd(p + st, z0);
    _mm_storel_pd(p + 2 * st, z1);
    _mm_storeh_pd(p + 3 * st, z1);
    _mm_storel_pd(p + 4 * st, z2);
    _mm_storeh_pd(p + 5 * st, z2);
    _mm_storel_pd(p + 6 * st, z3);
    _mm_storeh_pd(p + 7 * st, z3);
}

/* The reference skips a rotation that is exactly the identity. For finite data
   applying it anyway is exact, but it would turn -0.0 into +0.0 and 0*Inf into
   NaN, so the test is kept to stay bit-for-bit with dlasr.c. */
static inline logical fla_dlasr_is_identity(doublereal ct, doublereal st)
{
    return ct == 1. && st == 0.;
}

/* Rotations 1..m-1 in increasing order (direct='F') down one block of 8
   columns at col, where col[j + k*a_dim1] is row j of column k. Two rotations
   retire per iteration; the row they share and the row the next iteration
   starts from stay in register. A skipped rotation still stores -- the rows it
   covers are stale. */
static void fla_dlasr_col8_fwd_avx512(doublereal *col, aocl_int64_t m, aocl_int64_t a_dim1,
                                      const doublereal *c__, const doublereal *s)
{
    __m512d t0 = fla_dlasr_row8_avx512(col + 1, a_dim1);
    aocl_int64_t j;

    for(j = 1; j < m - 1; j += 2)
    {
        const __m512d vct0 = _mm512_set1_pd(c__[j]);
        const __m512d vst0 = _mm512_set1_pd(s[j]);
        const __m512d vct1 = _mm512_set1_pd(c__[j + 1]);
        const __m512d vst1 = _mm512_set1_pd(s[j + 1]);
        __m512d t1 = fla_dlasr_row8_avx512(col + j + 1, a_dim1);
        __m512d t2 = fla_dlasr_row8_avx512(col + j + 2, a_dim1);
        __m512d r0 = t0, m1 = t1, r1;

        if(!fla_dlasr_is_identity(c__[j], s[j]))
        {
            r0 = _mm512_fmadd_pd(vst0, t1, _mm512_mul_pd(vct0, t0));
            m1 = _mm512_fmsub_pd(vct0, t1, _mm512_mul_pd(vst0, t0));
        }
        r1 = m1;
        /* row j+2 is still live for the next pair, so keep it in register */
        t0 = t2;
        if(!fla_dlasr_is_identity(c__[j + 1], s[j + 1]))
        {
            r1 = _mm512_fmadd_pd(vst1, t2, _mm512_mul_pd(vct1, m1));
            t0 = _mm512_fmsub_pd(vct1, t2, _mm512_mul_pd(vst1, m1));
        }

        fla_dlasr_store8x2_avx512(col + j, a_dim1, r0, r1);
    }

    if(j == m - 1)
    {
        /* even m: one unfused rotation closes the sweep */
        const __m512d vct = _mm512_set1_pd(c__[j]);
        const __m512d vst = _mm512_set1_pd(s[j]);
        __m512d t1 = fla_dlasr_row8_avx512(col + m, a_dim1);
        __m512d r0 = t0, r1 = t1;

        if(!fla_dlasr_is_identity(c__[j], s[j]))
        {
            r0 = _mm512_fmadd_pd(vst, t1, _mm512_mul_pd(vct, t0));
            r1 = _mm512_fmsub_pd(vct, t1, _mm512_mul_pd(vst, t0));
        }

        fla_dlasr_store8x2_avx512(col + j, a_dim1, r0, r1);
    }
    else
    {
        fla_dlasr_store8_avx512(col + m, a_dim1, t0);
    }
}

/* As above for direct='B': rotations m-1..1 in decreasing order. */
static void fla_dlasr_col8_bwd_avx512(doublereal *col, aocl_int64_t m, aocl_int64_t a_dim1,
                                      const doublereal *c__, const doublereal *s)
{
    __m512d t2 = fla_dlasr_row8_avx512(col + m, a_dim1);
    aocl_int64_t j;

    for(j = m - 1; j >= 2; j -= 2)
    {
        const __m512d vct1 = _mm512_set1_pd(c__[j]);
        const __m512d vst1 = _mm512_set1_pd(s[j]);
        const __m512d vct0 = _mm512_set1_pd(c__[j - 1]);
        const __m512d vst0 = _mm512_set1_pd(s[j - 1]);
        __m512d t1 = fla_dlasr_row8_avx512(col + j, a_dim1);
        __m512d t0 = fla_dlasr_row8_avx512(col + j - 1, a_dim1);
        __m512d r2 = t2, m1 = t1, r1;

        if(!fla_dlasr_is_identity(c__[j], s[j]))
        {
            r2 = _mm512_fmsub_pd(vct1, t2, _mm512_mul_pd(vst1, t1));
            m1 = _mm512_fmadd_pd(vst1, t2, _mm512_mul_pd(vct1, t1));
        }
        r1 = m1;
        /* row j-1 is still live for the next pair, so keep it in register */
        t2 = t0;
        if(!fla_dlasr_is_identity(c__[j - 1], s[j - 1]))
        {
            r1 = _mm512_fmsub_pd(vct0, m1, _mm512_mul_pd(vst0, t0));
            t2 = _mm512_fmadd_pd(vst0, m1, _mm512_mul_pd(vct0, t0));
        }

        fla_dlasr_store8x2_avx512(col + j, a_dim1, r1, r2);
    }

    if(j == 1)
    {
        const __m512d vct = _mm512_set1_pd(c__[1]);
        const __m512d vst = _mm512_set1_pd(s[1]);
        __m512d t1 = fla_dlasr_row8_avx512(col + 1, a_dim1);
        __m512d r2 = t2, r1 = t1;

        if(!fla_dlasr_is_identity(c__[1], s[1]))
        {
            r2 = _mm512_fmsub_pd(vct, t2, _mm512_mul_pd(vst, t1));
            r1 = _mm512_fmadd_pd(vst, t2, _mm512_mul_pd(vct, t1));
        }

        fla_dlasr_store8x2_avx512(col + 1, a_dim1, r1, r2);
    }
    else
    {
        fla_dlasr_store8_avx512(col + 1, a_dim1, t2);
    }
}

/* Apply a sequence of plane rotations from the left with a variable pivot
 * (side='L', pivot='V'); forward != 0 selects direct='F', else direct='B'.
 * c__, s and a are the f2c-adjusted one-based pointers, so the callers'
 * --c__, --s and a -= a_offset must already have been applied.
 * */
void fla_dlasr_left_pivotv_avx512(logical forward, aocl_int64_t m, aocl_int64_t n, doublereal *c__,
                                  doublereal *s, doublereal *a, aocl_int64_t a_dim1)
{
    aocl_int64_t i = 1;

    if(m < 2)
    {
        return;
    }

    if(forward)
    {
        for(; i + 7 <= n; i += 8)
        {
            fla_dlasr_col8_fwd_avx512(a + i * a_dim1, m, a_dim1, c__, s);
        }
    }
    else
    {
        for(; i + 7 <= n; i += 8)
        {
            fla_dlasr_col8_bwd_avx512(a + i * a_dim1, m, a_dim1, c__, s);
        }
    }

    if(i <= n)
    {
        fla_dlasr_left_pivotv_avx2(forward, m, n - i + 1, c__, s, a + (i - 1) * a_dim1, a_dim1);
    }
    return;
}
#endif
