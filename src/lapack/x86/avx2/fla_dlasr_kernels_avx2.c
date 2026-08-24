/******************************************************************************
 * * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* Kernels for side='L', pivot='V'.

   The reference nest runs one rotation across all n columns before the next,
   so its working set is one cache line per column -- 64n bytes, outgrowing L1
   at n = 512 -- and the row two consecutive rotations share is written and
   read a full pass apart. That pass also collapses onto a handful of L1 sets
   when a_dim1 is a multiple of 64 doubles (the power-of-2 leading dimension
   cliff, 10-40x).

   These kernels invert the nest and run the whole rotation sequence down a
   block of columns: the working set is the block's 8 lines whatever n is, the
   shared row stays in register, and a fused pair loads 2 rows and stores 2
   rather than 3 and 3. Runtime is then flat in a_dim1, so no path here
   branches on it.

   A 4-column row slice is one ymm, assembled by hand from 8-byte loads --
   letting the compiler vectorize the row pair costs a cross-lane transpose,
   and hardware gather/scatter measured slower. Blocks are two such lanes so
   each rotation's broadcasts, now inside the column loop, amortize over 8
   columns instead of 4.

   fla_dlasr_kernels_avx512.c has the zmm version and takes over from
   n >= FLA_DLASR_L_SIMD_AVX512_THRESH_N; only its tail columns reach here. */

/* The reference skips a rotation that is exactly the identity. For finite data
   applying it anyway is exact, but it would turn -0.0 into +0.0 and 0*Inf into
   NaN, so the test is kept to stay bit-for-bit with dlasr.c. */
static inline logical fla_dlasr_is_identity(doublereal ct, doublereal st)
{
    return ct == 1. && st == 0.;
}

/* Gather one row across the 4 columns p, p+st, p+2*st, p+3*st. */
static inline __m256d fla_dlasr_row4_avx2(const doublereal *p, aocl_int64_t st)
{
    __m128d x0 = _mm_loadh_pd(_mm_load_sd(p), p + st);
    __m128d x1 = _mm_loadh_pd(_mm_load_sd(p + 2 * st), p + 3 * st);

    return _mm256_set_m128d(x1, x0);
}

/* Scatter two adjacent rows, r0 at p and r1 at p+1, across the same 4 columns.
   The pair is contiguous within a column, so one 16-byte store per column. */
static inline void fla_dlasr_store4x2_avx2(doublereal *p, aocl_int64_t st, __m256d r0, __m256d r1)
{
    __m256d lo = _mm256_unpacklo_pd(r0, r1);
    __m256d hi = _mm256_unpackhi_pd(r0, r1);

    _mm_storeu_pd(p, _mm256_castpd256_pd128(lo));
    _mm_storeu_pd(p + st, _mm256_castpd256_pd128(hi));
    _mm_storeu_pd(p + 2 * st, _mm256_extractf128_pd(lo, 1));
    _mm_storeu_pd(p + 3 * st, _mm256_extractf128_pd(hi, 1));
}

/* Scatter a single row across the same 4 columns. */
static inline void fla_dlasr_store4_avx2(doublereal *p, aocl_int64_t st, __m256d v)
{
    __m128d z0 = _mm256_castpd256_pd128(v);
    __m128d z1 = _mm256_extractf128_pd(v, 1);

    _mm_storel_pd(p, z0);
    _mm_storeh_pd(p + st, z0);
    _mm_storel_pd(p + 2 * st, z1);
    _mm_storeh_pd(p + 3 * st, z1);
}

/* Fused rotations j and j+1 (direct='F') on the 4-column lane at p = col + j,
   where col[j + k*a_dim1] is row j of column k. Rows j+1 and j+2 are read,
   rows j and j+1 written; row j+2 is returned so the next pair starts from
   register. A skipped rotation still stores -- the rows it covers are stale. */
static inline __m256d fla_dlasr_pair4_fwd_avx2(doublereal *p, aocl_int64_t st, __m256d t0,
                                               __m256d vct0, __m256d vst0, __m256d vct1,
                                               __m256d vst1, logical id0, logical id1)
{
    __m256d t1 = fla_dlasr_row4_avx2(p + 1, st);
    __m256d t2 = fla_dlasr_row4_avx2(p + 2, st);
    __m256d r0 = t0, m1 = t1, r1;

    if(!id0)
    {
        r0 = _mm256_fmadd_pd(vst0, t1, _mm256_mul_pd(vct0, t0));
        m1 = _mm256_fmsub_pd(vct0, t1, _mm256_mul_pd(vst0, t0));
    }
    r1 = m1;
    if(!id1)
    {
        r1 = _mm256_fmadd_pd(vst1, t2, _mm256_mul_pd(vct1, m1));
        t2 = _mm256_fmsub_pd(vct1, t2, _mm256_mul_pd(vst1, m1));
    }

    fla_dlasr_store4x2_avx2(p, st, r0, r1);
    return t2;
}

/* Lone rotation m-1 closing an even-m forward sweep, at p = col + m - 1. */
static inline void fla_dlasr_last4_fwd_avx2(doublereal *p, aocl_int64_t st, __m256d t0, __m256d vct,
                                            __m256d vst, logical id)
{
    __m256d t1 = fla_dlasr_row4_avx2(p + 1, st);
    __m256d r0 = t0, r1 = t1;

    if(!id)
    {
        r0 = _mm256_fmadd_pd(vst, t1, _mm256_mul_pd(vct, t0));
        r1 = _mm256_fmsub_pd(vct, t1, _mm256_mul_pd(vst, t0));
    }

    fla_dlasr_store4x2_avx2(p, st, r0, r1);
}

/* Backward counterpart: rotations j and j-1 at p = col + j read rows j and
   j-1, write rows j+1 and j, and return the value left in row j-1. */
static inline __m256d fla_dlasr_pair4_bwd_avx2(doublereal *p, aocl_int64_t st, __m256d t2,
                                               __m256d vct1, __m256d vst1, __m256d vct0,
                                               __m256d vst0, logical id1, logical id0)
{
    __m256d t1 = fla_dlasr_row4_avx2(p, st);
    __m256d t0 = fla_dlasr_row4_avx2(p - 1, st);
    __m256d r2 = t2, m1 = t1, r1;

    if(!id1)
    {
        r2 = _mm256_fmsub_pd(vct1, t2, _mm256_mul_pd(vst1, t1));
        m1 = _mm256_fmadd_pd(vst1, t2, _mm256_mul_pd(vct1, t1));
    }
    r1 = m1;
    if(!id0)
    {
        r1 = _mm256_fmsub_pd(vct0, m1, _mm256_mul_pd(vst0, t0));
        t0 = _mm256_fmadd_pd(vst0, m1, _mm256_mul_pd(vct0, t0));
    }

    fla_dlasr_store4x2_avx2(p, st, r1, r2);
    return t0;
}

/* Lone rotation 1 closing an even-m backward sweep, at p = col + 1. */
static inline void fla_dlasr_last4_bwd_avx2(doublereal *p, aocl_int64_t st, __m256d t2, __m256d vct,
                                            __m256d vst, logical id)
{
    __m256d t1 = fla_dlasr_row4_avx2(p, st);
    __m256d r2 = t2, r1 = t1;

    if(!id)
    {
        r2 = _mm256_fmsub_pd(vct, t2, _mm256_mul_pd(vst, t1));
        r1 = _mm256_fmadd_pd(vst, t2, _mm256_mul_pd(vct, t1));
    }

    fla_dlasr_store4x2_avx2(p, st, r1, r2);
}

/* Run rotations 1..m-1 in increasing order (direct='F') down one block of 8
   columns, as two ymm lanes sharing one set of broadcasts. */
static void fla_dlasr_col8_fwd_avx2(doublereal *col, aocl_int64_t m, aocl_int64_t a_dim1,
                                    const doublereal *c__, const doublereal *s)
{
    doublereal *colb = col + 4 * a_dim1;
    __m256d t0a = fla_dlasr_row4_avx2(col + 1, a_dim1);
    __m256d t0b = fla_dlasr_row4_avx2(colb + 1, a_dim1);
    aocl_int64_t j;

    for(j = 1; j < m - 1; j += 2)
    {
        const logical id0 = fla_dlasr_is_identity(c__[j], s[j]);
        const logical id1 = fla_dlasr_is_identity(c__[j + 1], s[j + 1]);
        const __m256d vct0 = _mm256_set1_pd(c__[j]);
        const __m256d vst0 = _mm256_set1_pd(s[j]);
        const __m256d vct1 = _mm256_set1_pd(c__[j + 1]);
        const __m256d vst1 = _mm256_set1_pd(s[j + 1]);

        t0a = fla_dlasr_pair4_fwd_avx2(col + j, a_dim1, t0a, vct0, vst0, vct1, vst1, id0, id1);
        t0b = fla_dlasr_pair4_fwd_avx2(colb + j, a_dim1, t0b, vct0, vst0, vct1, vst1, id0, id1);
    }

    if(j == m - 1)
    {
        const logical id = fla_dlasr_is_identity(c__[j], s[j]);
        const __m256d vct = _mm256_set1_pd(c__[j]);
        const __m256d vst = _mm256_set1_pd(s[j]);

        fla_dlasr_last4_fwd_avx2(col + j, a_dim1, t0a, vct, vst, id);
        fla_dlasr_last4_fwd_avx2(colb + j, a_dim1, t0b, vct, vst, id);
    }
    else
    {
        fla_dlasr_store4_avx2(col + m, a_dim1, t0a);
        fla_dlasr_store4_avx2(colb + m, a_dim1, t0b);
    }
}

/* As above for direct='B': rotations m-1..1 in decreasing order. */
static void fla_dlasr_col8_bwd_avx2(doublereal *col, aocl_int64_t m, aocl_int64_t a_dim1,
                                    const doublereal *c__, const doublereal *s)
{
    doublereal *colb = col + 4 * a_dim1;
    __m256d t2a = fla_dlasr_row4_avx2(col + m, a_dim1);
    __m256d t2b = fla_dlasr_row4_avx2(colb + m, a_dim1);
    aocl_int64_t j;

    for(j = m - 1; j >= 2; j -= 2)
    {
        const logical id1 = fla_dlasr_is_identity(c__[j], s[j]);
        const logical id0 = fla_dlasr_is_identity(c__[j - 1], s[j - 1]);
        const __m256d vct1 = _mm256_set1_pd(c__[j]);
        const __m256d vst1 = _mm256_set1_pd(s[j]);
        const __m256d vct0 = _mm256_set1_pd(c__[j - 1]);
        const __m256d vst0 = _mm256_set1_pd(s[j - 1]);

        t2a = fla_dlasr_pair4_bwd_avx2(col + j, a_dim1, t2a, vct1, vst1, vct0, vst0, id1, id0);
        t2b = fla_dlasr_pair4_bwd_avx2(colb + j, a_dim1, t2b, vct1, vst1, vct0, vst0, id1, id0);
    }

    if(j == 1)
    {
        const logical id = fla_dlasr_is_identity(c__[1], s[1]);
        const __m256d vct = _mm256_set1_pd(c__[1]);
        const __m256d vst = _mm256_set1_pd(s[1]);

        fla_dlasr_last4_bwd_avx2(col + 1, a_dim1, t2a, vct, vst, id);
        fla_dlasr_last4_bwd_avx2(colb + 1, a_dim1, t2b, vct, vst, id);
    }
    else
    {
        fla_dlasr_store4_avx2(col + 1, a_dim1, t2a);
        fla_dlasr_store4_avx2(colb + 1, a_dim1, t2b);
    }
}

/* Single-lane form of fla_dlasr_col8_fwd_avx2 for a trailing block of 4. */
static void fla_dlasr_col4_fwd_avx2(doublereal *col, aocl_int64_t m, aocl_int64_t a_dim1,
                                    const doublereal *c__, const doublereal *s)
{
    __m256d t0 = fla_dlasr_row4_avx2(col + 1, a_dim1);
    aocl_int64_t j;

    for(j = 1; j < m - 1; j += 2)
    {
        const logical id0 = fla_dlasr_is_identity(c__[j], s[j]);
        const logical id1 = fla_dlasr_is_identity(c__[j + 1], s[j + 1]);
        const __m256d vct0 = _mm256_set1_pd(c__[j]);
        const __m256d vst0 = _mm256_set1_pd(s[j]);
        const __m256d vct1 = _mm256_set1_pd(c__[j + 1]);
        const __m256d vst1 = _mm256_set1_pd(s[j + 1]);

        t0 = fla_dlasr_pair4_fwd_avx2(col + j, a_dim1, t0, vct0, vst0, vct1, vst1, id0, id1);
    }

    if(j == m - 1)
    {
        const logical id = fla_dlasr_is_identity(c__[j], s[j]);
        const __m256d vct = _mm256_set1_pd(c__[j]);
        const __m256d vst = _mm256_set1_pd(s[j]);

        fla_dlasr_last4_fwd_avx2(col + j, a_dim1, t0, vct, vst, id);
    }
    else
    {
        fla_dlasr_store4_avx2(col + m, a_dim1, t0);
    }
}

/* Single-lane form of fla_dlasr_col8_bwd_avx2 for a trailing block of 4. */
static void fla_dlasr_col4_bwd_avx2(doublereal *col, aocl_int64_t m, aocl_int64_t a_dim1,
                                    const doublereal *c__, const doublereal *s)
{
    __m256d t2 = fla_dlasr_row4_avx2(col + m, a_dim1);
    aocl_int64_t j;

    for(j = m - 1; j >= 2; j -= 2)
    {
        const logical id1 = fla_dlasr_is_identity(c__[j], s[j]);
        const logical id0 = fla_dlasr_is_identity(c__[j - 1], s[j - 1]);
        const __m256d vct1 = _mm256_set1_pd(c__[j]);
        const __m256d vst1 = _mm256_set1_pd(s[j]);
        const __m256d vct0 = _mm256_set1_pd(c__[j - 1]);
        const __m256d vst0 = _mm256_set1_pd(s[j - 1]);

        t2 = fla_dlasr_pair4_bwd_avx2(col + j, a_dim1, t2, vct1, vst1, vct0, vst0, id1, id0);
    }

    if(j == 1)
    {
        const logical id = fla_dlasr_is_identity(c__[1], s[1]);
        const __m256d vct = _mm256_set1_pd(c__[1]);
        const __m256d vst = _mm256_set1_pd(s[1]);

        fla_dlasr_last4_bwd_avx2(col + 1, a_dim1, t2, vct, vst, id);
    }
    else
    {
        fla_dlasr_store4_avx2(col + 1, a_dim1, t2);
    }
}

/* Scalar form of fla_dlasr_col4_fwd_avx2 for the n % 4 trailing columns. */
static void fla_dlasr_col1_fwd_avx2(doublereal *col, aocl_int64_t m, const doublereal *c__,
                                    const doublereal *s)
{
    doublereal t0 = col[1];
    aocl_int64_t j;

    for(j = 1; j < m - 1; j += 2)
    {
        doublereal ct0 = c__[j], st0 = s[j], ct1 = c__[j + 1], st1 = s[j + 1];
        doublereal t1 = col[j + 1], t2 = col[j + 2];
        doublereal m1 = t1;

        if(!fla_dlasr_is_identity(ct0, st0))
        {
            m1 = ct0 * t1 - st0 * t0;
            col[j] = st0 * t1 + ct0 * t0;
        }
        else
        {
            col[j] = t0;
        }
        if(!fla_dlasr_is_identity(ct1, st1))
        {
            col[j + 1] = st1 * t2 + ct1 * m1;
            t0 = ct1 * t2 - st1 * m1;
        }
        else
        {
            col[j + 1] = m1;
            t0 = t2;
        }
    }

    if(j == m - 1)
    {
        doublereal ct = c__[j], st = s[j], t1 = col[m];

        if(!fla_dlasr_is_identity(ct, st))
        {
            col[j] = st * t1 + ct * t0;
            col[m] = ct * t1 - st * t0;
        }
        else
        {
            col[j] = t0;
        }
    }
    else
    {
        col[m] = t0;
    }
}

/* Scalar form of fla_dlasr_col4_bwd_avx2 for the n % 4 trailing columns. */
static void fla_dlasr_col1_bwd_avx2(doublereal *col, aocl_int64_t m, const doublereal *c__,
                                    const doublereal *s)
{
    doublereal t2 = col[m];
    aocl_int64_t j;

    for(j = m - 1; j >= 2; j -= 2)
    {
        doublereal ct1 = c__[j], st1 = s[j], ct0 = c__[j - 1], st0 = s[j - 1];
        doublereal t1 = col[j], t0 = col[j - 1];
        doublereal m1 = t1;

        if(!fla_dlasr_is_identity(ct1, st1))
        {
            m1 = st1 * t2 + ct1 * t1;
            col[j + 1] = ct1 * t2 - st1 * t1;
        }
        else
        {
            col[j + 1] = t2;
        }
        if(!fla_dlasr_is_identity(ct0, st0))
        {
            col[j] = ct0 * m1 - st0 * t0;
            t2 = st0 * m1 + ct0 * t0;
        }
        else
        {
            col[j] = m1;
            t2 = t0;
        }
    }

    if(j == 1)
    {
        doublereal ct = c__[1], st = s[1], t1 = col[1];

        if(!fla_dlasr_is_identity(ct, st))
        {
            col[2] = ct * t2 - st * t1;
            col[1] = st * t2 + ct * t1;
        }
        else
        {
            col[2] = t2;
        }
    }
    else
    {
        col[1] = t2;
    }
}

/* Apply a sequence of plane rotations from the left with a variable pivot
 * (side='L', pivot='V'); forward != 0 selects direct='F', else direct='B'.
 * c__, s and a are the f2c-adjusted one-based pointers, so the callers'
 * --c__, --s and a -= a_offset must already have been applied.
 * */
void fla_dlasr_left_pivotv_avx2(logical forward, aocl_int64_t m, aocl_int64_t n, doublereal *c__,
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
            fla_dlasr_col8_fwd_avx2(a + i * a_dim1, m, a_dim1, c__, s);
        }
        if(i + 3 <= n)
        {
            fla_dlasr_col4_fwd_avx2(a + i * a_dim1, m, a_dim1, c__, s);
            i += 4;
        }
        for(; i <= n; ++i)
        {
            fla_dlasr_col1_fwd_avx2(a + i * a_dim1, m, c__, s);
        }
    }
    else
    {
        for(; i + 7 <= n; i += 8)
        {
            fla_dlasr_col8_bwd_avx2(a + i * a_dim1, m, a_dim1, c__, s);
        }
        if(i + 3 <= n)
        {
            fla_dlasr_col4_bwd_avx2(a + i * a_dim1, m, a_dim1, c__, s);
            i += 4;
        }
        for(; i <= n; ++i)
        {
            fla_dlasr_col1_bwd_avx2(a + i * a_dim1, m, c__, s);
        }
    }
    return;
}
#endif
