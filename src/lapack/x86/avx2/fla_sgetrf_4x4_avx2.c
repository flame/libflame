/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 ******************************************************************************/

/*! @file fla_sgetrf_4x4_avx2.c
 *  @brief LU with partial pivoting for exactly 4x4 single-precision matrices.
 */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* Absolute value of the lowest single-precision lane. */
#define XMM_ABS_SS(v) _mm_andnot_ps(_mm_set_ss(-0.0f), (v))

/* Extract scalar from lane 0 to 3 of an __m128 register. */
#define XMM_ELEM(v, lane)                                                \
    _mm_cvtss_f32((lane) == 0 ? (v)                                      \
                            : _mm_permute_ps((v), (lane) == 1   ? 0x55 \
                                                    : (lane) == 2 ? 0xAA \
                                                                : 0xFF))

/* Permute all four column registers to swap matrix rows in-place. */
#define SWAP_COL4(c0, c1, c2, c3, imm)      \
    do                                      \
    {                                       \
        (c0) = _mm_permute_ps((c0), (imm)); \
        (c1) = _mm_permute_ps((c1), (imm)); \
        (c2) = _mm_permute_ps((c2), (imm)); \
        (c3) = _mm_permute_ps((c3), (imm)); \
    } while(0)

/* c - a*b; AVX2-only dispatch must not assume FMA. */
#define XMM_NMSUB(a, b, c) _mm_sub_ps((c), _mm_mul_ps((a), (b)))

aocl_int64_t fla_sgetrf_4x4_avx2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda,
                                aocl_int_t *ipiv, aocl_int64_t *info)
{
    aocl_int64_t lda_val = *lda;
    real av0, av1, av2, av3; /* Absolute pivot candidates. */
    real mx, piv;            /* Max pivot value and current diagonal. */
    real a32, a23, a33, l32, s33; /* Stage-3 Schur complement scalars. */
    aocl_int64_t pi;         /* Pivot row index within the active panel. */
    __m128 c0, c1, c2, c3;   /* Matrix columns (one row per lane). */
    __m128 l, u;             /* Elimination multiplier and pivot row broadcast. */

    *info = 0;
    (void)m;
    (void)n;

    // Load four rows as column vectors; each lane holds one matrix row.
    c0 = _mm_loadu_ps(a);
    c1 = _mm_loadu_ps(a + lda_val);
    c2 = _mm_loadu_ps(a + 2 * lda_val);
    c3 = _mm_loadu_ps(a + 3 * lda_val);

    // Stage 1: pivot on column 0, swap rows, scale L and update trailing cols.
    av0 = XMM_ELEM(XMM_ABS_SS(c0), 0);
    av1 = XMM_ELEM(XMM_ABS_SS(_mm_permute_ps(c0, 0x55)), 0);
    av2 = XMM_ELEM(XMM_ABS_SS(_mm_permute_ps(c0, 0xAA)), 0);
    av3 = XMM_ELEM(XMM_ABS_SS(_mm_permute_ps(c0, 0xFF)), 0);
    pi = 0;
    mx = av0;
    if(av1 > mx)
    {
        mx = av1;
        pi = 1;
    }
    if(av2 > mx)
    {
        mx = av2;
        pi = 2;
    }
    if(av3 > mx)
    {
        mx = av3;
        pi = 3;
    }
    ipiv[0] = (aocl_int_t)(pi + 1);
    if(mx != 0.0f)
    {
        // Swap rows 0 and pi via column permutes (no scalar memory traffic).
        if(pi == 1)
        {
            SWAP_COL4(c0, c1, c2, c3, 0xE1);
        }
        else if(pi == 2)
        {
            SWAP_COL4(c0, c1, c2, c3, 0xC6);
        }
        else if(pi == 3)
        {
            SWAP_COL4(c0, c1, c2, c3, 0x27);
        }
        // L[:,0] = A[:,0] / A[0,0]; U[0,:] unchanged in c0.
        piv = _mm_cvtss_f32(c0);
        l = _mm_mul_ps(c0, _mm_set1_ps(1.0f / piv));
        c0 = _mm_move_ss(l, c0);
        // Trailing update: col j -= L[:,0] * U[0,j] for j = 1..3.
        u = _mm_permute_ps(c1, 0x00);
        c1 = _mm_move_ss(XMM_NMSUB(c0, u, c1), c1);
        u = _mm_permute_ps(c2, 0x00);
        c2 = _mm_move_ss(XMM_NMSUB(c0, u, c2), c2);
        u = _mm_permute_ps(c3, 0x00);
        c3 = _mm_move_ss(XMM_NMSUB(c0, u, c3), c3);
    }
    else
    {
        if(*info == 0)
            *info = 1;
    }

    // Stage 2: pivot on column 1 (rows 1..3), eliminate below A[1,1].
    av1 = XMM_ELEM(XMM_ABS_SS(_mm_permute_ps(c1, 0x55)), 0);
    av2 = XMM_ELEM(XMM_ABS_SS(_mm_permute_ps(c1, 0xAA)), 0);
    av3 = XMM_ELEM(XMM_ABS_SS(_mm_permute_ps(c1, 0xFF)), 0);
    pi = 1;
    mx = av1;
    if(av2 > mx)
    {
        mx = av2;
        pi = 2;
    }
    if(av3 > mx)
    {
        mx = av3;
        pi = 3;
    }
    ipiv[1] = (aocl_int_t)(pi + 1);
    if(mx != 0.0f)
    {
        if(pi == 2)
        {
            SWAP_COL4(c0, c1, c2, c3, 0xD8);
        }
        else if(pi == 3)
        {
            SWAP_COL4(c0, c1, c2, c3, 0x6C);
        }
        // Scale L[1:3,1] and update cols 2..3; preserve U[0:1,1] in c1.
        piv = XMM_ELEM(_mm_permute_ps(c1, 0x55), 0);
        l = _mm_mul_ps(c1, _mm_set1_ps(1.0f / piv));
        c1 = _mm_blend_ps(c1, l, 0xC);
        u = _mm_permute_ps(c2, 0x55);
        c2 = _mm_blend_ps(c2, XMM_NMSUB(c1, u, c2), 0xC);
        u = _mm_permute_ps(c3, 0x55);
        c3 = _mm_blend_ps(c3, XMM_NMSUB(c1, u, c3), 0xC);
    }
    else
    {
        if(*info == 0)
            *info = 2;
    }

    // Stage 3: pivot on column 2 (rows 2..3), finish 2x2 trailing block.
    av2 = XMM_ELEM(XMM_ABS_SS(_mm_permute_ps(c2, 0xAA)), 0);
    av3 = XMM_ELEM(XMM_ABS_SS(_mm_permute_ps(c2, 0xFF)), 0);
    pi = (av3 > av2) ? 3 : 2;
    mx = (pi == 2) ? av2 : av3;
    ipiv[2] = (aocl_int_t)(pi + 1);
    if(mx != 0.0f)
    {
        if(pi != 2)
        {
            SWAP_COL4(c0, c1, c2, c3, 0xB4);
        }
        // Schur update on the bottom-right 2x2: L[3,2] and U[3,3].
        piv = XMM_ELEM(_mm_permute_ps(c2, 0xAA), 0);
        a32 = XMM_ELEM(_mm_permute_ps(c2, 0xFF), 0);
        a23 = XMM_ELEM(_mm_permute_ps(c3, 0xAA), 0);
        a33 = XMM_ELEM(_mm_permute_ps(c3, 0xFF), 0);
        l32 = a32 / piv;
        s33 = a33 - l32 * a23;
        c2 = _mm_insert_ps(c2, _mm_set_ss(l32), 0x30);
        c3 = _mm_insert_ps(c3, _mm_set_ss(s33), 0x30);
    }
    else
    {
        if(*info == 0)
            *info = 3;
    }

    // Stage 4: record final pivot; flag singularity if A[3,3] is zero.
    a33 = XMM_ELEM(_mm_permute_ps(c3, 0xFF), 0);
    ipiv[3] = 4;
    if((a33 == 0.0f) && (*info == 0))
        *info = 4;

    // Store column registers back to column-major matrix layout.
    _mm_storeu_ps(a, c0);
    _mm_storeu_ps(a + lda_val, c1);
    _mm_storeu_ps(a + 2 * lda_val, c2);
    _mm_storeu_ps(a + 3 * lda_val, c3);

    return *info;
}

#undef XMM_ABS_SS
#undef XMM_ELEM
#undef SWAP_COL4
#undef XMM_NMSUB

#endif
