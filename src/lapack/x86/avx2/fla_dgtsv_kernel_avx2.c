/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/*
 * AVX2 optimized DGTSV.
 *
 * Solves A * X = B where A is an n-by-n tridiagonal matrix, using Gaussian
 * elimination with partial pivoting (identical algorithm to reference LAPACK
 * DGTSV).
 *
 * The tridiagonal factorization (arrays dl, d, du) carries a strict data
 * dependency from one row to the next, so it is computed with scalar code.
 * The update of the right-hand-side matrix B, however, is independent across
 * the NRHS columns, so it is vectorized across columns.  Because B is stored
 * column-major, elements of a fixed row across consecutive columns are strided
 * by ldb; AVX2 gather instructions are used to load them 4 columns at a time.
 * AVX2 has no scatter instruction, so results are written back with scalar
 * stores.
 */

/* Store the 4 packed doubles of VEC into b[BASE], b[BASE+ldb], b[BASE+2*ldb],
   b[BASE+3*ldb] (column-major write-back, emulates the missing AVX2 scatter). */
#define FLA_DGTSV_STORE4(BASE, VEC)  \
    do                               \
    {                                \
        double _t[4];                \
        _mm256_storeu_pd(_t, (VEC)); \
        (BASE)[0] = _t[0];           \
        (BASE)[b_dim1] = _t[1];      \
        (BASE)[2 * b_dim1] = _t[2];  \
        (BASE)[3 * b_dim1] = _t[3];  \
    } while(0)

int fla_dgtsv_kernel_avx2(aocl_int64_t *n, aocl_int64_t *nrhs, doublereal *dl, doublereal *d__,
                          doublereal *du, doublereal *b, aocl_int64_t *ldb, aocl_int64_t *info)
{
    aocl_int64_t b_dim1, b_offset, i__, j, nn, nr;
    doublereal fact, temp;
    __m256i vidx;

    *info = 0;
    nn = *n;
    nr = *nrhs;

    /* Parameter adjustments to use 1-based indexing (as in reference LAPACK) */
    --dl;
    --d__;
    --du;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    if(nn == 0)
    {
        return 0;
    }

    /* Column-offset index vector: { 0, ldb, 2*ldb, 3*ldb } (in elements). */
    vidx = _mm256_set_epi64x(3 * b_dim1, 2 * b_dim1, b_dim1, 0);

    /* Forward elimination with partial pivoting. */
    for(i__ = 1; i__ <= nn - 1; ++i__)
    {
        if(fabs(d__[i__]) >= fabs(dl[i__]))
        {
            /* No row interchange required */
            if(d__[i__] != 0.)
            {
                fact = dl[i__] / d__[i__];
                d__[i__ + 1] -= fact * du[i__];

                /* B(i+1,:) := B(i+1,:) - fact * B(i,:) */
                {
                    __m256d fact_vec = _mm256_set1_pd(fact);
                    doublereal *bi = &b[i__ + b_dim1];
                    doublereal *bi1 = &b[i__ + 1 + b_dim1];
                    for(j = 0; j + 4 <= nr; j += 4)
                    {
                        aocl_int64_t off = j * b_dim1;
                        __m256d vbi = _mm256_i64gather_pd(bi + off, vidx, 8);
                        __m256d vbi1 = _mm256_i64gather_pd(bi1 + off, vidx, 8);
                        vbi1 = _mm256_fnmadd_pd(fact_vec, vbi, vbi1);
                        FLA_DGTSV_STORE4(bi1 + off, vbi1);
                    }
                    for(; j < nr; ++j)
                    {
                        aocl_int64_t offset = (j + 1) * b_dim1;
                        b[i__ + 1 + offset] -= fact * b[i__ + offset];
                    }
                }
            }
            else
            {
                *info = i__;
                return 0;
            }
            if(i__ < nn - 1)
            {
                dl[i__] = 0.;
            }
        }
        else
        {
            /* Interchange rows I and I+1 */
            fact = d__[i__] / dl[i__];
            d__[i__] = dl[i__];
            temp = d__[i__ + 1];
            d__[i__ + 1] = du[i__] - fact * temp;
            if(i__ < nn - 1)
            {
                dl[i__] = du[i__ + 1];
                du[i__ + 1] = -fact * dl[i__];
            }
            du[i__] = temp;

            /* Swap rows I and I+1 of B and eliminate:
               new B(i,:)   = old B(i+1,:)
               new B(i+1,:) = old B(i,:) - fact * old B(i+1,:) */
            {
                __m256d fact_vec = _mm256_set1_pd(fact);
                doublereal *bi = &b[i__ + b_dim1];
                doublereal *bi1 = &b[i__ + 1 + b_dim1];
                for(j = 0; j + 4 <= nr; j += 4)
                {
                    aocl_int64_t off = j * b_dim1;
                    __m256d vbi = _mm256_i64gather_pd(bi + off, vidx, 8);
                    __m256d vbi1 = _mm256_i64gather_pd(bi1 + off, vidx, 8);
                    __m256d newbi1 = _mm256_fnmadd_pd(fact_vec, vbi1, vbi);
                    FLA_DGTSV_STORE4(bi + off, vbi1);
                    FLA_DGTSV_STORE4(bi1 + off, newbi1);
                }
                for(; j < nr; ++j)
                {
                    aocl_int64_t offset = (j + 1) * b_dim1;
                    temp = b[i__ + offset];
                    b[i__ + offset] = b[i__ + 1 + offset];
                    b[i__ + 1 + offset] = temp - fact * b[i__ + 1 + offset];
                }
            }
        }
    }

    if(d__[nn] == 0.)
    {
        *info = nn;
        return 0;
    }

    /* Back substitution with the upper triangular factor U.
       Sequential over rows, vectorized across NRHS columns. */
    for(j = 0; j + 4 <= nr; j += 4)
    {
        /* Base offset points at column j (0-based) => 1-based column (j+1). */
        aocl_int64_t off = b_dim1 + j * b_dim1;
        __m256d bip1, bip2, vi;

        /* Row n */
        vi = _mm256_i64gather_pd(&b[nn + off], vidx, 8);
        vi = _mm256_div_pd(vi, _mm256_set1_pd(d__[nn]));
        FLA_DGTSV_STORE4(&b[nn + off], vi);
        bip2 = vi;
        bip1 = vi;

        if(nn > 1)
        {
            /* Row n-1 */
            vi = _mm256_i64gather_pd(&b[nn - 1 + off], vidx, 8);
            vi = _mm256_fnmadd_pd(_mm256_set1_pd(du[nn - 1]), bip2, vi);
            vi = _mm256_div_pd(vi, _mm256_set1_pd(d__[nn - 1]));
            FLA_DGTSV_STORE4(&b[nn - 1 + off], vi);
            bip1 = vi;
        }

        for(i__ = nn - 2; i__ >= 1; --i__)
        {
            vi = _mm256_i64gather_pd(&b[i__ + off], vidx, 8);
            vi = _mm256_fnmadd_pd(_mm256_set1_pd(du[i__]), bip1, vi);
            vi = _mm256_fnmadd_pd(_mm256_set1_pd(dl[i__]), bip2, vi);
            vi = _mm256_div_pd(vi, _mm256_set1_pd(d__[i__]));
            FLA_DGTSV_STORE4(&b[i__ + off], vi);
            bip2 = bip1;
            bip1 = vi;
        }
    }

    /* Remaining columns (nrhs not a multiple of 4) with scalar back substitution. */
    for(; j < nr; ++j)
    {
        aocl_int64_t offset = (j + 1) * b_dim1;
        b[nn + offset] /= d__[nn];
        if(nn > 1)
        {
            b[nn - 1 + offset] = (b[nn - 1 + offset] - du[nn - 1] * b[nn + offset]) / d__[nn - 1];
        }
        for(i__ = nn - 2; i__ >= 1; --i__)
        {
            b[i__ + offset]
                = (b[i__ + offset] - du[i__] * b[i__ + 1 + offset] - dl[i__] * b[i__ + 2 + offset])
                  / d__[i__];
        }
    }

    return 0;
}
#endif
