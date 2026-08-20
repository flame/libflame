/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

/*
 * AVX512 optimized DGTSV.
 *
 * Solves A * X = B where A is an n-by-n tridiagonal matrix, using Gaussian
 * elimination with partial pivoting (identical algorithm to reference LAPACK
 * DGTSV).
 *
 * The tridiagonal factorization (arrays dl, d, du) carries a strict data
 * dependency from one row to the next, so it is computed with scalar code.
 * The update of the right-hand-side matrix B is independent across the NRHS
 * columns and is vectorized 8 columns at a time.  Because B is stored
 * column-major, a fixed row across consecutive columns is strided by ldb, so
 * AVX512 gather/scatter instructions are used for the load/store.
 */

int fla_dgtsv_kernel_avx512(aocl_int64_t *n, aocl_int64_t *nrhs, doublereal *dl, doublereal *d__,
                            doublereal *du, doublereal *b, aocl_int64_t *ldb, aocl_int64_t *info)
{
    aocl_int64_t b_dim1, b_offset, i__, j, nn, nr;
    doublereal fact, temp;
    __m512i vidx;

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

    /* Column-offset index vector: { 0, ldb, 2*ldb, ..., 7*ldb } (in elements). */
    vidx = _mm512_set_epi64(7 * b_dim1, 6 * b_dim1, 5 * b_dim1, 4 * b_dim1, 3 * b_dim1, 2 * b_dim1,
                            b_dim1, 0);

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
                    __m512d fact_vec = _mm512_set1_pd(fact);
                    doublereal *bi = &b[i__ + b_dim1];
                    doublereal *bi1 = &b[i__ + 1 + b_dim1];
                    for(j = 0; j + 8 <= nr; j += 8)
                    {
                        aocl_int64_t off = j * b_dim1;
                        __m512d vbi = _mm512_i64gather_pd(vidx, bi + off, 8);
                        __m512d vbi1 = _mm512_i64gather_pd(vidx, bi1 + off, 8);
                        vbi1 = _mm512_fnmadd_pd(fact_vec, vbi, vbi1);
                        _mm512_i64scatter_pd(bi1 + off, vidx, vbi1, 8);
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
                __m512d fact_vec = _mm512_set1_pd(fact);
                doublereal *bi = &b[i__ + b_dim1];
                doublereal *bi1 = &b[i__ + 1 + b_dim1];
                for(j = 0; j + 8 <= nr; j += 8)
                {
                    aocl_int64_t off = j * b_dim1;
                    __m512d vbi = _mm512_i64gather_pd(vidx, bi + off, 8);
                    __m512d vbi1 = _mm512_i64gather_pd(vidx, bi1 + off, 8);
                    __m512d newbi1 = _mm512_fnmadd_pd(fact_vec, vbi1, vbi);
                    _mm512_i64scatter_pd(bi + off, vidx, vbi1, 8);
                    _mm512_i64scatter_pd(bi1 + off, vidx, newbi1, 8);
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
    for(j = 0; j + 8 <= nr; j += 8)
    {
        /* Base offset points at column j (0-based) => 1-based column (j+1). */
        aocl_int64_t off = b_dim1 + j * b_dim1;
        __m512d bip1, bip2, vi;

        /* Row n */
        vi = _mm512_i64gather_pd(vidx, &b[nn + off], 8);
        vi = _mm512_div_pd(vi, _mm512_set1_pd(d__[nn]));
        _mm512_i64scatter_pd(&b[nn + off], vidx, vi, 8);
        bip2 = vi;
        bip1 = vi;

        if(nn > 1)
        {
            /* Row n-1 */
            vi = _mm512_i64gather_pd(vidx, &b[nn - 1 + off], 8);
            vi = _mm512_fnmadd_pd(_mm512_set1_pd(du[nn - 1]), bip2, vi);
            vi = _mm512_div_pd(vi, _mm512_set1_pd(d__[nn - 1]));
            _mm512_i64scatter_pd(&b[nn - 1 + off], vidx, vi, 8);
            bip1 = vi;
        }

        for(i__ = nn - 2; i__ >= 1; --i__)
        {
            vi = _mm512_i64gather_pd(vidx, &b[i__ + off], 8);
            vi = _mm512_fnmadd_pd(_mm512_set1_pd(du[i__]), bip1, vi);
            vi = _mm512_fnmadd_pd(_mm512_set1_pd(dl[i__]), bip2, vi);
            vi = _mm512_div_pd(vi, _mm512_set1_pd(d__[i__]));
            _mm512_i64scatter_pd(&b[i__ + off], vidx, vi, 8);
            bip2 = bip1;
            bip1 = vi;
        }
    }

    /* Remaining columns (nrhs not a multiple of 8) with scalar back substitution. */
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
