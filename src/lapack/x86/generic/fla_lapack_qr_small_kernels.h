/******************************************************************************
 ** Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 ******************************************************************************/

/*! @file fla_lapack_qr_small_kernels.h
 *  @brief Common front-end functions of QR factorization for small sizes
 *         for double precision to choose optimized paths.
 */

#ifndef FLA_LAPACK_QR_SMALL_KERNELS_H
#define FLA_LAPACK_QR_SMALL_KERNELS_H

#if FLA_ENABLE_AMD_OPT

#include "FLAME.h"

/*
 * DORGQR argument check (same rules as dorgqr_check.c / lapack_dorgqr): sets work[0], validates
 * arguments, xerbla on error.
 *
 * IMPORTANT: This macro uses `return;` to exit early (error / lquery / quick return).
 * It MUST be expanded only inside a `void` function. Do NOT use from non-void functions
 */
#define LAPACK_DORGQR_CHECK_INLINED(m, n, k, lda, work, lwork, info) \
    {                                                                \
        aocl_int64_t _dqc_lquery = ((lwork) == -1);                  \
        *(info) = 0;                                                 \
        if((m) < 0)                                                  \
        {                                                            \
            *(info) = -1;                                            \
        }                                                            \
        else if((n) < 0 || (n) > (m))                                \
        {                                                            \
            *(info) = -2;                                            \
        }                                                            \
        else if((k) < 0 || (k) > (n))                                \
        {                                                            \
            *(info) = -3;                                            \
        }                                                            \
        else if((lda) < fla_max(1, (m)))                             \
        {                                                            \
            *(info) = -5;                                            \
        }                                                            \
        else if((lwork) < fla_max(1, (n)) && !_dqc_lquery)           \
        {                                                            \
            *(info) = -8;                                            \
        }                                                            \
        if(*(info) != 0)                                             \
        {                                                            \
            aocl_int64_t _dqc_i1 = -(*(info));                       \
            aocl_blas_xerbla("DORGQR", &_dqc_i1, (ftnlen)6);         \
            AOCL_DTL_TRACE_LOG_EXIT;                                 \
            return;                                                  \
        }                                                            \
        else if(_dqc_lquery)                                         \
        {                                                            \
            LAPACK_DORGQR_LWKOPT_INLINED(m, n, k, work);             \
            AOCL_DTL_TRACE_LOG_EXIT;                                 \
            return;                                                  \
        }                                                            \
        else if((n) <= 0)                                            \
        {                                                            \
            (work)[0] = 1.;                                          \
            AOCL_DTL_TRACE_LOG_EXIT;                                 \
            return;                                                  \
        }                                                            \
        LAPACK_DORGQR_LWKOPT_INLINED(m, n, k, work);                 \
    }

/* To get lwork value for work buffer */
#define LAPACK_DORGQR_LWKOPT_INLINED(m, n, k, work)                 \
    {                                                               \
        if((n) == 1)                                                \
        {                                                           \
            (work)[0] = 32.;                                        \
        }                                                           \
        else if((n) <= 6)                                           \
        {                                                           \
            (work)[0] = 192.;                                       \
        }                                                           \
        else                                                        \
        {                                                           \
            /* Small-kernel dispatch skips ilaenv on every call. */ \
            (work)[0] = ((n)*32);                                   \
        }                                                           \
    }

/*
 * LAPACK_ORGQR_SMALL_COMMON
 * -------------------------
 * Common code for the small DORGQR kernel.
 *
 * Why the `type` parameter:
 * - Allows the same kernel structure to be reused for other real precisions by
 *   passing `real` or `doublereal` for temporaries and pointers, while preserving the
 *   same loop structure and dispatch logic.
 */
#define LAPACK_ORGQR_SMALL_COMMON(m, n, k, a, lda, tau, work, type)                              \
    {                                                                                            \
        aocl_int64_t _ii, _jj;                                                                   \
        aocl_int64_t _base_j;                                                                    \
        /* Columns k..n-1 exist only when k < n. */                                              \
        /* Init columns k..n-1 to unit vectors e_j (zero column, then 1 on diagonal). */         \
        for(_jj = (k); _jj < (n); ++_jj)                                                         \
        {                                                                                        \
            _base_j = (lda)*_jj;                                                                 \
            for(_ii = 0; _ii < (m); ++_ii)                                                       \
                (a)[_ii + _base_j] = 0.0;                                                        \
            (a)[_jj + _base_j] = 1.0;                                                            \
        }                                                                                        \
        /* Apply reflectors in reverse order j=k-1..0 if k>0. else identity Q already formed. */ \
        if((k) > 0)                                                                              \
        {                                                                                        \
            aocl_int64_t _cc, _base_c, _len_v;                                                   \
            type _tau_j, _neg_tau, _s, _dot;                                                     \
            type *_v, *_col;                                                                     \
            /* Reflectors _jj = k-1 .. 0. Tall m>n uses a two-pass update. */                    \
            if((m) > (n))                                                                        \
            {                                                                                    \
                aocl_int64_t _nt, _trail_base, _row_lin_base;                                    \
                type _wk0_save;                                                                  \
                /* Preserve WORK(1); this branch reuses work[0..n-j-2] as temporary dot storage. \
                 */                                                                              \
                _wk0_save = (work)[0];                                                           \
                for(_jj = (k)-1; _jj >= 0; --_jj)                                                \
                {                                                                                \
                    _tau_j = (tau)[_jj];                                                         \
                    _base_j = (lda)*_jj;                                                         \
                    if(_tau_j == 0.0)                                                            \
                    {                                                                            \
                        /* H = I: column j of Q is e_j (DORG2R / DLARF when TAU = 0). */         \
                        for(_ii = 0; _ii < (m); ++_ii)                                           \
                            (a)[_ii + _base_j] = 0.0;                                            \
                        (a)[_jj + _base_j] = 1.0;                                                \
                    }                                                                            \
                    else                                                                         \
                    {                                                                            \
                        _v = &(a)[_jj + _base_j]; /* reflector col _jj; v[0] at row _jj */       \
                        _v[0] = 1.0;                                                             \
                        _len_v = (m) - (_jj); /* Householder length: rows _jj .. m-1 */          \
                        /* Tall m>n: two-pass; work[c] = tau * v'*col_c. */                      \
                        _nt = (n) - (_jj)-1;                                                     \
                        /* Base offset into column-major A for trailing cols j+1..n-1. */        \
                        _trail_base = (lda) * ((_jj) + 1);                                       \
                        for(_cc = 0; _cc < _nt; ++_cc)                                           \
                        {                                                                        \
                            _base_c = _trail_base + (lda) * (_cc);                               \
                            _dot = 0.0;                                                          \
                            _col = &(a)[_jj + _base_c];                                          \
                            for(_ii = 0; _ii < _len_v; ++_ii)                                    \
                                _dot += _v[_ii] * _col[_ii];                                     \
                            (work)[_cc] = _tau_j * _dot;                                         \
                        }                                                                        \
                        for(_ii = 0; _ii < _len_v; ++_ii)                                        \
                        {                                                                        \
                            _s = _v[_ii];                                                        \
                            /* Linear index of A(_jj+_ii, _jj+1); step by lda per col. */        \
                            _row_lin_base = _jj + _ii + _trail_base;                             \
                            for(_cc = 0; _cc < _nt; ++_cc)                                       \
                                (a)[_row_lin_base + (lda) * (_cc)] -= (work)[_cc] * _s;          \
                        }                                                                        \
                        _neg_tau = -_tau_j;                                                      \
                        /* Column _jj of Q: scale v below diag; v(_jj,_jj)=1-tau. */             \
                        for(_ii = 1; _ii < _len_v; ++_ii)                                        \
                            _v[_ii] *= _neg_tau;                                                 \
                        _v[0] = 1.0 - _tau_j;                                                    \
                        for(_ii = 0; _ii < _jj; ++_ii)                                           \
                            (a)[_ii + _base_j] = 0.0;                                            \
                    }                                                                            \
                }                                                                                \
                (work)[0] = _wk0_save;                                                           \
            }                                                                                    \
            else                                                                                 \
            {                                                                                    \
                /* Square case: apply one trailing column at a time, matching DORG2R order. */   \
                for(_jj = (k)-1; _jj >= 0; --_jj)                                                \
                {                                                                                \
                    _tau_j = (tau)[_jj];                                                         \
                    _base_j = (lda)*_jj;                                                         \
                    if(_tau_j == 0.0)                                                            \
                    {                                                                            \
                        /* H = I: column j of Q is e_j (DORG2R / DLARF when TAU = 0). */         \
                        for(_ii = 0; _ii < (m); ++_ii)                                           \
                            (a)[_ii + _base_j] = 0.0;                                            \
                        (a)[_jj + _base_j] = 1.0;                                                \
                    }                                                                            \
                    else                                                                         \
                    {                                                                            \
                        _v = &(a)[_jj + _base_j]; /* reflector col _jj; v[0] at row _jj */       \
                        _v[0] = 1.0;                                                             \
                        _len_v = (m) - (_jj); /* Householder length: rows _jj .. m-1 */          \
                        /* Apply one trailing column at a time (DORG2R order). */                \
                        for(_cc = _jj + 1; _cc < (n); ++_cc)                                     \
                        {                                                                        \
                            _base_c = (lda)*_cc;                                                 \
                            _col = &(a)[_jj + _base_c];                                          \
                            _dot = 0.0;                                                          \
                            for(_ii = 0; _ii < _len_v; ++_ii)                                    \
                                _dot += _v[_ii] * _col[_ii];                                     \
                            _s = _tau_j * _dot;                                                  \
                            for(_ii = 0; _ii < _len_v; ++_ii)                                    \
                                _col[_ii] -= _s * _v[_ii];                                       \
                        }                                                                        \
                        _neg_tau = -_tau_j;                                                      \
                        /* Column _jj of Q: scale v below diag; v(_jj,_jj)=1-tau. */             \
                        for(_ii = 1; _ii < _len_v; ++_ii)                                        \
                            _v[_ii] *= _neg_tau;                                                 \
                        _v[0] = 1.0 - _tau_j;                                                    \
                        for(_ii = 0; _ii < _jj; ++_ii)                                           \
                            (a)[_ii + _base_j] = 0.0;                                            \
                    }                                                                            \
                }                                                                                \
            }                                                                                    \
        }                                                                                        \
    }

/*
 * LAPACK_DORGQR_SMALL
 * ---------------------
 * Generic “small” DORGQR kernel with a scalar type parameter (`type`).
 *
 * Why the `type` parameter:
 * - Allows the same kernel structure to be reused for other real precisions by
 *   passing `real` or `doublereal` for temporaries and pointers, while preserving the
 *   same loop structure and dispatch logic.
 *
 * Key behaviors:
 * - Initializes columns k..n-1 to unit vectors e_j (DORG2R “columns k+1:n” in 1-based).
 *   For k==0 this is the entire Q block; there are no reflectors to apply (same as reference).
 * - Applies reflectors in reverse order j=k-1..0 when k>0.
 * - Tall (`m > n`) uses a two-pass update for trailing columns and preserves `work[0]`
 *   while using WORK as temporary dot storage.
 * - Square (`m == n`) uses a DORG2R-style column-wise update.
 * - Restores `work[0]` after use if it is repurposed as scratch.
 */
#define LAPACK_DORGQR_SMALL(m, n, k, a, lda, tau, work, lwork, info)       \
    {                                                                      \
        LAPACK_DORGQR_CHECK_INLINED(m, n, k, lda, work, lwork, info);      \
        LAPACK_ORGQR_SMALL_COMMON(m, n, k, a, lda, tau, work, doublereal); \
    }
#endif // #if FLA_ENABLE_AMD_OPT
#endif // #ifndef FLA_LAPACK_QR_SMALL_KERNELS_H
