/*
 *    Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

/* ./dlarf1f.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */

#ifdef FLA_ENABLE_AMD_OPT
#include "fla_lapack_x86_common.h"

void fla_dlarf_left_tuning_params(aocl_int64_t m, aocl_int64_t n, FLA_Bool *use_blocked_flag,
                                  aocl_int64_t *nthreads);
void fla_dlarf_right_tuning_params(aocl_int64_t m, aocl_int64_t n, aocl_int64_t *block_size,
                                   aocl_int64_t *nthreads);

#endif

#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#else
extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
#endif

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static aocl_int64_t c__1 = 1;
/* > \brief \b DLARF1F applies an elementary reflector to a general rectangular */
/* matrix assuming v(1) = 1. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARF1F + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf1f
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf1f
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf1f
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE */
/* INTEGER INCV, LDC, M, N */
/* DOUBLE PRECISION TAU */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION C( LDC, * ), V( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARF1F applies a real elementary reflector H to a real m by n matrix */
/* > C, from either the left or the right. H is represented in the form */
/* > */
/* > H = I - tau * v * v**T */
/* > */
/* > where tau is a real scalar and v is a real vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': form H * C */
/* > = 'R': form C * H */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is DOUBLE PRECISION array, dimension */
/* > (1 + (M-1)*abs(INCV)) if SIDE = 'L' */
/* > or (1 + (N-1)*abs(INCV)) if SIDE = 'R' */
/* > The vector v in the representation of H. V is not used if */
/* > TAU = 0. V(1) is not referenced or modified. */
/* > \endverbatim */
/* > */
/* > \param[in] INCV */
/* > \verbatim */
/* > INCV is INTEGER */
/* > The increment between elements of v. INCV <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION */
/* > The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (LDC,N) */
/* > On entry, the m by n matrix C. */
/* > On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
/* > or C * H if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension */
/* > (N) if SIDE = 'L' */
/* > or (M) if SIDE = 'R' */
/* > \endverbatim */
/* To take advantage of the fact that v(1) = 1, we do the following */
/* v = [ 1 v_2 ]**T */
/* If SIDE='L' */
/* |-----| */
/* | C_1 | */
/* C =| C_2 | */
/* |-----| */
/* C_1\in\mathbb{
R}
^{
1\times n}
, C_2\in\mathbb{
R}
^{
m-1\times n}
*/
/* So we compute: */
/* C = HC = (I - \tau vv**T)C */
/* = C - \tau vv**T C */
/* w = C**T v = [ C_1**T C_2**T ] [ 1 v_2 ]**T */
/* = C_1**T + C_2**T v ( DGEMM then DAXPY ) */
/* C = C - \tau vv**T C */
/* = C - \tau vw**T */
/* Giving us C_1 = C_1 - \tau w**T ( DAXPY ) */
/* and */
/* C_2 = C_2 - \tau v_2w**T ( DGER ) */
/* If SIDE='R' */
/* C = [ C_1 C_2 ] */
/* C_1\in\mathbb{
R}
^{
m\times 1}
, C_2\in\mathbb{
R}
^{
m\times n-1}
*/
/* So we compute: */
/* C = CH = C(I - \tau vv**T) */
/* = C - \tau Cvv**T */
/* w = Cv = [ C_1 C_2 ] [ 1 v_2 ]**T */
/* = C_1 + C_2v_2 ( DGEMM then DAXPY ) */
/* C = C - \tau Cvv**T */
/* = C - \tau wv**T */
/* Giving us C_1 = C_1 - \tau w ( DAXPY ) */
/* and */
/* C_2 = C_2 - \tau wv_2**T ( DGER ) */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup larf */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void dlarf1f_(char *side, aocl_int_t *m, aocl_int_t *n, doublereal *v, aocl_int_t *incv,
              doublereal *tau, doublereal *c__, aocl_int_t *ldc, doublereal *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dlarf1f(side, m, n, v, incv, tau, c__, ldc, work);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t incv_64 = *incv;
    aocl_int64_t ldc_64 = *ldc;

    aocl_lapack_dlarf1f(side, &m_64, &n_64, v, &incv_64, tau, c__, &ldc_64, work);

    *ldc = (aocl_int_t)ldc_64;
#endif
}

void aocl_lapack_dlarf1f(char *side, aocl_int64_t *m, aocl_int64_t *n, doublereal *v,
                         aocl_int64_t *incv, doublereal *tau, doublereal *c__, aocl_int64_t *ldc,
                         doublereal *work)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlarf1f inputs: side %c, m %" FLA_IS ", n %" FLA_IS ", incv %" FLA_IS
                      ", ldc %" FLA_IS "",
                      *side, *m, *n, *incv, *ldc);

    /* System generated locals */
    aocl_int64_t c_dim1, c_offset;
#ifndef FLA_ENABLE_AMD_OPT
    aocl_int64_t i__1;
#endif
    doublereal d__1;
    /* Local variables */
    aocl_int64_t i__;
    logical applyleft;
    aocl_int64_t lastc;
    aocl_int64_t lastv;
    aocl_int64_t istart;
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    applyleft = lsame_(side, "L", 1, 1);
    lastv = 1;
    lastc = 0;
    if(*tau != 0.)
    {
        /* Set up variables for scanning V. LASTV begins pointing to the end */
        /* of V. */
        if(applyleft)
        {
            lastv = *m;
        }
        else
        {
            lastv = *n;
        }
        if(*incv > 0)
        {
            i__ = (lastv - 1) * *incv + 1;
        }
        else
        {
            i__ = 1;
        }
        /* Look for the last non-zero row in V. */
        /* Since we are assuming that V(1) = 1, and it is not stored, so we */
        /* shouldn't access it. */
        while(lastv > 1 && v[i__] == 0.)
        {
            --lastv;
            i__ -= *incv;
        }

        /* If incv < 0, zeros at the end of v appear at the start address.
           V is adjusted to point to the last non-zero element. */
        istart = *incv < 0 ? i__ : *incv + 1;

        if(applyleft)
        {
            /* Scan for the last non-zero column in C(1:lastv,:). */
            lastc = aocl_lapack_iladlc(&lastv, n, &c__[c_offset], ldc);
        }
        else
        {
            /* Scan for the last non-zero row in C(:,1:lastv). */
            lastc = aocl_lapack_iladlr(m, &lastv, &c__[c_offset], ldc);
        }
    }
    if(lastc == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(applyleft)
    {
        /* Form H * C */
        /* Check if lastv = 1. This means v = 1, So we just need to comp */
        /* C := HC = (1-\tau)C. */
        if(lastv == 1)
        {
            /* C(1,1:lastc) := ( 1 - tau ) * C(1,1:lastc) */
            d__1 = 1. - *tau;
#ifndef FLA_ENABLE_AMD_OPT
            aocl_blas_dscal(&lastc, &d__1, &c__[c_offset], ldc);
#else
            fla_dscal(&lastc, &d__1, &c__[c_offset], ldc);
#endif
        }
        else
        {
#ifndef FLA_ENABLE_AMD_OPT
            /* w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1) */
            /* w(1:lastc,1) := C(2:lastv,1:lastc)**T * v(2:lastv,1) */
            i__1 = lastv - 1;
            aocl_blas_dgemv("Transpose", &i__1, &lastc, &c_b4, &c__[c_dim1 + 2], ldc, &v[istart],
                            incv, &c_b5, &work[1], &c__1);
            /* w(1:lastc,1) += C(1,1:lastc)**T * v(1,1) = C(1,1:lastc)**T */
            aocl_blas_daxpy(&lastc, &c_b4, &c__[c_offset], ldc, &work[1], &c__1);
            /* C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**T */
            /* C(1, 1:lastc) := C(...) - tau * v(1,1) * w(1:lastc,1)**T */
            /* = C(...) - tau * w(1:lastc,1)**T */
            d__1 = -(*tau);
            aocl_blas_daxpy(&lastc, &d__1, &work[1], &c__1, &c__[c_offset], ldc);
            /* C(2:lastv,1:lastc) := C(...) - tau * v(2:lastv,1)*w(1:last */
            i__1 = lastv - 1;
            d__1 = -(*tau);
            aocl_blas_dger(&i__1, &lastc, &d__1, &v[istart], incv, &work[1], &c__1,
                           &c__[c_dim1 + 2], ldc);
#else
            FLA_Bool use_blocked = 0;
            aocl_int64_t opt_nthreads = 1;

            /* Initialize global context data */
            aocl_fla_init();
            d__1 = -(*tau);

            fla_dlarf_left_tuning_params(lastv, lastc, &use_blocked, &opt_nthreads);

            /* The fused SIMD kernel is single threaded, so it is only used when no
               threading is selected, and against the blocked path, which is itself
               fused per column, only while C is small enough to stay inside L3 */
            FLA_Bool use_fused_simd
                = (opt_nthreads == 1)
                  && (!use_blocked || lastv * lastc <= FLA_DLARF1_L_FUSED_SIMD_THRESH);

            /* If the size of the matrix is small and incv =1, use the optimized path */
            if(use_fused_simd && *incv == c__1 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
            {
                /* Call optimized routine */
                fla_dlarf1f_small_incv1_simd(lastv, lastc, &c__[c_offset], *ldc, &v[1], d__1,
                                             &work[1]);
            }
            else
            {

                const aocl_int64_t lastv_eff = lastv - 1;

                /* If use_blocked is 1, process in blocks */
                if(use_blocked)
                {
                    /* Process in blocks */
#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp parallel for num_threads(opt_nthreads) private(i__)
#endif
                    /* Loop for each column of C */
                    for(i__ = 1; i__ <= lastc; ++i__)
                    {
                        /* W(i) =  C(1:lastv,i)**T * v(1:lastv,1)  */
                        work[i__] = c__[i__ * *ldc + 1];
                        work[i__] += aocl_blas_ddot(&lastv_eff, &v[istart], incv,
                                                    &c__[i__ * *ldc + 2], &c__1);
                        /* C(1:lastv,i) = C(1:lastv,i) - v(1:lastv,1) * -tau * W(i) */
                        doublereal d__2 = d__1 * work[i__];
                        c__[i__ * *ldc + 1] += d__2;
                        aocl_blas_daxpy(&lastv_eff, &d__2, &v[istart], incv, &c__[i__ * *ldc + 2],
                                        &c__1);
                    }
                }
                else
                {
                    /* Process in a single call */
                    /* w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1) */
#if FLA_ENABLE_AOCL_BLAS && defined(BLIS_KERNELS_ZEN4)
                    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512) && *incv > 0)
                    {
                        /* Use direct single threaded BLIS kernel */
                        bli_dgemv_t_zen4_int(BLIS_CONJUGATE, BLIS_NO_CONJUGATE, lastv_eff, lastc,
                                             &c_b4, &c__[c_offset + 1], 1, *ldc, &v[istart], *incv,
                                             &c_b5, &work[1], c__1, NULL);
                    }
                    else
                    {
#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp teams num_teams(1) thread_limit(1)
#endif
                        {
                            aocl_blas_dgemv("Transpose", &lastv_eff, &lastc, &c_b4,
                                            &c__[c_offset + 1], ldc, &v[istart], incv, &c_b5,
                                            &work[1], &c__1);
                        }
                    }
#else
                    aocl_blas_dgemv("Transpose", &lastv_eff, &lastc, &c_b4, &c__[c_offset + 1], ldc,
                                    &v[istart], incv, &c_b5, &work[1], &c__1);
#endif
                    /* v(1) = 1, so row 1 of C contributes directly to w and is  */
                    /* updated directly from w. Both passes over the (strided)   */
                    /* first row of C are fused into a single loop.              */
                    /* w(1:lastc,1) += C(1,1:lastc)**T * v(1,1) = C(1,1:lastc)**T */
                    /* C(1,1:lastc)  := C(...) - tau * v(1,1) * w(1:lastc,1)**T   */
                    /*               =  C(...) - tau * w(1:lastc,1)**T            */
                    {
                        doublereal *crow = &c__[c_offset];
                        for(i__ = 0; i__ < lastc; ++i__)
                        {
                            doublereal wj = work[i__ + 1] + crow[i__ * *ldc];
                            work[i__ + 1] = wj;
                            crow[i__ * *ldc] += d__1 * wj;
                        }
                    }
                    /* C(2:lastv,1:lastc) := C(...) - tau * v(2:lastv,1)*w(1:last */
                    aocl_blas_dger(&lastv_eff, &lastc, &d__1, &v[istart], incv, &work[1], &c__1,
                                   &c__[c_offset + 1], ldc);
                }
            }
#endif
        }
    }
    else
    {
        /* Form C * H */
        /* Check if n = 1. This means v = 1, so we just need to compute */
        /* C := CH = C(1-\tau). */
        if(lastv == 1)
        {
            /* C(1:lastc,1) := ( 1 - tau ) * C(1:lastc,1) */
            d__1 = 1. - *tau;
#ifndef FLA_ENABLE_AMD_OPT
            aocl_blas_dscal(&lastc, &d__1, &c__[c_offset], &c__1);
#else
            fla_dscal(&lastc, &d__1, &c__[c_offset], &c__1);
#endif
        }
        else
        {
#ifndef FLA_ENABLE_AMD_OPT
            /* w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */
            /* w(1:lastc,1) := C(1:lastc,2:lastv) * v(2:lastv,1) */
            i__1 = lastv - 1;
            aocl_blas_dgemv("No transpose", &lastc, &i__1, &c_b4, &c__[(c_dim1 << 1) + 1], ldc,
                            &v[istart], incv, &c_b5, &work[1], &c__1);
            /* w(1:lastc,1) += C(1:lastc,1) v(1,1) = C(1:lastc,1) */
            aocl_blas_daxpy(&lastc, &c_b4, &c__[c_offset], &c__1, &work[1], &c__1);
            /* C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T */
            /* C(1:lastc,1) := C(...) - tau * w(1:lastc,1) * v(1,1)** */
            /* = C(...) - tau * w(1:lastc,1) */
            d__1 = -(*tau);
            aocl_blas_daxpy(&lastc, &d__1, &work[1], &c__1, &c__[c_offset], &c__1);
            /* C(1:lastc,2:lastv) := C(...) - tau * w(1:lastc,1) * v(2:la */
            i__1 = lastv - 1;
            d__1 = -(*tau);
            aocl_blas_dger(&lastc, &i__1, &d__1, &work[1], &c__1, &v[istart], incv,
                           &c__[(c_dim1 << 1) + 1], ldc);
#else
            aocl_int64_t opt_nthreads = 1;
            aocl_int64_t nb = 0;

            /* Initialize global context data */
            aocl_fla_init();
            d__1 = -(*tau);

            fla_dlarf_right_tuning_params(lastc, lastv, &nb, &opt_nthreads);

            const aocl_int64_t lastv_eff = lastv - 1;

            /* For small sizes the fused SIMD kernel is faster than GEMV + GER,
               whose call overhead dominates at these sizes */
            if(opt_nthreads == 1 && *incv == c__1 && (*ldc & FLA_DLARF1_R_SIMD_LDC_ALIAS_MASK) != 0
               && lastv_eff <= FLA_DLARF1_R_SIMD_COLS_THRESH
               && lastc * lastv_eff <= FLA_DLARF1_R_SIMD_THRESH
               && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
            {
                fla_dlarf1_small_right_incv1_simd(lastc, lastv_eff, &c__[(c_dim1 << 1) + 1], *ldc,
                                                  &v[istart], &c__[c_offset], d__1, &work[1]);
            }
            /* If nb is non zero, process in blocks */
            else if(nb && opt_nthreads > 1)
            {
                /* The first panel will process starting unaligned elements
                 * to ensure that all other panels aligned memory addresses
                 */
                uint64_t unaligned_bytes
                    = ((uint64_t)(c__ + c_offset)) % ((nb * sizeof(doublereal)));
                aocl_int64_t first_thread_rows
                    = fla_min((nb - (unaligned_bytes / sizeof(doublereal))), lastc);
                aocl_int64_t panels
                    = (!!first_thread_rows) + (((lastc - first_thread_rows) + (nb - 1)) / nb);

#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp parallel for num_threads(opt_nthreads) private(i__)
#endif
                for(i__ = 1; i__ <= panels; i__ += 1)
                {
                    aocl_int64_t completed_rows = i__ == 1 ? 0 : (i__ - 2) * nb + first_thread_rows;
                    aocl_int64_t current_block_size
                        = i__ == 1 ? first_thread_rows : fla_min(nb, lastc - completed_rows);
                    aocl_int64_t cur_idx = completed_rows + 1;

                    /* w(1:lastc,1) := C(1:lastc,2:lastv) * v(2:lastv,1) */
                    aocl_blas_dgemv("No transpose", &current_block_size, &lastv_eff, &c_b4,
                                    &c__[(c_dim1 << 1) + cur_idx], ldc, &v[istart], incv, &c_b5,
                                    &work[cur_idx], &c__1);
                    /* v(1) = 1, so column 1 of C contributes directly to w and is */
                    /* updated directly from w. Both passes over the first column  */
                    /* of C are fused into a single kernel call.                   */
                    fla_dlarf1_right_update_c1_simd(current_block_size, &c__[c_dim1 + cur_idx],
                                                    d__1, &work[cur_idx]);
                    /* C(1:lastc,2:lastv) := C(...) - tau * w(1:lastc,1) * v(2:lastv,1)**T */
                    aocl_blas_dger(&current_block_size, &lastv_eff, &d__1, &work[cur_idx], &c__1,
                                   &v[istart], incv, &c__[(c_dim1 << 1) + cur_idx], ldc);
                }
            }
            else
            {
                /* w(1:lastc,1) := C(1:lastc,2:lastv) * v(2:lastv,1) */
#if FLA_ENABLE_AOCL_BLAS && defined(BLIS_KERNELS_ZEN4)
                if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512) && *incv > 0)
                {
                    bli_dgemv_n_zen4_int_40x2_st(BLIS_NO_TRANSPOSE, BLIS_NO_CONJUGATE, lastc,
                                                 lastv_eff, &c_b4, &c__[(c_dim1 << 1) + 1], c__1,
                                                 c_dim1, &v[istart], *incv, &c_b5, &work[1], c__1,
                                                 NULL);
                }
                else
                {
                    aocl_blas_dgemv("No transpose", &lastc, &lastv_eff, &c_b4,
                                    &c__[(c_dim1 << 1) + 1], ldc, &v[istart], incv, &c_b5, &work[1],
                                    &c__1);
                }
#else
                aocl_blas_dgemv("No transpose", &lastc, &lastv_eff, &c_b4, &c__[(c_dim1 << 1) + 1],
                                ldc, &v[istart], incv, &c_b5, &work[1], &c__1);
#endif
                /* v(1) = 1, so column 1 of C contributes directly to w and is */
                /* updated directly from w. Both passes over the first column  */
                /* of C are fused into a single kernel call.                   */
                /* w(1:lastc,1) += C(1:lastc,1) * v(1,1) = C(1:lastc,1)        */
                /* C(1:lastc,1)  = C(1:lastc,1) - tau * w(1:lastc,1) * v(1,1)  */
                fla_dlarf1_right_update_c1_simd(lastc, &c__[c_offset], d__1, &work[1]);

                /* C(1:lastc,2:lastv) := C(...) - w(1:lastc,1) * v(2:lastv,1)**T */
                aocl_blas_dger(&lastc, &lastv_eff, &d__1, &work[1], &c__1, &v[istart], incv,
                               &c__[(c_dim1 << 1) + 1], ldc);
            }
#endif
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLARF1F */
}
/* dlarf1f_ */
