/*
 *  Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

/* ./zlarf1f.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
#if FLA_ENABLE_AMD_OPT
#include "fla_lapack_x86_common.h"
#endif
static dcomplex c_b1 = {1., 0.};
static dcomplex c_b2 = {0., 0.};
static aocl_int64_t c__1 = 1;
/* > \brief \b ZLARF1F applies an elementary reflector to a general rectangular */
/* matrix assuming v(1) = 1. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLARF1F + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarf1f
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarf1f
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarf1f
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE */
/* INTEGER INCV, LDC, M, N */
/* COMPLEX*16 TAU */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 C( LDC, * ), V( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARF1F applies a complex elementary reflector H to a real m by n matrix */
/* > C, from either the left or the right. H is represented in the form */
/* > */
/* > H = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar and v is a complex vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > */
/* > To apply H**H, supply conjg(tau) instead */
/* > tau. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': form H * C */
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
/* > V is COMPLEX*16 array, dimension */
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
/* > TAU is COMPLEX*16 */
/* > The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX*16 array, dimension (LDC,N) */
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
/* > WORK is COMPLEX*16 array, dimension */
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
C}
^{
1\times n}
, C_2\in\mathbb{
C}
^{
m-1\times n}
*/
/* So we compute: */
/* C = HC = (I - \tau vv**T)C */
/* = C - \tau vv**T C */
/* w = C**T v = [ C_1**T C_2**T ] [ 1 v_2 ]**T */
/* = C_1**T + C_2**T v ( ZGEMM then ZAXPYC-like ) */
/* C = C - \tau vv**T C */
/* = C - \tau vw**T */
/* Giving us C_1 = C_1 - \tau w**T ( ZAXPYC-like ) */
/* and */
/* C_2 = C_2 - \tau v_2w**T ( ZGERC ) */
/* If SIDE='R' */
/* C = [ C_1 C_2 ] */
/* C_1\in\mathbb{
C}
^{
m\times 1}
, C_2\in\mathbb{
C}
^{
m\times n-1}
*/
/* So we compute: */
/* C = CH = C(I - \tau vv**T) */
/* = C - \tau Cvv**T */
/* w = Cv = [ C_1 C_2 ] [ 1 v_2 ]**T */
/* = C_1 + C_2v_2 ( ZGEMM then ZAXPYC-like ) */
/* C = C - \tau Cvv**T */
/* = C - \tau wv**T */
/* Giving us C_1 = C_1 - \tau w ( ZAXPYC-like ) */
/* and */
/* C_2 = C_2 - \tau wv_2**T ( ZGERC ) */
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
void zlarf1f_(char *side, aocl_int_t *m, aocl_int_t *n, dcomplex *v, aocl_int_t *incv,
              dcomplex *tau, dcomplex *c__, aocl_int_t *ldc, dcomplex *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlarf1f(side, m, n, v, incv, tau, c__, ldc, work);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t incv_64 = *incv;
    aocl_int64_t ldc_64 = *ldc;

    aocl_lapack_zlarf1f(side, &m_64, &n_64, v, &incv_64, tau, c__, &ldc_64, work);

    *ldc = (aocl_int_t)ldc_64;
#endif
}

void aocl_lapack_zlarf1f(char *side, aocl_int64_t *m, aocl_int64_t *n, dcomplex *v,
                         aocl_int64_t *incv, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc,
                         dcomplex *work)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlarf1f inputs: side %c, m %" FLA_IS ", n %" FLA_IS ", incv %" FLA_IS
                      ", ldc %" FLA_IS "",
                      *side, *m, *n, *incv, *ldc);

    /* System generated locals */
    aocl_int64_t c_dim1, c_offset, i__1, i__2, i__3;
    dcomplex z__1, z__2, z__3;
    /* Builtin functions */
    void d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t i__;
    logical applyleft;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
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
    /* .. Intrinsic Functions .. */
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
    if(tau->real != 0. || tau->imag != 0.)
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
        for(;;)
        {
            /* while(complicated condition) */
            i__1 = i__;
            if(!(lastv > 1 && (v[i__1].real == 0. && v[i__1].imag == 0.)))
                break;
            --lastv;
            i__ -= *incv;
        }

        /* If incv < 0, zeros at the end of v appear at the start address.
           V is adjusted to point to the last non-zero element. */
        istart = *incv < 0 ? i__ : *incv + 1;

        if(applyleft)
        {
            /* Scan for the last non-zero column in C(1:lastv,:). */
            lastc = aocl_lapack_ilazlc(&lastv, n, &c__[c_offset], ldc);
        }
        else
        {
            /* Scan for the last non-zero row in C(:,1:lastv). */
            lastc = aocl_lapack_ilazlr(m, &lastv, &c__[c_offset], ldc);
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
        /* Check if m = 1. This means v = 1, So we just need to compu */
        /* C := HC = (1-\tau)C. */
        if(lastv == 1)
        {
            z__1.real = 1. - tau->real;
            z__1.imag = 0. - tau->imag; // , expr subst
            aocl_blas_zscal(&lastc, &z__1, &c__[c_offset], ldc);
        }
        else
        {
#if FLA_ENABLE_AMD_OPT
            aocl_fla_init();
            if(lastv <= FLA_ZLARF1_SMALL_DIM_THRESH && lastc <= FLA_ZLARF1_SMALL_DIM_THRESH
               && *incv == c__1 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
            {
                z__1.real = -tau->real;
                z__1.imag = -tau->imag; // , expr subst
                fla_zlarf1f_small_left_incv1_simd(lastv, lastc, c__, *ldc, v, &z__1, work);
            }
            else
            {
#endif
            /* w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1) */
            /* (I - tvv**H)C = C - tvv**H C */
            /* First compute w**H = v**H c -> w = C**H v */
            /* C = [ C_1 C_2 ]**T, v = [1 v_2]**T */
            /* w = C_1**H + C_2**Hv_2 */
            /* w = C_2**Hv_2 */
            i__1 = lastv - 1;
            aocl_blas_zgemv("Conjugate transpose", &i__1, &lastc, &c_b1, &c__[c_dim1 + 2], ldc,
                            &v[istart], incv, &c_b2, &work[1], &c__1);
            /* w(1:lastc,1) += v(1,1) * C(1,1:lastc)**H */
            i__1 = lastc;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = i__;
                i__3 = i__;
                d_cnjg(&z__2, &c__[i__ * c_dim1 + 1]);
                z__1.real = work[i__3].real + z__2.real;
                z__1.imag = work[i__3].imag + z__2.imag; // , expr subst
                work[i__2].real = z__1.real;
                work[i__2].imag = z__1.imag; // , expr subst
            }
            /* C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**H */
            /* C(1, 1:lastc) := C(...) - tau * v(1,1) * w(1:lastc,1)**H */
            /* = C(...) - tau * Conj(w(1:lastc,1)) */
            /* This is essentially a zaxpyc */
            i__1 = lastc;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = i__ * c_dim1 + 1;
                i__3 = i__ * c_dim1 + 1;
                d_cnjg(&z__3, &work[i__]);
                z__2.real = tau->real * z__3.real - tau->imag * z__3.imag;
                z__2.imag = tau->real * z__3.imag + tau->imag * z__3.real; // , expr subst
                z__1.real = c__[i__3].real - z__2.real;
                z__1.imag = c__[i__3].imag - z__2.imag; // , expr subst
                c__[i__2].real = z__1.real;
                c__[i__2].imag = z__1.imag; // , expr subst
            }
            /* C(2:lastv,1:lastc) += - tau * v(2:lastv,1) * w(1:lastc,1)**H */
            i__1 = lastv - 1;
            z__1.real = -tau->real;
            z__1.imag = -tau->imag; // , expr subst
            aocl_blas_zgerc(&i__1, &lastc, &z__1, &v[istart], incv, &work[1], &c__1,
                            &c__[c_dim1 + 2], ldc);
#if FLA_ENABLE_AMD_OPT
            }
#endif
        }
    }
    else
    {
        /* Form C * H */
        /* Check if n = 1. This means v = 1, so we just need to compu */
        /* C := CH = C(1-\tau). */
        if(lastv == 1)
        {
            z__1.real = 1. - tau->real;
            z__1.imag = 0. - tau->imag; // , expr subst
            aocl_blas_zscal(&lastc, &z__1, &c__[c_offset], &c__1);
        }
        else
        {
#if FLA_ENABLE_AMD_OPT
            aocl_fla_init();
            if(lastv <= FLA_ZLARF1_SMALL_DIM_THRESH && lastc <= FLA_ZLARF1_SMALL_DIM_THRESH
               && *incv == c__1 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
            {
                aocl_int64_t lastv_eff = lastv - 1;

                z__1.real = -tau->real;
                z__1.imag = -tau->imag; // , expr subst
                fla_zlarf1_small_right_incv1_simd(lastc, lastv_eff, &c__[(c_dim1 << 1) + 1],
                                                  *ldc, &v[istart], &c__[c_offset], &z__1,
                                                  &work[1]);
            }
            else
            {
#endif
            /* w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */
            /* w(1:lastc,1) := C(1:lastc,2:lastv) * v(2:lastv,1) */
            i__1 = lastv - 1;
            aocl_blas_zgemv("No transpose", &lastc, &i__1, &c_b1, &c__[(c_dim1 << 1) + 1], ldc,
                            &v[istart], incv, &c_b2, &work[1], &c__1);
            /* w(1:lastc,1) += C(1:lastc,1) v(1,1) = C(1:lastc,1) */
            aocl_blas_zaxpy(&lastc, &c_b1, &c__[c_offset], &c__1, &work[1], &c__1);
            /* C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T */
            /* C(1:lastc,1) := C(...) - tau * w(1:lastc,1) * v(1,1 */
            /* = C(...) - tau * w(1:lastc,1) */
            z__1.real = -tau->real;
            z__1.imag = -tau->imag; // , expr subst
            aocl_blas_zaxpy(&lastc, &z__1, &work[1], &c__1, &c__[c_offset], &c__1);
            /* C(1:lastc,2:lastv) := C(...) - tau * w(1:lastc,1) * v(2 */
            i__1 = lastv - 1;
            z__1.real = -tau->real;
            z__1.imag = -tau->imag; // , expr subst
            aocl_blas_zgerc(&lastc, &i__1, &z__1, &work[1], &c__1, &v[istart], incv,
                            &c__[(c_dim1 << 1) + 1], ldc);
#if FLA_ENABLE_AMD_OPT
            }
#endif
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLARF1F */
}
/* zlarf1f_ */
