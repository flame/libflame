/* ../netlib/zrot.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */

/*
    Copyright (c) 2022 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLA_f2c.h" /* > \brief \b ZROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors. */
#ifdef FLA_ENABLE_AMD_OPT
#include "immintrin.h"
#endif

/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZROT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zrot.f" > */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zrot.f" > */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zrot.f" > */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, INCY, N */
/* DOUBLE PRECISION C */
/* COMPLEX*16 S */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 CX( * ), CY( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZROT applies a plane rotation, where the cos (C) is real and the */
/* > sin (S) is complex, and the vectors CX and CY are complex. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of elements in the vectors CX and CY. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CX */
/* > \verbatim */
/* > CX is COMPLEX*16 array, dimension (N) */
/* > On input, the vector X. */
/* > On output, CX is overwritten with C*X + S*Y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive values of CY. INCX <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CY */
/* > \verbatim */
/* > CY is COMPLEX*16 array, dimension (N) */
/* > On input, the vector Y. */
/* > On output, CY is overwritten with -CONJG(S)*X + C*Y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > The increment between successive values of CY. INCX <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is COMPLEX*16 */
/* > C and S define a rotation */
/* > [ C S ] */
/* > [ -conjg(S) C ] */
/* > where C*C + S*CONJG(S) = 1.0. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int zrot_(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zrot inputs: n %" FLA_IS ", incx %" FLA_IS ", incy %" FLA_IS "",*n, *incx, *incy);
    extern fla_context global_context;
    extern int fla_zrot_avx2(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s);
    extern int fla_zrot_native(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s);

    /* Initialize global context data */
    aocl_fla_init();

    int retval = 0;
#ifdef FLA_ENABLE_AMD_OPT
    if (global_context.is_avx2)
    {
      retval = fla_zrot_avx2(n, cx, incx, cy, incy, c__, s);
    }
    else
    {
      retval = fla_zrot_native(n, cx, incx, cy, incy, c__, s);
    }
#else
      retval = fla_zrot_native(n, cx, incx, cy, incy, c__, s);
#endif

    AOCL_DTL_TRACE_LOG_EXIT
    return retval;
}

int fla_zrot_avx2(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s)
{
    /* System generated locals */
    integer i__1;
    doublecomplex z__1, z__2, z__3;
    /* Local variables */
    integer i__, ix, iy;
    doublereal lc, sr, si, msi;

    __m256d cmm, srmm, simm, sinm;
    __m256d sirmm, srimm, msirmm, msrimm;
    __m256d xmm0, ymm0, xmm1, ymm1;
    __m256d xrmm0, yrmm0, ximm0, yimm0;
    __m256d xrmm1, yrmm1, ximm1, yimm1;
    __m256d oxm0, oym0, oxm1, oym1;

    __m128d cm, srm, sim, sin;
    __m128d sirm, srim, msirm, msrim;
    __m128d xmm, ymm;
    __m128d xrmm, yrmm, ximm, yimm;
    __m128d oxm, oym;
    __m128d hxmm0, hxmm1, hymm0, hymm1;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --cy;
    --cx;
    /* Function Body */

    if (*n <= 0)
    {
        return 0;
    }
    lc  = *c__;
    sr  = s->r;
    si  = s->i;
    msi = -si;

    cmm  = _mm256_broadcast_sd((double const *) &lc);
    srmm = _mm256_broadcast_sd((double const *) &sr);
    simm = _mm256_broadcast_sd((double const *) &si);
    sinm = _mm256_broadcast_sd((double const *) &msi);

    sirmm = _mm256_shuffle_pd(srmm, simm, 0xA);  
    srimm = _mm256_shuffle_pd(simm, srmm, 0x5);  
    msirmm = _mm256_shuffle_pd(srmm, sinm, 0xA); 
    msrimm = _mm256_shuffle_pd(sinm, srmm, 0x5); 

    cm  = _mm_loaddup_pd ((double const *) &lc);
    srm = _mm_loaddup_pd ((double const *) &sr);
    sim = _mm_loaddup_pd ((double const *) &si);
    sin = _mm_loaddup_pd ((double const *) &msi);

    sirm = _mm_shuffle_pd(srm, sim, 0x2); 
    srim = _mm_shuffle_pd(sim, srm, 0x1);
    msirm = _mm_shuffle_pd(srm, sin, 0x2);
    msrim = _mm_shuffle_pd(sin, srm, 0x1);


    if (*incx == 1 && *incy == 1)
    {
        goto L20;
    }
    /* Code for unequal increments or equal increments not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }

    i__1 = *n;
    if (*incx != *incy)
    {
        for (i__ = 1; i__ <= i__1; ++i__)
        {
            z__2.r = lc * cx[ix].r;
            z__2.i = lc * cx[ix].i; // , expr subst
            z__3.r = sr * cy[iy].r - si * cy[iy].i;
            z__3.i = sr * cy[iy].i + si * cy[iy].r; // , expr subst
            z__1.r = z__2.r + z__3.r;
            z__1.i = z__2.i + z__3.i; // , expr subst

            z__2.r = lc * cy[iy].r;
            z__2.i = lc * cy[iy].i; // , expr subst
            z__3.r = sr * cx[ix].r + si * cx[ix].i;
            z__3.i = sr * cx[ix].i - si * cx[ix].r; // , expr subst

            cy[iy].r = z__2.r - z__3.r;
            cy[iy].i = z__2.i - z__3.i; // , expr subst
            cx[ix].r = z__1.r;
            cx[ix].i = z__1.i; // , expr subst
            ix += *incx;
            iy += *incy;
        }
    }
    else
    {
        for (i__ = 1; i__ <= (i__1 - 1); i__ += 2)
        {
            /* load complex inputs from x & y */
            xmm0   = _mm256_loadu_pd((double const *) &cx[ix]);
            hxmm1 = _mm_loadu_pd((double const *) &cx[ix + *incx]);
            ymm0   = _mm256_loadu_pd((double const *) &cy[ix]);
            hymm1 = _mm_loadu_pd((double const *) &cy[ix + *incx]);

            /* pack the inputs into 256-bit registers */
            xmm0 = _mm256_insertf128_pd(xmm0, hxmm1, 0x1);
            ymm0 = _mm256_insertf128_pd(ymm0, hymm1, 0x1);

            /* shuffle the loaded inputs */
            xrmm0 = _mm256_movedup_pd(xmm0);
            ximm0 = _mm256_unpackhi_pd(xmm0, xmm0);
            yrmm0 = _mm256_movedup_pd(ymm0);
            yimm0 = _mm256_unpackhi_pd(ymm0, ymm0);

            /* compute x outputs */
            oxm0 = _mm256_mul_pd(srimm, yimm0);
            oxm0 = _mm256_fmaddsub_pd(sirmm, yrmm0, oxm0);
            oxm0 = _mm256_fmadd_pd(cmm, xmm0, oxm0);

            /* compute y outputs */
            oym0 = _mm256_mul_pd(msrimm, ximm0);
            oym0 = _mm256_fmaddsub_pd(msirmm, xrmm0, oym0);
            oym0 = _mm256_fmsub_pd(cmm, ymm0, oym0);

            /* extract the results */
            hxmm0 = _mm256_extractf128_pd(oxm0, 0x0);
            hxmm1 = _mm256_extractf128_pd(oxm0, 0x1);
            hymm0 = _mm256_extractf128_pd(oym0, 0x0);
            hymm1 = _mm256_extractf128_pd(oym0, 0x1);

            /* store the results */
            _mm_storeu_pd((double *) &cx[ix], hxmm0);
            _mm_storeu_pd((double *) &cx[ix + *incx], hxmm1);
            _mm_storeu_pd((double *) &cy[ix], hymm0);
            _mm_storeu_pd((double *) &cy[ix + *incx], hymm1);

            ix += 2 * *incx;
        }
        for ( ; i__ <= i__1; ++i__)
        {
            /* load complex inputs from x & y */
            xmm  = _mm_loadu_pd((double const *) &cx[ix]);
            ymm  = _mm_loadu_pd((double const *) &cy[ix]);

            /* shuffle the loaded inputs */
            xrmm = _mm_movedup_pd(xmm);
            ximm = _mm_unpackhi_pd(xmm, xmm);
            yrmm = _mm_movedup_pd(ymm);
            yimm = _mm_unpackhi_pd(ymm, ymm);   

            /* compute x outputs */
            oxm = _mm_mul_pd(srim, yimm);
            oxm = _mm_fmaddsub_pd(sirm, yrmm, oxm);
            oxm = _mm_fmadd_pd(cm, xmm, oxm);

            /* compute y outputs */
            oym = _mm_mul_pd(msrim, ximm);
            oym = _mm_fmaddsub_pd(msirm, xrmm, oym);
            oym = _mm_fmsub_pd(cm, ymm, oym);

            /* store the results */
            _mm_storeu_pd((double *) &cx[ix], oxm);
            _mm_storeu_pd((double *) &cy[ix], oym);

            ix += *incx;
        }
    }
    return 0;
    /* Code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= (i__1 - 3); i__ += 4)
    {
        /* load complex inputs from x & y */
        xmm0 = _mm256_loadu_pd((double const *) &cx[i__]);
        ymm0 = _mm256_loadu_pd((double const *) &cy[i__]);
        xmm1 = _mm256_loadu_pd((double const *) &cx[i__ + 2]);
        ymm1 = _mm256_loadu_pd((double const *) &cy[i__ + 2]);

        /* shuffle the loaded inputs */
        xrmm0 = _mm256_movedup_pd(xmm0);
        ximm0 = _mm256_unpackhi_pd(xmm0, xmm0);
        yrmm0 = _mm256_movedup_pd(ymm0);
        yimm0 = _mm256_unpackhi_pd(ymm0, ymm0);

        xrmm1 = _mm256_movedup_pd(xmm1);
        ximm1 = _mm256_unpackhi_pd(xmm1, xmm1);
        yrmm1 = _mm256_movedup_pd(ymm1);
        yimm1 = _mm256_unpackhi_pd(ymm1, ymm1);

        /* compute x outputs */
        oxm0 = _mm256_mul_pd(srimm, yimm0);
        oxm0 = _mm256_fmaddsub_pd(sirmm, yrmm0, oxm0);
        oxm0 = _mm256_fmadd_pd(cmm, xmm0, oxm0);

        oxm1 = _mm256_mul_pd(srimm, yimm1);
        oxm1 = _mm256_fmaddsub_pd(sirmm, yrmm1, oxm1);
        oxm1 = _mm256_fmadd_pd(cmm, xmm1, oxm1);

        /* compute y outputs */
        oym0 = _mm256_mul_pd(msrimm, ximm0);
        oym0 = _mm256_fmaddsub_pd(msirmm, xrmm0, oym0);
        oym0 = _mm256_fmsub_pd(cmm, ymm0, oym0);

        oym1 = _mm256_mul_pd(msrimm, ximm1);
        oym1 = _mm256_fmaddsub_pd(msirmm, xrmm1, oym1);
        oym1 = _mm256_fmsub_pd(cmm, ymm1, oym1);

        /* store the results */
        _mm256_storeu_pd((double *) &cx[i__], oxm0);
        _mm256_storeu_pd((double *) &cy[i__], oym0);
        _mm256_storeu_pd((double *) &cx[i__ + 2], oxm1);
        _mm256_storeu_pd((double *) &cy[i__ + 2], oym1);
    }

    for ( ; i__ <= i__1; ++i__)
    {
        /* load complex inputs from x & y */
        xmm  = _mm_loadu_pd((double const *) &cx[i__]);
        ymm  = _mm_loadu_pd((double const *) &cy[i__]);

        /* shuffle the loaded inputs */
        xrmm = _mm_movedup_pd(xmm);
        ximm = _mm_unpackhi_pd(xmm, xmm);
        yrmm = _mm_movedup_pd(ymm);
        yimm = _mm_unpackhi_pd(ymm, ymm);   

        /* compute x outputs */
        oxm = _mm_mul_pd(srim, yimm);
        oxm = _mm_fmaddsub_pd(sirm, yrmm, oxm);
        oxm = _mm_fmadd_pd(cm, xmm, oxm);

        /* compute y outputs */
        oym = _mm_mul_pd(msrim, ximm);
        oym = _mm_fmaddsub_pd(msirm, xrmm, oym);
        oym = _mm_fmsub_pd(cm, ymm, oym);

        /* store the results */
        _mm_storeu_pd((double *) &cx[i__], oxm);
        _mm_storeu_pd((double *) &cy[i__], oym);
    }

    return 0;
}

int fla_zrot_native(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s)
{
    /* System generated locals */
    integer i__1;
    doublecomplex z__1, z__2, z__3;
    /* Local variables */
    integer i__, ix, iy;
    doublereal lc, sr, si;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --cy;
    --cx;
    /* Function Body */

    if (*n <= 0)
    {
        return 0;
    }
    lc  = *c__;
    sr = s->r;
    si = s->i;

    if (*incx == 1 && *incy == 1)
    {
        goto L20;
    }
    /* Code for unequal increments or equal increments not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }

    i__1 = *n;
    if (*incx != *incy)
    {
        for (i__ = 1; i__ <= i__1; ++i__)
        {
            z__2.r = lc * cx[ix].r;
            z__2.i = lc * cx[ix].i; // , expr subst
            z__3.r = sr * cy[iy].r - si * cy[iy].i;
            z__3.i = sr * cy[iy].i + si * cy[iy].r; // , expr subst
            z__1.r = z__2.r + z__3.r;
            z__1.i = z__2.i + z__3.i; // , expr subst

            z__2.r = lc * cy[iy].r;
            z__2.i = lc * cy[iy].i; // , expr subst
            z__3.r = sr * cx[ix].r + si * cx[ix].i;
            z__3.i = sr * cx[ix].i - si * cx[ix].r; // , expr subst

            cy[iy].r = z__2.r - z__3.r;
            cy[iy].i = z__2.i - z__3.i; // , expr subst
            cx[ix].r = z__1.r;
            cx[ix].i = z__1.i; // , expr subst
            ix += *incx;
            iy += *incy;
        }
    }
    else
    {
        for ( i__ = 1; i__ <= i__1; ++i__)
        {
            z__2.r = lc * cx[ix].r;
            z__2.i = lc * cx[ix].i; // , expr subst
            z__3.r = sr * cy[ix].r - si * cy[ix].i;
            z__3.i = sr * cy[ix].i + si * cy[ix].r; // , expr subst
            z__1.r = z__2.r + z__3.r;
            z__1.i = z__2.i + z__3.i; // , expr subst

            z__2.r = lc * cy[ix].r;
            z__2.i = lc * cy[ix].i; // , expr subst
            z__3.r = sr * cx[ix].r + si * cx[ix].i;
            z__3.i = sr * cx[ix].i - si * cx[ix].r; // , expr subst

            cy[ix].r = z__2.r - z__3.r;
            cy[ix].i = z__2.i - z__3.i; // , expr subst
            cx[ix].r = z__1.r;
            cx[ix].i = z__1.i; // , expr subst
            ix += *incx;
        }
    }
    return 0;
    /* Code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        z__2.r = lc * cx[i__].r;
        z__2.i = lc * cx[i__].i; // , expr subst
        z__3.r = sr * cy[i__].r - si * cy[i__].i;
        z__3.i = sr * cy[i__].i + si * cy[i__].r; // , expr subst
        z__1.r = z__2.r + z__3.r;
        z__1.i = z__2.i + z__3.i; // , expr subst

        z__2.r = lc * cy[i__].r;
        z__2.i = lc * cy[i__].i; // , expr subst
        z__3.r = sr * cx[i__].r + si * cx[i__].i;
        z__3.i = sr * cx[i__].i - si * cx[i__].r; // , expr subst

        cy[i__].r = z__2.r - z__3.r;
        cy[i__].i = z__2.i - z__3.i; // , expr subst
        cx[i__].r = z__1.r;
        cx[i__].i = z__1.i; // , expr subst
    }
    return 0;
}
/* zrot_ */
