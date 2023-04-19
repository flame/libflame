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
      retval = fla_zrot(n, cx, incx, cy, incy, c__, s);
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
