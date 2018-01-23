/* ../netlib/clahqr.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c__2 = 2;
/* > \brief \b CLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using th e double-shift/single-shift QR algorithm. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAHQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahqr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahqr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahqr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, */
/* IHIZ, Z, LDZ, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N */
/* LOGICAL WANTT, WANTZ */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX H( LDH, * ), W( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAHQR is an auxiliary routine called by CHSEQR to update the */
/* > eigenvalues and Schur decomposition already computed by CHSEQR, by */
/* > dealing with the Hessenberg submatrix in rows and columns ILO to */
/* > IHI. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] WANTT */
/* > \verbatim */
/* > WANTT is LOGICAL */
/* > = .TRUE. : the full Schur form T is required;
*/
/* > = .FALSE.: only eigenvalues are required. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is LOGICAL */
/* > = .TRUE. : the matrix of Schur vectors Z is required;
*/
/* > = .FALSE.: Schur vectors are not required. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix H. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > It is assumed that H is already upper triangular in rows and */
/* > columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1). */
/* > CLAHQR works primarily with the Hessenberg submatrix in rows */
/* > and columns ILO to IHI, but applies transformations to all of */
/* > H if WANTT is .TRUE.. */
/* > 1 <= ILO <= max(1,IHI);
IHI <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* > H is COMPLEX array, dimension (LDH,N) */
/* > On entry, the upper Hessenberg matrix H. */
/* > On exit, if INFO is zero and if WANTT is .TRUE., then H */
/* > is upper triangular in rows and columns ILO:IHI. If INFO */
/* > is zero and if WANTT is .FALSE., then the contents of H */
/* > are unspecified on exit. The output state of H in case */
/* > INF is positive is below under the description of INFO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is INTEGER */
/* > The leading dimension of the array H. LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is COMPLEX array, dimension (N) */
/* > The computed eigenvalues ILO to IHI are stored in the */
/* > corresponding elements of W. If WANTT is .TRUE., the */
/* > eigenvalues are stored in the same order as on the diagonal */
/* > of the Schur form returned in H, with W(i) = H(i,i). */
/* > \endverbatim */
/* > */
/* > \param[in] ILOZ */
/* > \verbatim */
/* > ILOZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHIZ */
/* > \verbatim */
/* > IHIZ is INTEGER */
/* > Specify the rows of Z to which transformations must be */
/* > applied if WANTZ is .TRUE.. */
/* > 1 <= ILOZ <= ILO;
IHI <= IHIZ <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ,N) */
/* > If WANTZ is .TRUE., on entry Z must contain the current */
/* > matrix Z of transformations accumulated by CHSEQR, and on */
/* > exit Z has been updated;
transformations are applied only to */
/* > the submatrix Z(ILOZ:IHIZ,ILO:IHI). */
/* > If WANTZ is .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > .GT. 0: if INFO = i, CLAHQR failed to compute all the */
/* > eigenvalues ILO to IHI in a total of 30 iterations */
/* > per eigenvalue;
elements i+1:ihi of W contain */
/* > those eigenvalues which have been successfully */
/* > computed. */
/* > */
/* > If INFO .GT. 0 and WANTT is .FALSE., then on exit, */
/* > the remaining unconverged eigenvalues are the */
/* > eigenvalues of the upper Hessenberg matrix */
/* > rows and columns ILO thorugh INFO of the final, */
/* > output value of H. */
/* > */
/* > If INFO .GT. 0 and WANTT is .TRUE., then on exit */
/* > (*) (initial value of H)*U = U*(final value of H) */
/* > where U is an orthognal matrix. The final */
/* > value of H is upper Hessenberg and triangular in */
/* > rows and columns INFO+1 through IHI. */
/* > */
/* > If INFO .GT. 0 and WANTZ is .TRUE., then on exit */
/* > (final value of Z) = (initial value of Z)*U */
/* > where U is the orthogonal matrix in (*) */
/* > (regardless of the value of WANTT.) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > 02-96 Based on modifications by */
/* > David Day, Sandia National Laboratory, USA */
/* > */
/* > 12-04 Further modifications by */
/* > Ralph Byers, University of Kansas, USA */
/* > This is a modified version of CLAHQR from LAPACK version 3.0. */
/* > It is (1) more robust against overflow and underflow and */
/* > (2) adopts the more conservative Ahues & Tisseur stopping */
/* > criterion (LAWN 122, 1997). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int clahqr_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, complex *h__, integer *ldh, complex *w, integer *iloz, integer *ihiz, complex *z__, integer *ldz, integer * info)
{
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7;
    /* Builtin functions */
    double r_imag(complex *);
    void r_cnjg(complex *, complex *);
    double c_abs(complex *);
    void c_sqrt(complex *, complex *), pow_ci(complex *, complex *, integer *) ;
    /* Local variables */
    integer i__, j, k, l, m;
    real s;
    complex t, u, v[2], x, y;
    integer i1, i2;
    complex t1;
    real t2;
    complex v2;
    real aa, ab, ba, bb, h10;
    complex h11;
    real h21;
    complex h22, sc;
    integer nh, nz;
    real sx;
    integer jhi;
    complex h11s;
    integer jlo, its;
    real ulp;
    complex sum;
    real tst;
    complex temp;
    extern /* Subroutine */
    int cscal_(integer *, complex *, complex *, integer *), ccopy_(integer *, complex *, integer *, complex *, integer *);
    real rtemp;
    extern /* Subroutine */
    int slabad_(real *, real *), clarfg_(integer *, complex *, complex *, integer *, complex *);
    extern /* Complex */
    VOID cladiv_(complex *, complex *, complex *);
    extern real slamch_(char *);
    real safmin, safmax, smlnum;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ========================================================= */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    /* Function Body */
    *info = 0;
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    if (*ilo == *ihi)
    {
        i__1 = *ilo;
        i__2 = *ilo + *ilo * h_dim1;
        w[i__1].r = h__[i__2].r;
        w[i__1].i = h__[i__2].i; // , expr subst
        return 0;
    }
    /* ==== clear out the trash ==== */
    i__1 = *ihi - 3;
    for (j = *ilo;
            j <= i__1;
            ++j)
    {
        i__2 = j + 2 + j * h_dim1;
        h__[i__2].r = 0.f;
        h__[i__2].i = 0.f; // , expr subst
        i__2 = j + 3 + j * h_dim1;
        h__[i__2].r = 0.f;
        h__[i__2].i = 0.f; // , expr subst
        /* L10: */
    }
    if (*ilo <= *ihi - 2)
    {
        i__1 = *ihi + (*ihi - 2) * h_dim1;
        h__[i__1].r = 0.f;
        h__[i__1].i = 0.f; // , expr subst
    }
    /* ==== ensure that subdiagonal entries are real ==== */
    if (*wantt)
    {
        jlo = 1;
        jhi = *n;
    }
    else
    {
        jlo = *ilo;
        jhi = *ihi;
    }
    i__1 = *ihi;
    for (i__ = *ilo + 1;
            i__ <= i__1;
            ++i__)
    {
        if (r_imag(&h__[i__ + (i__ - 1) * h_dim1]) != 0.f)
        {
            /* ==== The following redundant normalization */
            /* . avoids problems with both gradual and */
            /* . sudden underflow in ABS(H(I,I-1)) ==== */
            i__2 = i__ + (i__ - 1) * h_dim1;
            i__3 = i__ + (i__ - 1) * h_dim1;
            r__3 = (r__1 = h__[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[i__ + (i__ - 1) * h_dim1]), f2c_abs(r__2));
            q__1.r = h__[i__2].r / r__3;
            q__1.i = h__[i__2].i / r__3; // , expr subst
            sc.r = q__1.r;
            sc.i = q__1.i; // , expr subst
            r_cnjg(&q__2, &sc);
            r__1 = c_abs(&sc);
            q__1.r = q__2.r / r__1;
            q__1.i = q__2.i / r__1; // , expr subst
            sc.r = q__1.r;
            sc.i = q__1.i; // , expr subst
            i__2 = i__ + (i__ - 1) * h_dim1;
            r__1 = c_abs(&h__[i__ + (i__ - 1) * h_dim1]);
            h__[i__2].r = r__1;
            h__[i__2].i = 0.f; // , expr subst
            i__2 = jhi - i__ + 1;
            cscal_(&i__2, &sc, &h__[i__ + i__ * h_dim1], ldh);
            /* Computing MIN */
            i__3 = jhi;
            i__4 = i__ + 1; // , expr subst
            i__2 = min(i__3,i__4) - jlo + 1;
            r_cnjg(&q__1, &sc);
            cscal_(&i__2, &q__1, &h__[jlo + i__ * h_dim1], &c__1);
            if (*wantz)
            {
                i__2 = *ihiz - *iloz + 1;
                r_cnjg(&q__1, &sc);
                cscal_(&i__2, &q__1, &z__[*iloz + i__ * z_dim1], &c__1);
            }
        }
        /* L20: */
    }
    nh = *ihi - *ilo + 1;
    nz = *ihiz - *iloz + 1;
    /* Set machine-dependent constants for the stopping criterion. */
    safmin = slamch_("SAFE MINIMUM");
    safmax = 1.f / safmin;
    slabad_(&safmin, &safmax);
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real) nh / ulp);
    /* I1 and I2 are the indices of the first row and last column of H */
    /* to which transformations must be applied. If eigenvalues only are */
    /* being computed, I1 and I2 are set inside the main loop. */
    if (*wantt)
    {
        i1 = 1;
        i2 = *n;
    }
    /* The main loop begins here. I is the loop index and decreases from */
    /* IHI to ILO in steps of 1. Each iteration of the loop works */
    /* with the active submatrix in rows and columns L to I. */
    /* Eigenvalues I+1 to IHI have already converged. Either L = ILO, or */
    /* H(L,L-1) is negligible so that the matrix splits. */
    i__ = *ihi;
L30:
    if (i__ < *ilo)
    {
        goto L150;
    }
    /* Perform QR iterations on rows and columns ILO to I until a */
    /* submatrix of order 1 splits off at the bottom because a */
    /* subdiagonal element has become negligible. */
    l = *ilo;
    for (its = 0;
            its <= 30;
            ++its)
    {
        /* Look for a single small subdiagonal element. */
        i__1 = l + 1;
        for (k = i__;
                k >= i__1;
                --k)
        {
            i__2 = k + (k - 1) * h_dim1;
            if ((r__1 = h__[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + (k - 1) * h_dim1]), f2c_abs(r__2)) <= smlnum)
            {
                goto L50;
            }
            i__2 = k - 1 + (k - 1) * h_dim1;
            i__3 = k + k * h_dim1;
            tst = (r__1 = h__[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k - 1 + (k - 1) * h_dim1]), f2c_abs(r__2)) + ((r__3 = h__[i__3].r, f2c_abs(r__3)) + (r__4 = r_imag(&h__[k + k * h_dim1]), f2c_abs( r__4)));
            if (tst == 0.f)
            {
                if (k - 2 >= *ilo)
                {
                    i__2 = k - 1 + (k - 2) * h_dim1;
                    tst += (r__1 = h__[i__2].r, f2c_abs(r__1));
                }
                if (k + 1 <= *ihi)
                {
                    i__2 = k + 1 + k * h_dim1;
                    tst += (r__1 = h__[i__2].r, f2c_abs(r__1));
                }
            }
            /* ==== The following is a conservative small subdiagonal */
            /* . deflation criterion due to Ahues & Tisseur (LAWN 122, */
            /* . 1997). It has better mathematical foundation and */
            /* . improves accuracy in some examples. ==== */
            i__2 = k + (k - 1) * h_dim1;
            if ((r__1 = h__[i__2].r, f2c_abs(r__1)) <= ulp * tst)
            {
                /* Computing MAX */
                i__2 = k + (k - 1) * h_dim1;
                i__3 = k - 1 + k * h_dim1;
                r__5 = (r__1 = h__[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[ k + (k - 1) * h_dim1]), f2c_abs(r__2));
                r__6 = (r__3 = h__[i__3].r, f2c_abs(r__3)) + (r__4 = r_imag(&h__[k - 1 + k * h_dim1]), f2c_abs(r__4)); // , expr subst
                ab = max(r__5,r__6);
                /* Computing MIN */
                i__2 = k + (k - 1) * h_dim1;
                i__3 = k - 1 + k * h_dim1;
                r__5 = (r__1 = h__[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[ k + (k - 1) * h_dim1]), f2c_abs(r__2));
                r__6 = (r__3 = h__[i__3].r, f2c_abs(r__3)) + (r__4 = r_imag(&h__[k - 1 + k * h_dim1]), f2c_abs(r__4)); // , expr subst
                ba = min(r__5,r__6);
                i__2 = k - 1 + (k - 1) * h_dim1;
                i__3 = k + k * h_dim1;
                q__2.r = h__[i__2].r - h__[i__3].r;
                q__2.i = h__[i__2].i - h__[i__3].i; // , expr subst
                q__1.r = q__2.r;
                q__1.i = q__2.i; // , expr subst
                /* Computing MAX */
                i__4 = k + k * h_dim1;
                r__5 = (r__1 = h__[i__4].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[ k + k * h_dim1]), f2c_abs(r__2));
                r__6 = (r__3 = q__1.r, f2c_abs(r__3)) + (r__4 = r_imag(&q__1), f2c_abs(r__4)); // , expr subst
                aa = max(r__5,r__6);
                i__2 = k - 1 + (k - 1) * h_dim1;
                i__3 = k + k * h_dim1;
                q__2.r = h__[i__2].r - h__[i__3].r;
                q__2.i = h__[i__2].i - h__[i__3].i; // , expr subst
                q__1.r = q__2.r;
                q__1.i = q__2.i; // , expr subst
                /* Computing MIN */
                i__4 = k + k * h_dim1;
                r__5 = (r__1 = h__[i__4].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[ k + k * h_dim1]), f2c_abs(r__2));
                r__6 = (r__3 = q__1.r, f2c_abs(r__3)) + (r__4 = r_imag(&q__1), f2c_abs(r__4)); // , expr subst
                bb = min(r__5,r__6);
                s = aa + ab;
                /* Computing MAX */
                r__1 = smlnum;
                r__2 = ulp * (bb * (aa / s)); // , expr subst
                if (ba * (ab / s) <= max(r__1,r__2))
                {
                    goto L50;
                }
            }
            /* L40: */
        }
L50:
        l = k;
        if (l > *ilo)
        {
            /* H(L,L-1) is negligible */
            i__1 = l + (l - 1) * h_dim1;
            h__[i__1].r = 0.f;
            h__[i__1].i = 0.f; // , expr subst
        }
        /* Exit from loop if a submatrix of order 1 has split off. */
        if (l >= i__)
        {
            goto L140;
        }
        /* Now the active submatrix is in rows and columns L to I. If */
        /* eigenvalues only are being computed, only the active submatrix */
        /* need be transformed. */
        if (! (*wantt))
        {
            i1 = l;
            i2 = i__;
        }
        if (its == 10)
        {
            /* Exceptional shift. */
            i__1 = l + 1 + l * h_dim1;
            s = (r__1 = h__[i__1].r, f2c_abs(r__1)) * .75f;
            i__1 = l + l * h_dim1;
            q__1.r = s + h__[i__1].r;
            q__1.i = h__[i__1].i; // , expr subst
            t.r = q__1.r;
            t.i = q__1.i; // , expr subst
        }
        else if (its == 20)
        {
            /* Exceptional shift. */
            i__1 = i__ + (i__ - 1) * h_dim1;
            s = (r__1 = h__[i__1].r, f2c_abs(r__1)) * .75f;
            i__1 = i__ + i__ * h_dim1;
            q__1.r = s + h__[i__1].r;
            q__1.i = h__[i__1].i; // , expr subst
            t.r = q__1.r;
            t.i = q__1.i; // , expr subst
        }
        else
        {
            /* Wilkinson's shift. */
            i__1 = i__ + i__ * h_dim1;
            t.r = h__[i__1].r;
            t.i = h__[i__1].i; // , expr subst
            c_sqrt(&q__2, &h__[i__ - 1 + i__ * h_dim1]);
            c_sqrt(&q__3, &h__[i__ + (i__ - 1) * h_dim1]);
            q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
            q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
            u.r = q__1.r;
            u.i = q__1.i; // , expr subst
            s = (r__1 = u.r, f2c_abs(r__1)) + (r__2 = r_imag(&u), f2c_abs(r__2));
            if (s != 0.f)
            {
                i__1 = i__ - 1 + (i__ - 1) * h_dim1;
                q__2.r = h__[i__1].r - t.r;
                q__2.i = h__[i__1].i - t.i; // , expr subst
                q__1.r = q__2.r * .5f;
                q__1.i = q__2.i * .5f; // , expr subst
                x.r = q__1.r;
                x.i = q__1.i; // , expr subst
                sx = (r__1 = x.r, f2c_abs(r__1)) + (r__2 = r_imag(&x), f2c_abs(r__2));
                /* Computing MAX */
                r__3 = s;
                r__4 = (r__1 = x.r, f2c_abs(r__1)) + (r__2 = r_imag(&x), f2c_abs(r__2)); // , expr subst
                s = max(r__3,r__4);
                q__5.r = x.r / s;
                q__5.i = x.i / s; // , expr subst
                pow_ci(&q__4, &q__5, &c__2);
                q__7.r = u.r / s;
                q__7.i = u.i / s; // , expr subst
                pow_ci(&q__6, &q__7, &c__2);
                q__3.r = q__4.r + q__6.r;
                q__3.i = q__4.i + q__6.i; // , expr subst
                c_sqrt(&q__2, &q__3);
                q__1.r = s * q__2.r;
                q__1.i = s * q__2.i; // , expr subst
                y.r = q__1.r;
                y.i = q__1.i; // , expr subst
                if (sx > 0.f)
                {
                    q__1.r = x.r / sx;
                    q__1.i = x.i / sx; // , expr subst
                    q__2.r = x.r / sx;
                    q__2.i = x.i / sx; // , expr subst
                    if (q__1.r * y.r + r_imag(&q__2) * r_imag(&y) < 0.f)
                    {
                        q__3.r = -y.r;
                        q__3.i = -y.i; // , expr subst
                        y.r = q__3.r;
                        y.i = q__3.i; // , expr subst
                    }
                }
                q__4.r = x.r + y.r;
                q__4.i = x.i + y.i; // , expr subst
                cladiv_(&q__3, &u, &q__4);
                q__2.r = u.r * q__3.r - u.i * q__3.i;
                q__2.i = u.r * q__3.i + u.i * q__3.r; // , expr subst
                q__1.r = t.r - q__2.r;
                q__1.i = t.i - q__2.i; // , expr subst
                t.r = q__1.r;
                t.i = q__1.i; // , expr subst
            }
        }
        /* Look for two consecutive small subdiagonal elements. */
        i__1 = l + 1;
        for (m = i__ - 1;
                m >= i__1;
                --m)
        {
            /* Determine the effect of starting the single-shift QR */
            /* iteration at row M, and see if this would make H(M,M-1) */
            /* negligible. */
            i__2 = m + m * h_dim1;
            h11.r = h__[i__2].r;
            h11.i = h__[i__2].i; // , expr subst
            i__2 = m + 1 + (m + 1) * h_dim1;
            h22.r = h__[i__2].r;
            h22.i = h__[i__2].i; // , expr subst
            q__1.r = h11.r - t.r;
            q__1.i = h11.i - t.i; // , expr subst
            h11s.r = q__1.r;
            h11s.i = q__1.i; // , expr subst
            i__2 = m + 1 + m * h_dim1;
            h21 = h__[i__2].r;
            s = (r__1 = h11s.r, f2c_abs(r__1)) + (r__2 = r_imag(&h11s), f2c_abs(r__2)) + f2c_abs(h21);
            q__1.r = h11s.r / s;
            q__1.i = h11s.i / s; // , expr subst
            h11s.r = q__1.r;
            h11s.i = q__1.i; // , expr subst
            h21 /= s;
            v[0].r = h11s.r;
            v[0].i = h11s.i; // , expr subst
            v[1].r = h21;
            v[1].i = 0.f; // , expr subst
            i__2 = m + (m - 1) * h_dim1;
            h10 = h__[i__2].r;
            if (f2c_abs(h10) * f2c_abs(h21) <= ulp * (((r__1 = h11s.r, f2c_abs(r__1)) + ( r__2 = r_imag(&h11s), f2c_abs(r__2))) * ((r__3 = h11.r, f2c_abs( r__3)) + (r__4 = r_imag(&h11), f2c_abs(r__4)) + ((r__5 = h22.r, f2c_abs(r__5)) + (r__6 = r_imag(&h22), f2c_abs(r__6))))))
            {
                goto L70;
            }
            /* L60: */
        }
        i__1 = l + l * h_dim1;
        h11.r = h__[i__1].r;
        h11.i = h__[i__1].i; // , expr subst
        i__1 = l + 1 + (l + 1) * h_dim1;
        h22.r = h__[i__1].r;
        h22.i = h__[i__1].i; // , expr subst
        q__1.r = h11.r - t.r;
        q__1.i = h11.i - t.i; // , expr subst
        h11s.r = q__1.r;
        h11s.i = q__1.i; // , expr subst
        i__1 = l + 1 + l * h_dim1;
        h21 = h__[i__1].r;
        s = (r__1 = h11s.r, f2c_abs(r__1)) + (r__2 = r_imag(&h11s), f2c_abs(r__2)) + f2c_abs(h21);
        q__1.r = h11s.r / s;
        q__1.i = h11s.i / s; // , expr subst
        h11s.r = q__1.r;
        h11s.i = q__1.i; // , expr subst
        h21 /= s;
        v[0].r = h11s.r;
        v[0].i = h11s.i; // , expr subst
        v[1].r = h21;
        v[1].i = 0.f; // , expr subst
L70: /* Single-shift QR step */
        i__1 = i__ - 1;
        for (k = m;
                k <= i__1;
                ++k)
        {
            /* The first iteration of this loop determines a reflection G */
            /* from the vector V and applies it from left and right to H, */
            /* thus creating a nonzero bulge below the subdiagonal. */
            /* Each subsequent iteration determines a reflection G to */
            /* restore the Hessenberg form in the (K-1)th column, and thus */
            /* chases the bulge one step toward the bottom of the active */
            /* submatrix. */
            /* V(2) is always real before the call to CLARFG, and hence */
            /* after the call T2 ( = T1*V(2) ) is also real. */
            if (k > m)
            {
                ccopy_(&c__2, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
            }
            clarfg_(&c__2, v, &v[1], &c__1, &t1);
            if (k > m)
            {
                i__2 = k + (k - 1) * h_dim1;
                h__[i__2].r = v[0].r;
                h__[i__2].i = v[0].i; // , expr subst
                i__2 = k + 1 + (k - 1) * h_dim1;
                h__[i__2].r = 0.f;
                h__[i__2].i = 0.f; // , expr subst
            }
            v2.r = v[1].r;
            v2.i = v[1].i; // , expr subst
            q__1.r = t1.r * v2.r - t1.i * v2.i;
            q__1.i = t1.r * v2.i + t1.i * v2.r; // , expr subst
            t2 = q__1.r;
            /* Apply G from the left to transform the rows of the matrix */
            /* in columns K to I2. */
            i__2 = i2;
            for (j = k;
                    j <= i__2;
                    ++j)
            {
                r_cnjg(&q__3, &t1);
                i__3 = k + j * h_dim1;
                q__2.r = q__3.r * h__[i__3].r - q__3.i * h__[i__3].i;
                q__2.i = q__3.r * h__[i__3].i + q__3.i * h__[i__3].r; // , expr subst
                i__4 = k + 1 + j * h_dim1;
                q__4.r = t2 * h__[i__4].r;
                q__4.i = t2 * h__[i__4].i; // , expr subst
                q__1.r = q__2.r + q__4.r;
                q__1.i = q__2.i + q__4.i; // , expr subst
                sum.r = q__1.r;
                sum.i = q__1.i; // , expr subst
                i__3 = k + j * h_dim1;
                i__4 = k + j * h_dim1;
                q__1.r = h__[i__4].r - sum.r;
                q__1.i = h__[i__4].i - sum.i; // , expr subst
                h__[i__3].r = q__1.r;
                h__[i__3].i = q__1.i; // , expr subst
                i__3 = k + 1 + j * h_dim1;
                i__4 = k + 1 + j * h_dim1;
                q__2.r = sum.r * v2.r - sum.i * v2.i;
                q__2.i = sum.r * v2.i + sum.i * v2.r; // , expr subst
                q__1.r = h__[i__4].r - q__2.r;
                q__1.i = h__[i__4].i - q__2.i; // , expr subst
                h__[i__3].r = q__1.r;
                h__[i__3].i = q__1.i; // , expr subst
                /* L80: */
            }
            /* Apply G from the right to transform the columns of the */
            /* matrix in rows I1 to min(K+2,I). */
            /* Computing MIN */
            i__3 = k + 2;
            i__2 = min(i__3,i__);
            for (j = i1;
                    j <= i__2;
                    ++j)
            {
                i__3 = j + k * h_dim1;
                q__2.r = t1.r * h__[i__3].r - t1.i * h__[i__3].i;
                q__2.i = t1.r * h__[i__3].i + t1.i * h__[i__3].r; // , expr subst
                i__4 = j + (k + 1) * h_dim1;
                q__3.r = t2 * h__[i__4].r;
                q__3.i = t2 * h__[i__4].i; // , expr subst
                q__1.r = q__2.r + q__3.r;
                q__1.i = q__2.i + q__3.i; // , expr subst
                sum.r = q__1.r;
                sum.i = q__1.i; // , expr subst
                i__3 = j + k * h_dim1;
                i__4 = j + k * h_dim1;
                q__1.r = h__[i__4].r - sum.r;
                q__1.i = h__[i__4].i - sum.i; // , expr subst
                h__[i__3].r = q__1.r;
                h__[i__3].i = q__1.i; // , expr subst
                i__3 = j + (k + 1) * h_dim1;
                i__4 = j + (k + 1) * h_dim1;
                r_cnjg(&q__3, &v2);
                q__2.r = sum.r * q__3.r - sum.i * q__3.i;
                q__2.i = sum.r * q__3.i + sum.i * q__3.r; // , expr subst
                q__1.r = h__[i__4].r - q__2.r;
                q__1.i = h__[i__4].i - q__2.i; // , expr subst
                h__[i__3].r = q__1.r;
                h__[i__3].i = q__1.i; // , expr subst
                /* L90: */
            }
            if (*wantz)
            {
                /* Accumulate transformations in the matrix Z */
                i__2 = *ihiz;
                for (j = *iloz;
                        j <= i__2;
                        ++j)
                {
                    i__3 = j + k * z_dim1;
                    q__2.r = t1.r * z__[i__3].r - t1.i * z__[i__3].i;
                    q__2.i = t1.r * z__[i__3].i + t1.i * z__[i__3].r; // , expr subst
                    i__4 = j + (k + 1) * z_dim1;
                    q__3.r = t2 * z__[i__4].r;
                    q__3.i = t2 * z__[i__4].i; // , expr subst
                    q__1.r = q__2.r + q__3.r;
                    q__1.i = q__2.i + q__3.i; // , expr subst
                    sum.r = q__1.r;
                    sum.i = q__1.i; // , expr subst
                    i__3 = j + k * z_dim1;
                    i__4 = j + k * z_dim1;
                    q__1.r = z__[i__4].r - sum.r;
                    q__1.i = z__[i__4].i - sum.i; // , expr subst
                    z__[i__3].r = q__1.r;
                    z__[i__3].i = q__1.i; // , expr subst
                    i__3 = j + (k + 1) * z_dim1;
                    i__4 = j + (k + 1) * z_dim1;
                    r_cnjg(&q__3, &v2);
                    q__2.r = sum.r * q__3.r - sum.i * q__3.i;
                    q__2.i = sum.r * q__3.i + sum.i * q__3.r; // , expr subst
                    q__1.r = z__[i__4].r - q__2.r;
                    q__1.i = z__[i__4].i - q__2.i; // , expr subst
                    z__[i__3].r = q__1.r;
                    z__[i__3].i = q__1.i; // , expr subst
                    /* L100: */
                }
            }
            if (k == m && m > l)
            {
                /* If the QR step was started at row M > L because two */
                /* consecutive small subdiagonals were found, then extra */
                /* scaling must be performed to ensure that H(M,M-1) remains */
                /* real. */
                q__1.r = 1.f - t1.r;
                q__1.i = 0.f - t1.i; // , expr subst
                temp.r = q__1.r;
                temp.i = q__1.i; // , expr subst
                r__1 = c_abs(&temp);
                q__1.r = temp.r / r__1;
                q__1.i = temp.i / r__1; // , expr subst
                temp.r = q__1.r;
                temp.i = q__1.i; // , expr subst
                i__2 = m + 1 + m * h_dim1;
                i__3 = m + 1 + m * h_dim1;
                r_cnjg(&q__2, &temp);
                q__1.r = h__[i__3].r * q__2.r - h__[i__3].i * q__2.i;
                q__1.i = h__[i__3].r * q__2.i + h__[i__3].i * q__2.r; // , expr subst
                h__[i__2].r = q__1.r;
                h__[i__2].i = q__1.i; // , expr subst
                if (m + 2 <= i__)
                {
                    i__2 = m + 2 + (m + 1) * h_dim1;
                    i__3 = m + 2 + (m + 1) * h_dim1;
                    q__1.r = h__[i__3].r * temp.r - h__[i__3].i * temp.i;
                    q__1.i = h__[i__3].r * temp.i + h__[i__3].i * temp.r; // , expr subst
                    h__[i__2].r = q__1.r;
                    h__[i__2].i = q__1.i; // , expr subst
                }
                i__2 = i__;
                for (j = m;
                        j <= i__2;
                        ++j)
                {
                    if (j != m + 1)
                    {
                        if (i2 > j)
                        {
                            i__3 = i2 - j;
                            cscal_(&i__3, &temp, &h__[j + (j + 1) * h_dim1], ldh);
                        }
                        i__3 = j - i1;
                        r_cnjg(&q__1, &temp);
                        cscal_(&i__3, &q__1, &h__[i1 + j * h_dim1], &c__1);
                        if (*wantz)
                        {
                            r_cnjg(&q__1, &temp);
                            cscal_(&nz, &q__1, &z__[*iloz + j * z_dim1], & c__1);
                        }
                    }
                    /* L110: */
                }
            }
            /* L120: */
        }
        /* Ensure that H(I,I-1) is real. */
        i__1 = i__ + (i__ - 1) * h_dim1;
        temp.r = h__[i__1].r;
        temp.i = h__[i__1].i; // , expr subst
        if (r_imag(&temp) != 0.f)
        {
            rtemp = c_abs(&temp);
            i__1 = i__ + (i__ - 1) * h_dim1;
            h__[i__1].r = rtemp;
            h__[i__1].i = 0.f; // , expr subst
            q__1.r = temp.r / rtemp;
            q__1.i = temp.i / rtemp; // , expr subst
            temp.r = q__1.r;
            temp.i = q__1.i; // , expr subst
            if (i2 > i__)
            {
                i__1 = i2 - i__;
                r_cnjg(&q__1, &temp);
                cscal_(&i__1, &q__1, &h__[i__ + (i__ + 1) * h_dim1], ldh);
            }
            i__1 = i__ - i1;
            cscal_(&i__1, &temp, &h__[i1 + i__ * h_dim1], &c__1);
            if (*wantz)
            {
                cscal_(&nz, &temp, &z__[*iloz + i__ * z_dim1], &c__1);
            }
        }
        /* L130: */
    }
    /* Failure to converge in remaining number of iterations */
    *info = i__;
    return 0;
L140: /* H(I,I-1) is negligible: one eigenvalue has converged. */
    i__1 = i__;
    i__2 = i__ + i__ * h_dim1;
    w[i__1].r = h__[i__2].r;
    w[i__1].i = h__[i__2].i; // , expr subst
    /* return to start of the main loop with new value of I. */
    i__ = l - 1;
    goto L30;
L150:
    return 0;
    /* End of CLAHQR */
}
/* clahqr_ */
