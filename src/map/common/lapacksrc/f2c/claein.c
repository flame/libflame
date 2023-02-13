/* ../netlib/claein.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse iteration. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claein. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claein. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claein. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAEIN( RIGHTV, NOINIT, N, H, LDH, W, V, B, LDB, RWORK, */
/* EPS3, SMLNUM, INFO ) */
/* .. Scalar Arguments .. */
/* LOGICAL NOINIT, RIGHTV */
/* INTEGER INFO, LDB, LDH, N */
/* REAL EPS3, SMLNUM */
/* COMPLEX W */
/* .. */
/* .. Array Arguments .. */
/* REAL RWORK( * ) */
/* COMPLEX B( LDB, * ), H( LDH, * ), V( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAEIN uses inverse iteration to find a right or left eigenvector */
/* > corresponding to the eigenvalue W of a complex upper Hessenberg */
/* > matrix H. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] RIGHTV */
/* > \verbatim */
/* > RIGHTV is LOGICAL */
/* > = .TRUE. : compute right eigenvector;
*/
/* > = .FALSE.: compute left eigenvector. */
/* > \endverbatim */
/* > */
/* > \param[in] NOINIT */
/* > \verbatim */
/* > NOINIT is LOGICAL */
/* > = .TRUE. : no initial vector supplied in V */
/* > = .FALSE.: initial vector supplied in V. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix H. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* > H is COMPLEX array, dimension (LDH,N) */
/* > The upper Hessenberg matrix H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is INTEGER */
/* > The leading dimension of the array H. LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* > W is COMPLEX */
/* > The eigenvalue of H whose corresponding right or left */
/* > eigenvector is to be computed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* > V is COMPLEX array, dimension (N) */
/* > On entry, if NOINIT = .FALSE., V must contain a starting */
/* > vector for inverse iteration;
otherwise V need not be set. */
/* > On exit, V contains the computed eigenvector, normalized so */
/* > that the component of largest magnitude has magnitude 1;
here */
/* > the magnitude of a complex number (x,y) is taken to be */
/* > |x| + |y|. */
/* > \endverbatim */
/* > */
/* > \param[out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in] EPS3 */
/* > \verbatim */
/* > EPS3 is REAL */
/* > A small machine-dependent value which is used to perturb */
/* > close eigenvalues, and to replace zero pivots. */
/* > \endverbatim */
/* > */
/* > \param[in] SMLNUM */
/* > \verbatim */
/* > SMLNUM is REAL */
/* > A machine-dependent value close to the underflow threshold. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > = 1: inverse iteration did not converge;
V is set to the */
/* > last iterate. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int claein_(logical *rightv, logical *noinit, integer *n, complex *h__, integer *ldh, complex *w, complex *v, complex *b, integer *ldb, real *rwork, real *eps3, real *smlnum, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, h_dim1, h_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2;
    /* Builtin functions */
    double sqrt(doublereal), r_imag(complex *);
    /* Local variables */
    integer i__, j;
    complex x, ei, ej;
    integer its, ierr;
    complex temp;
    real scale;
    char trans[1];
    real rtemp, rootn, vnorm;
    extern real scnrm2_(integer *, complex *, integer *);
    extern integer icamax_(integer *, complex *, integer *);
    extern /* Complex */
    VOID cladiv_(complex *, complex *, complex *);
    extern /* Subroutine */
    int csscal_(integer *, real *, complex *, integer *), clatrs_(char *, char *, char *, char *, integer *, complex *, integer *, complex *, real *, real *, integer *);
    extern real scasum_(integer *, complex *, integer *);
    char normin[1];
    real nrmsml, growto;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --v;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --rwork;
    /* Function Body */
    *info = 0;
    /* GROWTO is the threshold used in the acceptance test for an */
    /* eigenvector. */
    rootn = sqrt((real) (*n));
    growto = .1f / rootn;
    /* Computing MAX */
    r__1 = 1.f;
    r__2 = *eps3 * rootn; // , expr subst
    nrmsml = max(r__1,r__2) * *smlnum;
    /* Form B = H - W*I (except that the subdiagonal elements are not */
    /* stored). */
    i__1 = *n;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        i__2 = j - 1;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            i__3 = i__ + j * b_dim1;
            i__4 = i__ + j * h_dim1;
            b[i__3].r = h__[i__4].r;
            b[i__3].i = h__[i__4].i; // , expr subst
            /* L10: */
        }
        i__2 = j + j * b_dim1;
        i__3 = j + j * h_dim1;
        q__1.r = h__[i__3].r - w->r;
        q__1.i = h__[i__3].i - w->i; // , expr subst
        b[i__2].r = q__1.r;
        b[i__2].i = q__1.i; // , expr subst
        /* L20: */
    }
    if (*noinit)
    {
        /* Initialize V. */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__;
            v[i__2].r = *eps3;
            v[i__2].i = 0.f; // , expr subst
            /* L30: */
        }
    }
    else
    {
        /* Scale supplied initial vector. */
        vnorm = scnrm2_(n, &v[1], &c__1);
        r__1 = *eps3 * rootn / max(vnorm,nrmsml);
        csscal_(n, &r__1, &v[1], &c__1);
    }
    if (*rightv)
    {
        /* LU decomposition with partial pivoting of B, replacing zero */
        /* pivots by EPS3. */
        i__1 = *n - 1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__ + 1 + i__ * h_dim1;
            ei.r = h__[i__2].r;
            ei.i = h__[i__2].i; // , expr subst
            i__2 = i__ + i__ * b_dim1;
            if ((r__1 = b[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&b[i__ + i__ * b_dim1]), f2c_abs(r__2)) < (r__3 = ei.r, f2c_abs(r__3)) + (r__4 = r_imag(&ei), f2c_abs(r__4)))
            {
                /* Interchange rows and eliminate. */
                cladiv_(&q__1, &b[i__ + i__ * b_dim1], &ei);
                x.r = q__1.r;
                x.i = q__1.i; // , expr subst
                i__2 = i__ + i__ * b_dim1;
                b[i__2].r = ei.r;
                b[i__2].i = ei.i; // , expr subst
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = i__ + 1 + j * b_dim1;
                    temp.r = b[i__3].r;
                    temp.i = b[i__3].i; // , expr subst
                    i__3 = i__ + 1 + j * b_dim1;
                    i__4 = i__ + j * b_dim1;
                    q__2.r = x.r * temp.r - x.i * temp.i;
                    q__2.i = x.r * temp.i + x.i * temp.r; // , expr subst
                    q__1.r = b[i__4].r - q__2.r;
                    q__1.i = b[i__4].i - q__2.i; // , expr subst
                    b[i__3].r = q__1.r;
                    b[i__3].i = q__1.i; // , expr subst
                    i__3 = i__ + j * b_dim1;
                    b[i__3].r = temp.r;
                    b[i__3].i = temp.i; // , expr subst
                    /* L40: */
                }
            }
            else
            {
                /* Eliminate without interchange. */
                i__2 = i__ + i__ * b_dim1;
                if (b[i__2].r == 0.f && b[i__2].i == 0.f)
                {
                    i__3 = i__ + i__ * b_dim1;
                    b[i__3].r = *eps3;
                    b[i__3].i = 0.f; // , expr subst
                }
                cladiv_(&q__1, &ei, &b[i__ + i__ * b_dim1]);
                x.r = q__1.r;
                x.i = q__1.i; // , expr subst
                if (x.r != 0.f || x.i != 0.f)
                {
                    i__2 = *n;
                    for (j = i__ + 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = i__ + 1 + j * b_dim1;
                        i__4 = i__ + 1 + j * b_dim1;
                        i__5 = i__ + j * b_dim1;
                        q__2.r = x.r * b[i__5].r - x.i * b[i__5].i;
                        q__2.i = x.r * b[i__5].i + x.i * b[i__5].r; // , expr subst
                        q__1.r = b[i__4].r - q__2.r;
                        q__1.i = b[i__4].i - q__2.i; // , expr subst
                        b[i__3].r = q__1.r;
                        b[i__3].i = q__1.i; // , expr subst
                        /* L50: */
                    }
                }
            }
            /* L60: */
        }
        i__1 = *n + *n * b_dim1;
        if (b[i__1].r == 0.f && b[i__1].i == 0.f)
        {
            i__2 = *n + *n * b_dim1;
            b[i__2].r = *eps3;
            b[i__2].i = 0.f; // , expr subst
        }
        *(unsigned char *)trans = 'N';
    }
    else
    {
        /* UL decomposition with partial pivoting of B, replacing zero */
        /* pivots by EPS3. */
        for (j = *n;
                j >= 2;
                --j)
        {
            i__1 = j + (j - 1) * h_dim1;
            ej.r = h__[i__1].r;
            ej.i = h__[i__1].i; // , expr subst
            i__1 = j + j * b_dim1;
            if ((r__1 = b[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&b[j + j * b_dim1]), f2c_abs(r__2)) < (r__3 = ej.r, f2c_abs(r__3)) + (r__4 = r_imag(&ej), f2c_abs(r__4)))
            {
                /* Interchange columns and eliminate. */
                cladiv_(&q__1, &b[j + j * b_dim1], &ej);
                x.r = q__1.r;
                x.i = q__1.i; // , expr subst
                i__1 = j + j * b_dim1;
                b[i__1].r = ej.r;
                b[i__1].i = ej.i; // , expr subst
                i__1 = j - 1;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = i__ + (j - 1) * b_dim1;
                    temp.r = b[i__2].r;
                    temp.i = b[i__2].i; // , expr subst
                    i__2 = i__ + (j - 1) * b_dim1;
                    i__3 = i__ + j * b_dim1;
                    q__2.r = x.r * temp.r - x.i * temp.i;
                    q__2.i = x.r * temp.i + x.i * temp.r; // , expr subst
                    q__1.r = b[i__3].r - q__2.r;
                    q__1.i = b[i__3].i - q__2.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = i__ + j * b_dim1;
                    b[i__2].r = temp.r;
                    b[i__2].i = temp.i; // , expr subst
                    /* L70: */
                }
            }
            else
            {
                /* Eliminate without interchange. */
                i__1 = j + j * b_dim1;
                if (b[i__1].r == 0.f && b[i__1].i == 0.f)
                {
                    i__2 = j + j * b_dim1;
                    b[i__2].r = *eps3;
                    b[i__2].i = 0.f; // , expr subst
                }
                cladiv_(&q__1, &ej, &b[j + j * b_dim1]);
                x.r = q__1.r;
                x.i = q__1.i; // , expr subst
                if (x.r != 0.f || x.i != 0.f)
                {
                    i__1 = j - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        i__2 = i__ + (j - 1) * b_dim1;
                        i__3 = i__ + (j - 1) * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        q__2.r = x.r * b[i__4].r - x.i * b[i__4].i;
                        q__2.i = x.r * b[i__4].i + x.i * b[i__4].r; // , expr subst
                        q__1.r = b[i__3].r - q__2.r;
                        q__1.i = b[i__3].i - q__2.i; // , expr subst
                        b[i__2].r = q__1.r;
                        b[i__2].i = q__1.i; // , expr subst
                        /* L80: */
                    }
                }
            }
            /* L90: */
        }
        i__1 = b_dim1 + 1;
        if (b[i__1].r == 0.f && b[i__1].i == 0.f)
        {
            i__2 = b_dim1 + 1;
            b[i__2].r = *eps3;
            b[i__2].i = 0.f; // , expr subst
        }
        *(unsigned char *)trans = 'C';
    }
    *(unsigned char *)normin = 'N';
    i__1 = *n;
    for (its = 1;
            its <= i__1;
            ++its)
    {
        /* Solve U*x = scale*v for a right eigenvector */
        /* or U**H *x = scale*v for a left eigenvector, */
        /* overwriting x on v. */
        clatrs_("Upper", trans, "Nonunit", normin, n, &b[b_offset], ldb, &v[1] , &scale, &rwork[1], &ierr);
        *(unsigned char *)normin = 'Y';
        /* Test for sufficient growth in the norm of v. */
        vnorm = scasum_(n, &v[1], &c__1);
        if (vnorm >= growto * scale)
        {
            goto L120;
        }
        /* Choose new orthogonal starting vector and try again. */
        rtemp = *eps3 / (rootn + 1.f);
        v[1].r = *eps3;
        v[1].i = 0.f; // , expr subst
        i__2 = *n;
        for (i__ = 2;
                i__ <= i__2;
                ++i__)
        {
            i__3 = i__;
            v[i__3].r = rtemp;
            v[i__3].i = 0.f; // , expr subst
            /* L100: */
        }
        i__2 = *n - its + 1;
        i__3 = *n - its + 1;
        r__1 = *eps3 * rootn;
        q__1.r = v[i__3].r - r__1;
        q__1.i = v[i__3].i; // , expr subst
        v[i__2].r = q__1.r;
        v[i__2].i = q__1.i; // , expr subst
        /* L110: */
    }
    /* Failure to find eigenvector in N iterations. */
    *info = 1;
L120: /* Normalize eigenvector. */
    i__ = icamax_(n, &v[1], &c__1);
    i__1 = i__;
    r__3 = 1.f / ((r__1 = v[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&v[i__]), f2c_abs(r__2)));
    csscal_(n, &r__3, &v[1], &c__1);
    return 0;
    /* End of CLAEIN */
}
/* claein_ */
