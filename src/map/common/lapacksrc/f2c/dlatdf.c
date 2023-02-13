/* ../netlib/dlatdf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b23 = 1.;
static doublereal c_b37 = -1.;
/* > \brief \b DLATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contrib ution to the reciprocal Dif-estimate. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLATDF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatdf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatdf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatdf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, */
/* JPIV ) */
/* .. Scalar Arguments .. */
/* INTEGER IJOB, LDZ, N */
/* DOUBLE PRECISION RDSCAL, RDSUM */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), JPIV( * ) */
/* DOUBLE PRECISION RHS( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLATDF uses the LU factorization of the n-by-n matrix Z computed by */
/* > DGETC2 and computes a contribution to the reciprocal Dif-estimate */
/* > by solving Z * x = b for x, and choosing the r.h.s. b such that */
/* > the norm of x is as large as possible. On entry RHS = b holds the */
/* > contribution from earlier solved sub-systems, and on return RHS = x. */
/* > */
/* > The factorization of Z returned by DGETC2 has the form Z = P*L*U*Q, */
/* > where P and Q are permutation matrices. L is lower triangular with */
/* > unit diagonal elements and U is upper triangular. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] IJOB */
/* > \verbatim */
/* > IJOB is INTEGER */
/* > IJOB = 2: First compute an approximative null-vector e */
/* > of Z using DGECON, e is normalized and solve for */
/* > Zx = +-e - f with the sign giving the greater value */
/* > of 2-norm(x). About 5 times as expensive as Default. */
/* > IJOB .ne. 2: Local look ahead strategy where all entries of */
/* > the r.h.s. b is choosen as either +1 or -1 (Default). */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix Z. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* > On entry, the LU part of the factorization of the n-by-n */
/* > matrix Z computed by DGETC2: Z = P * L * U * Q */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDA >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] RHS */
/* > \verbatim */
/* > RHS is DOUBLE PRECISION array, dimension (N) */
/* > On entry, RHS contains contributions from other subsystems. */
/* > On exit, RHS contains the solution of the subsystem with */
/* > entries acoording to the value of IJOB (see above). */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSUM */
/* > \verbatim */
/* > RDSUM is DOUBLE PRECISION */
/* > On entry, the sum of squares of computed contributions to */
/* > the Dif-estimate under computation by DTGSYL, where the */
/* > scaling factor RDSCAL (see below) has been factored out. */
/* > On exit, the corresponding sum of squares updated with the */
/* > contributions from the current sub-system. */
/* > If TRANS = 'T' RDSUM is not touched. */
/* > NOTE: RDSUM only makes sense when DTGSY2 is called by STGSYL. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSCAL */
/* > \verbatim */
/* > RDSCAL is DOUBLE PRECISION */
/* > On entry, scaling factor used to prevent overflow in RDSUM. */
/* > On exit, RDSCAL is updated w.r.t. the current contributions */
/* > in RDSUM. */
/* > If TRANS = 'T', RDSCAL is not touched. */
/* > NOTE: RDSCAL only makes sense when DTGSY2 is called by */
/* > DTGSYL. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N). */
/* > The pivot indices;
for 1 <= i <= N, row i of the */
/* > matrix has been interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] JPIV */
/* > \verbatim */
/* > JPIV is INTEGER array, dimension (N). */
/* > The pivot indices;
for 1 <= j <= N, column j of the */
/* > matrix has been interchanged with column JPIV(j). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > This routine is a further developed implementation of algorithm */
/* > BSOLVE in [1] using complete pivoting in the LU factorization. */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* > Umea University, S-901 87 Umea, Sweden. */
/* > \par References: */
/* ================ */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > [1] Bo Kagstrom and Lars Westin, */
/* > Generalized Schur Methods with Condition Estimators for */
/* > Solving the Generalized Sylvester Equation, IEEE Transactions */
/* > on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751. */
/* > */
/* > [2] Peter Poromaa, */
/* > On Efficient and Robust Estimators for the Separation */
/* > between two Regular Matrix Pairs with Applications in */
/* > Condition Estimation. Report IMINF-95.05, Departement of */
/* > Computing Science, Umea University, S-901 87 Umea, Sweden, 1995. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int dlatdf_(integer *ijob, integer *n, doublereal *z__, integer *ldz, doublereal *rhs, doublereal *rdsum, doublereal *rdscal, integer *ipiv, integer *jpiv)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k;
    doublereal bm, bp, xm[8], xp[8];
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer info;
    doublereal temp, work[32];
    extern /* Subroutine */
    int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    doublereal pmone;
    extern /* Subroutine */
    int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *), daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    doublereal sminu;
    integer iwork[8];
    doublereal splus;
    extern /* Subroutine */
    int dgesc2_(integer *, doublereal *, integer *, doublereal *, integer *, integer *, doublereal *), dgecon_(char *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, integer *, integer *), dlassq_(integer *, doublereal *, integer *, doublereal *, doublereal *), dlaswp_( integer *, doublereal *, integer *, integer *, integer *, integer *, integer *);
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --rhs;
    --ipiv;
    --jpiv;
    /* Function Body */
    if (*ijob != 2)
    {
        /* Apply permutations IPIV to RHS */
        i__1 = *n - 1;
        dlaswp_(&c__1, &rhs[1], ldz, &c__1, &i__1, &ipiv[1], &c__1);
        /* Solve for L-part choosing RHS either to +1 or -1. */
        pmone = -1.;
        i__1 = *n - 1;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            bp = rhs[j] + 1.;
            bm = rhs[j] - 1.;
            splus = 1.;
            /* Look-ahead for L-part RHS(1:N-1) = + or -1, SPLUS and */
            /* SMIN computed more efficiently than in BSOLVE [1]. */
            i__2 = *n - j;
            splus += ddot_(&i__2, &z__[j + 1 + j * z_dim1], &c__1, &z__[j + 1 + j * z_dim1], &c__1);
            i__2 = *n - j;
            sminu = ddot_(&i__2, &z__[j + 1 + j * z_dim1], &c__1, &rhs[j + 1], &c__1);
            splus *= rhs[j];
            if (splus > sminu)
            {
                rhs[j] = bp;
            }
            else if (sminu > splus)
            {
                rhs[j] = bm;
            }
            else
            {
                /* In this case the updating sums are equal and we can */
                /* choose RHS(J) +1 or -1. The first time this happens */
                /* we choose -1, thereafter +1. This is a simple way to */
                /* get good estimates of matrices like Byers well-known */
                /* example (see [1]). (Not done in BSOLVE.) */
                rhs[j] += pmone;
                pmone = 1.;
            }
            /* Compute the remaining r.h.s. */
            temp = -rhs[j];
            i__2 = *n - j;
            daxpy_(&i__2, &temp, &z__[j + 1 + j * z_dim1], &c__1, &rhs[j + 1], &c__1);
            /* L10: */
        }
        /* Solve for U-part, look-ahead for RHS(N) = +-1. This is not done */
        /* in BSOLVE and will hopefully give us a better estimate because */
        /* any ill-conditioning of the original matrix is transfered to U */
        /* and not to L. U(N, N) is an approximation to sigma_min(LU). */
        i__1 = *n - 1;
        dcopy_(&i__1, &rhs[1], &c__1, xp, &c__1);
        xp[*n - 1] = rhs[*n] + 1.;
        rhs[*n] += -1.;
        splus = 0.;
        sminu = 0.;
        for (i__ = *n;
                i__ >= 1;
                --i__)
        {
            temp = 1. / z__[i__ + i__ * z_dim1];
            xp[i__ - 1] *= temp;
            rhs[i__] *= temp;
            i__1 = *n;
            for (k = i__ + 1;
                    k <= i__1;
                    ++k)
            {
                xp[i__ - 1] -= xp[k - 1] * (z__[i__ + k * z_dim1] * temp);
                rhs[i__] -= rhs[k] * (z__[i__ + k * z_dim1] * temp);
                /* L20: */
            }
            splus += (d__1 = xp[i__ - 1], f2c_abs(d__1));
            sminu += (d__1 = rhs[i__], f2c_abs(d__1));
            /* L30: */
        }
        if (splus > sminu)
        {
            dcopy_(n, xp, &c__1, &rhs[1], &c__1);
        }
        /* Apply the permutations JPIV to the computed solution (RHS) */
        i__1 = *n - 1;
        dlaswp_(&c__1, &rhs[1], ldz, &c__1, &i__1, &jpiv[1], &c_n1);
        /* Compute the sum of squares */
        dlassq_(n, &rhs[1], &c__1, rdscal, rdsum);
    }
    else
    {
        /* IJOB = 2, Compute approximate nullvector XM of Z */
        dgecon_("I", n, &z__[z_offset], ldz, &c_b23, &temp, work, iwork, & info);
        dcopy_(n, &work[*n], &c__1, xm, &c__1);
        /* Compute RHS */
        i__1 = *n - 1;
        dlaswp_(&c__1, xm, ldz, &c__1, &i__1, &ipiv[1], &c_n1);
        temp = 1. / sqrt(ddot_(n, xm, &c__1, xm, &c__1));
        dscal_(n, &temp, xm, &c__1);
        dcopy_(n, xm, &c__1, xp, &c__1);
        daxpy_(n, &c_b23, &rhs[1], &c__1, xp, &c__1);
        daxpy_(n, &c_b37, xm, &c__1, &rhs[1], &c__1);
        dgesc2_(n, &z__[z_offset], ldz, &rhs[1], &ipiv[1], &jpiv[1], &temp);
        dgesc2_(n, &z__[z_offset], ldz, xp, &ipiv[1], &jpiv[1], &temp);
        if (dasum_(n, xp, &c__1) > dasum_(n, &rhs[1], &c__1))
        {
            dcopy_(n, xp, &c__1, &rhs[1], &c__1);
        }
        /* Compute the sum of squares */
        dlassq_(n, &rhs[1], &c__1, rdscal, rdsum);
    }
    return 0;
    /* End of DLATDF */
}
/* dlatdf_ */
