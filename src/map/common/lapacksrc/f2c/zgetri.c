/* ../netlib/zgetri.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b2 =
{
    1.,0.
}
;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
/* > \brief \b ZGETRI */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGETRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetri. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetri. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetri. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGETRI computes the inverse of a matrix using the LU factorization */
/* > computed by ZGETRF. */
/* > */
/* > This method inverts U and then computes inv(A) by solving the system */
/* > inv(A)*L = inv(U) for inv(A). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the factors L and U from the factorization */
/* > A = P*L*U as computed by ZGETRF. */
/* > On exit, if INFO = 0, the inverse of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices from ZGETRF;
for 1<=i<=N, row i of the */
/* > matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO=0, then WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= max(1,N). */
/* > For optimal performance LWORK >= N*NB, where NB is */
/* > the optimal blocksize returned by ILAENV. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, U(i,i) is exactly zero;
the matrix is */
/* > singular and its inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complex16GEcomputational */
/* ===================================================================== */
/* Subroutine */
int zgetri_(integer *n, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;
    /* Local variables */
    integer i__, j, jb, nb, jj, jp, nn, iws, nbmin;
    extern /* Subroutine */
    int zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *), zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *), zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer ldwork, lwkopt;
    logical lquery;
    extern /* Subroutine */
    int ztrtri_(char *, char *, integer *, doublecomplex *, integer *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;
    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "ZGETRI", " ", n, &c_n1, &c_n1, &c_n1);
    lwkopt = *n * nb;
    work[1].r = (doublereal) lwkopt;
    work[1].i = 0.; // , expr subst
    lquery = *lwork == -1;
    if (*n < 0)
    {
        *info = -1;
    }
    else if (*lda < max(1,*n))
    {
        *info = -3;
    }
    else if (*lwork < max(1,*n) && ! lquery)
    {
        *info = -6;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGETRI", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Form inv(U). If INFO > 0 from ZTRTRI, then U is singular, */
    /* and the inverse is not computed. */
    ztrtri_("Upper", "Non-unit", n, &a[a_offset], lda, info);
    if (*info > 0)
    {
        return 0;
    }
    nbmin = 2;
    ldwork = *n;
    if (nb > 1 && nb < *n)
    {
        /* Computing MAX */
        i__1 = ldwork * nb;
        iws = max(i__1,1);
        if (*lwork < iws)
        {
            nb = *lwork / ldwork;
            /* Computing MAX */
            i__1 = 2;
            i__2 = ilaenv_(&c__2, "ZGETRI", " ", n, &c_n1, &c_n1, & c_n1); // , expr subst
            nbmin = max(i__1,i__2);
        }
    }
    else
    {
        iws = *n;
    }
    /* Solve the equation inv(A)*L = inv(U) for inv(A). */
    if (nb < nbmin || nb >= *n)
    {
        /* Use unblocked code. */
        for (j = *n;
                j >= 1;
                --j)
        {
            /* Copy current column of L to WORK and replace with zeros. */
            i__1 = *n;
            for (i__ = j + 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = i__;
                i__3 = i__ + j * a_dim1;
                work[i__2].r = a[i__3].r;
                work[i__2].i = a[i__3].i; // , expr subst
                i__2 = i__ + j * a_dim1;
                a[i__2].r = 0.;
                a[i__2].i = 0.; // , expr subst
                /* L10: */
            }
            /* Compute current column of inv(A). */
            if (j < *n)
            {
                i__1 = *n - j;
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemv_("No transpose", n, &i__1, &z__1, &a[(j + 1) * a_dim1 + 1], lda, &work[j + 1], &c__1, &c_b2, &a[j * a_dim1 + 1], &c__1);
            }
            /* L20: */
        }
    }
    else
    {
        /* Use blocked code. */
        nn = (*n - 1) / nb * nb + 1;
        i__1 = -nb;
        for (j = nn;
                i__1 < 0 ? j >= 1 : j <= 1;
                j += i__1)
        {
            /* Computing MIN */
            i__2 = nb;
            i__3 = *n - j + 1; // , expr subst
            jb = min(i__2,i__3);
            /* Copy current block column of L to WORK and replace with */
            /* zeros. */
            i__2 = j + jb - 1;
            for (jj = j;
                    jj <= i__2;
                    ++jj)
            {
                i__3 = *n;
                for (i__ = jj + 1;
                        i__ <= i__3;
                        ++i__)
                {
                    i__4 = i__ + (jj - j) * ldwork;
                    i__5 = i__ + jj * a_dim1;
                    work[i__4].r = a[i__5].r;
                    work[i__4].i = a[i__5].i; // , expr subst
                    i__4 = i__ + jj * a_dim1;
                    a[i__4].r = 0.;
                    a[i__4].i = 0.; // , expr subst
                    /* L30: */
                }
                /* L40: */
            }
            /* Compute current block column of inv(A). */
            if (j + jb <= *n)
            {
                i__2 = *n - j - jb + 1;
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemm_("No transpose", "No transpose", n, &jb, &i__2, &z__1, & a[(j + jb) * a_dim1 + 1], lda, &work[j + jb], &ldwork, &c_b2, &a[j * a_dim1 + 1], lda);
            }
            ztrsm_("Right", "Lower", "No transpose", "Unit", n, &jb, &c_b2, & work[j], &ldwork, &a[j * a_dim1 + 1], lda);
            /* L50: */
        }
    }
    /* Apply column interchanges. */
    for (j = *n - 1;
            j >= 1;
            --j)
    {
        jp = ipiv[j];
        if (jp != j)
        {
            zswap_(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
        }
        /* L60: */
    }
    work[1].r = (doublereal) iws;
    work[1].i = 0.; // , expr subst
    return 0;
    /* End of ZGETRI */
}
/* zgetri_ */
