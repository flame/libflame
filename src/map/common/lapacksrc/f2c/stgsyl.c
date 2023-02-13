/* ../netlib/stgsyl.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__5 = 5;
static real c_b14 = 0.f;
static integer c__1 = 1;
static real c_b51 = -1.f;
static real c_b52 = 1.f;
/* > \brief \b STGSYL */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STGSYL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsyl. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsyl. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsyl. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, */
/* LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK, */
/* IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, */
/* $ LWORK, M, N */
/* REAL DIF, SCALE */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL A( LDA, * ), B( LDB, * ), C( LDC, * ), */
/* $ D( LDD, * ), E( LDE, * ), F( LDF, * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGSYL solves the generalized Sylvester equation: */
/* > */
/* > A * R - L * B = scale * C (1) */
/* > D * R - L * E = scale * F */
/* > */
/* > where R and L are unknown m-by-n matrices, (A, D), (B, E) and */
/* > (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n, */
/* > respectively, with real entries. (A, D) and (B, E) must be in */
/* > generalized (real) Schur canonical form, i.e. A, B are upper quasi */
/* > triangular and D, E are upper triangular. */
/* > */
/* > The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output */
/* > scaling factor chosen to avoid overflow. */
/* > */
/* > In matrix notation (1) is equivalent to solve Zx = scale b, where */
/* > Z is defined as */
/* > */
/* > Z = [ kron(In, A) -kron(B**T, Im) ] (2) */
/* > [ kron(In, D) -kron(E**T, Im) ]. */
/* > */
/* > Here Ik is the identity matrix of size k and X**T is the transpose of */
/* > X. kron(X, Y) is the Kronecker product between the matrices X and Y. */
/* > */
/* > If TRANS = 'T', STGSYL solves the transposed system Z**T*y = scale*b, */
/* > which is equivalent to solve for R and L in */
/* > */
/* > A**T * R + D**T * L = scale * C (3) */
/* > R * B**T + L * E**T = scale * -F */
/* > */
/* > This case (TRANS = 'T') is used to compute an one-norm-based estimate */
/* > of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D) */
/* > and (B,E), using SLACON. */
/* > */
/* > If IJOB >= 1, STGSYL computes a Frobenius norm-based estimate */
/* > of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the */
/* > reciprocal of the smallest singular value of Z. See [1-2] for more */
/* > information. */
/* > */
/* > This is a level 3 BLAS algorithm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N', solve the generalized Sylvester equation (1). */
/* > = 'T', solve the 'transposed' system (3). */
/* > \endverbatim */
/* > */
/* > \param[in] IJOB */
/* > \verbatim */
/* > IJOB is INTEGER */
/* > Specifies what kind of functionality to be performed. */
/* > =0: solve (1) only. */
/* > =1: The functionality of 0 and 3. */
/* > =2: The functionality of 0 and 4. */
/* > =3: Only an estimate of Dif[(A,D), (B,E)] is computed. */
/* > (look ahead strategy IJOB = 1 is used). */
/* > =4: Only an estimate of Dif[(A,D), (B,E)] is computed. */
/* > ( SGECON on sub-systems is used ). */
/* > Not referenced if TRANS = 'T'. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The order of the matrices A and D, and the row dimension of */
/* > the matrices C, F, R and L. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices B and E, and the column dimension */
/* > of the matrices C, F, R and L. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, M) */
/* > The upper quasi triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB, N) */
/* > The upper quasi triangular matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is REAL array, dimension (LDC, N) */
/* > On entry, C contains the right-hand-side of the first matrix */
/* > equation in (1) or (3). */
/* > On exit, if IJOB = 0, 1 or 2, C has been overwritten by */
/* > the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R, */
/* > the solution achieved during the computation of the */
/* > Dif-estimate. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (LDD, M) */
/* > The upper triangular matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] LDD */
/* > \verbatim */
/* > LDD is INTEGER */
/* > The leading dimension of the array D. LDD >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is REAL array, dimension (LDE, N) */
/* > The upper triangular matrix E. */
/* > \endverbatim */
/* > */
/* > \param[in] LDE */
/* > \verbatim */
/* > LDE is INTEGER */
/* > The leading dimension of the array E. LDE >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] F */
/* > \verbatim */
/* > F is REAL array, dimension (LDF, N) */
/* > On entry, F contains the right-hand-side of the second matrix */
/* > equation in (1) or (3). */
/* > On exit, if IJOB = 0, 1 or 2, F has been overwritten by */
/* > the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L, */
/* > the solution achieved during the computation of the */
/* > Dif-estimate. */
/* > \endverbatim */
/* > */
/* > \param[in] LDF */
/* > \verbatim */
/* > LDF is INTEGER */
/* > The leading dimension of the array F. LDF >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[out] DIF */
/* > \verbatim */
/* > DIF is REAL */
/* > On exit DIF is the reciprocal of a lower bound of the */
/* > reciprocal of the Dif-function, i.e. DIF is an upper bound of */
/* > Dif[(A,D), (B,E)] = sigma_min(Z), where Z as in (2). */
/* > IF IJOB = 0 or TRANS = 'T', DIF is not touched. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is REAL */
/* > On exit SCALE is the scaling factor in (1) or (3). */
/* > If 0 < SCALE < 1, C and F hold the solutions R and L, resp., */
/* > to a slightly perturbed system but the input matrices A, B, D */
/* > and E have not been changed. If SCALE = 0, C and F hold the */
/* > solutions R and L, respectively, to the homogeneous system */
/* > with C = F = 0. Normally, SCALE = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK > = 1. */
/* > If IJOB = 1 or 2 and TRANS = 'N', LWORK >= max(1,2*M*N). */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (M+N+6) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > =0: successful exit */
/* > <0: If INFO = -i, the i-th argument had an illegal value. */
/* > >0: (A, D) and (B, E) have common or close eigenvalues. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realSYcomputational */
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
/* > [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software */
/* > for Solving the Generalized Sylvester Equation and Estimating the */
/* > Separation between Regular Matrix Pairs, Report UMINF - 93.23, */
/* > Department of Computing Science, Umea University, S-901 87 Umea, */
/* > Sweden, December 1993, Revised April 1994, Also as LAPACK Working */
/* > Note 75. To appear in ACM Trans. on Math. Software, Vol 22, */
/* > No 1, 1996. */
/* > */
/* > [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester */
/* > Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal. */
/* > Appl., 15(4):1045-1060, 1994 */
/* > */
/* > [3] B. Kagstrom and L. Westin, Generalized Schur Methods with */
/* > Condition Estimators for Solving the Generalized Sylvester */
/* > Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7, */
/* > July 1989, pp 745-751. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int stgsyl_(char *trans, integer *ijob, integer *m, integer * n, real *a, integer *lda, real *b, integer *ldb, real *c__, integer * ldc, real *d__, integer *ldd, real *e, integer *lde, real *f, integer *ldf, real *scale, real *dif, real *work, integer *lwork, integer * iwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, d_offset, e_dim1, e_offset, f_dim1, f_offset, i__1, i__2, i__3, i__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k, p, q, ie, je, mb, nb, is, js, pq;
    real dsum;
    integer ppqq;
    extern logical lsame_(char *, char *);
    integer ifunc;
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    integer linfo;
    extern /* Subroutine */
    int sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    integer lwmin;
    real scale2, dscale;
    extern /* Subroutine */
    int stgsy2_(char *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, real * , integer *, real *, integer *, real *, integer *, real *, real *, real *, integer *, integer *, integer *);
    real scaloc;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *), slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
    integer iround;
    logical notran;
    integer isolve;
    logical lquery;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* Replaced various illegal calls to SCOPY by calls to SLASET. */
    /* Sven Hammarling, 1/5/02. */
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
    /* Decode and test input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N");
    lquery = *lwork == -1;
    if (! notran && ! lsame_(trans, "T"))
    {
        *info = -1;
    }
    else if (notran)
    {
        if (*ijob < 0 || *ijob > 4)
        {
            *info = -2;
        }
    }
    if (*info == 0)
    {
        if (*m <= 0)
        {
            *info = -3;
        }
        else if (*n <= 0)
        {
            *info = -4;
        }
        else if (*lda < max(1,*m))
        {
            *info = -6;
        }
        else if (*ldb < max(1,*n))
        {
            *info = -8;
        }
        else if (*ldc < max(1,*m))
        {
            *info = -10;
        }
        else if (*ldd < max(1,*m))
        {
            *info = -12;
        }
        else if (*lde < max(1,*n))
        {
            *info = -14;
        }
        else if (*ldf < max(1,*m))
        {
            *info = -16;
        }
    }
    if (*info == 0)
    {
        if (notran)
        {
            if (*ijob == 1 || *ijob == 2)
            {
                /* Computing MAX */
                i__1 = 1;
                i__2 = (*m << 1) * *n; // , expr subst
                lwmin = max(i__1,i__2);
            }
            else
            {
                lwmin = 1;
            }
        }
        else
        {
            lwmin = 1;
        }
        work[1] = (real) lwmin;
        if (*lwork < lwmin && ! lquery)
        {
            *info = -20;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("STGSYL", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        *scale = 1.f;
        if (notran)
        {
            if (*ijob != 0)
            {
                *dif = 0.f;
            }
        }
        return 0;
    }
    /* Determine optimal block sizes MB and NB */
    mb = ilaenv_(&c__2, "STGSYL", trans, m, n, &c_n1, &c_n1);
    nb = ilaenv_(&c__5, "STGSYL", trans, m, n, &c_n1, &c_n1);
    isolve = 1;
    ifunc = 0;
    if (notran)
    {
        if (*ijob >= 3)
        {
            ifunc = *ijob - 2;
            slaset_("F", m, n, &c_b14, &c_b14, &c__[c_offset], ldc) ;
            slaset_("F", m, n, &c_b14, &c_b14, &f[f_offset], ldf);
        }
        else if (*ijob >= 1 && notran)
        {
            isolve = 2;
        }
    }
    if (mb <= 1 && nb <= 1 || mb >= *m && nb >= *n)
    {
        i__1 = isolve;
        for (iround = 1;
                iround <= i__1;
                ++iround)
        {
            /* Use unblocked Level 2 solver */
            dscale = 0.f;
            dsum = 1.f;
            pq = 0;
            stgsy2_(trans, &ifunc, m, n, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, &e[e_offset], lde, &f[f_offset], ldf, scale, &dsum, &dscale, &iwork[1], &pq, info);
            if (dscale != 0.f)
            {
                if (*ijob == 1 || *ijob == 3)
                {
                    *dif = sqrt((real) ((*m << 1) * *n)) / (dscale * sqrt( dsum));
                }
                else
                {
                    *dif = sqrt((real) pq) / (dscale * sqrt(dsum));
                }
            }
            if (isolve == 2 && iround == 1)
            {
                if (notran)
                {
                    ifunc = *ijob;
                }
                scale2 = *scale;
                slacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m);
                slacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m);
                slaset_("F", m, n, &c_b14, &c_b14, &c__[c_offset], ldc);
                slaset_("F", m, n, &c_b14, &c_b14, &f[f_offset], ldf);
            }
            else if (isolve == 2 && iround == 2)
            {
                slacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc);
                slacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf);
                *scale = scale2;
            }
            /* L30: */
        }
        return 0;
    }
    /* Determine block structure of A */
    p = 0;
    i__ = 1;
L40:
    if (i__ > *m)
    {
        goto L50;
    }
    ++p;
    iwork[p] = i__;
    i__ += mb;
    if (i__ >= *m)
    {
        goto L50;
    }
    if (a[i__ + (i__ - 1) * a_dim1] != 0.f)
    {
        ++i__;
    }
    goto L40;
L50:
    iwork[p + 1] = *m + 1;
    if (iwork[p] == iwork[p + 1])
    {
        --p;
    }
    /* Determine block structure of B */
    q = p + 1;
    j = 1;
L60:
    if (j > *n)
    {
        goto L70;
    }
    ++q;
    iwork[q] = j;
    j += nb;
    if (j >= *n)
    {
        goto L70;
    }
    if (b[j + (j - 1) * b_dim1] != 0.f)
    {
        ++j;
    }
    goto L60;
L70:
    iwork[q + 1] = *n + 1;
    if (iwork[q] == iwork[q + 1])
    {
        --q;
    }
    if (notran)
    {
        i__1 = isolve;
        for (iround = 1;
                iround <= i__1;
                ++iround)
        {
            /* Solve (I, J)-subsystem */
            /* A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J) */
            /* D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J) */
            /* for I = P, P - 1,..., 1;
            J = 1, 2,..., Q */
            dscale = 0.f;
            dsum = 1.f;
            pq = 0;
            *scale = 1.f;
            i__2 = q;
            for (j = p + 2;
                    j <= i__2;
                    ++j)
            {
                js = iwork[j];
                je = iwork[j + 1] - 1;
                nb = je - js + 1;
                for (i__ = p;
                        i__ >= 1;
                        --i__)
                {
                    is = iwork[i__];
                    ie = iwork[i__ + 1] - 1;
                    mb = ie - is + 1;
                    ppqq = 0;
                    stgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], lda, &b[js + js * b_dim1], ldb, &c__[is + js * c_dim1], ldc, &d__[is + is * d_dim1], ldd, &e[js + js * e_dim1], lde, &f[is + js * f_dim1], ldf, & scaloc, &dsum, &dscale, &iwork[q + 2], &ppqq, & linfo);
                    if (linfo > 0)
                    {
                        *info = linfo;
                    }
                    pq += ppqq;
                    if (scaloc != 1.f)
                    {
                        i__3 = js - 1;
                        for (k = 1;
                                k <= i__3;
                                ++k)
                        {
                            sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
                            sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                            /* L80: */
                        }
                        i__3 = je;
                        for (k = js;
                                k <= i__3;
                                ++k)
                        {
                            i__4 = is - 1;
                            sscal_(&i__4, &scaloc, &c__[k * c_dim1 + 1], & c__1);
                            i__4 = is - 1;
                            sscal_(&i__4, &scaloc, &f[k * f_dim1 + 1], &c__1);
                            /* L90: */
                        }
                        i__3 = je;
                        for (k = js;
                                k <= i__3;
                                ++k)
                        {
                            i__4 = *m - ie;
                            sscal_(&i__4, &scaloc, &c__[ie + 1 + k * c_dim1], &c__1);
                            i__4 = *m - ie;
                            sscal_(&i__4, &scaloc, &f[ie + 1 + k * f_dim1], & c__1);
                            /* L100: */
                        }
                        i__3 = *n;
                        for (k = je + 1;
                                k <= i__3;
                                ++k)
                        {
                            sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
                            sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                            /* L110: */
                        }
                        *scale *= scaloc;
                    }
                    /* Substitute R(I, J) and L(I, J) into remaining */
                    /* equation. */
                    if (i__ > 1)
                    {
                        i__3 = is - 1;
                        sgemm_("N", "N", &i__3, &nb, &mb, &c_b51, &a[is * a_dim1 + 1], lda, &c__[is + js * c_dim1], ldc, &c_b52, &c__[js * c_dim1 + 1], ldc);
                        i__3 = is - 1;
                        sgemm_("N", "N", &i__3, &nb, &mb, &c_b51, &d__[is * d_dim1 + 1], ldd, &c__[is + js * c_dim1], ldc, &c_b52, &f[js * f_dim1 + 1], ldf);
                    }
                    if (j < q)
                    {
                        i__3 = *n - je;
                        sgemm_("N", "N", &mb, &i__3, &nb, &c_b52, &f[is + js * f_dim1], ldf, &b[js + (je + 1) * b_dim1], ldb, &c_b52, &c__[is + (je + 1) * c_dim1], ldc);
                        i__3 = *n - je;
                        sgemm_("N", "N", &mb, &i__3, &nb, &c_b52, &f[is + js * f_dim1], ldf, &e[js + (je + 1) * e_dim1], lde, &c_b52, &f[is + (je + 1) * f_dim1], ldf);
                    }
                    /* L120: */
                }
                /* L130: */
            }
            if (dscale != 0.f)
            {
                if (*ijob == 1 || *ijob == 3)
                {
                    *dif = sqrt((real) ((*m << 1) * *n)) / (dscale * sqrt( dsum));
                }
                else
                {
                    *dif = sqrt((real) pq) / (dscale * sqrt(dsum));
                }
            }
            if (isolve == 2 && iround == 1)
            {
                if (notran)
                {
                    ifunc = *ijob;
                }
                scale2 = *scale;
                slacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m);
                slacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m);
                slaset_("F", m, n, &c_b14, &c_b14, &c__[c_offset], ldc);
                slaset_("F", m, n, &c_b14, &c_b14, &f[f_offset], ldf);
            }
            else if (isolve == 2 && iround == 2)
            {
                slacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc);
                slacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf);
                *scale = scale2;
            }
            /* L150: */
        }
    }
    else
    {
        /* Solve transposed (I, J)-subsystem */
        /* A(I, I)**T * R(I, J) + D(I, I)**T * L(I, J) = C(I, J) */
        /* R(I, J) * B(J, J)**T + L(I, J) * E(J, J)**T = -F(I, J) */
        /* for I = 1,2,..., P;
        J = Q, Q-1,..., 1 */
        *scale = 1.f;
        i__1 = p;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            is = iwork[i__];
            ie = iwork[i__ + 1] - 1;
            mb = ie - is + 1;
            i__2 = p + 2;
            for (j = q;
                    j >= i__2;
                    --j)
            {
                js = iwork[j];
                je = iwork[j + 1] - 1;
                nb = je - js + 1;
                stgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], lda, & b[js + js * b_dim1], ldb, &c__[is + js * c_dim1], ldc, &d__[is + is * d_dim1], ldd, &e[js + js * e_dim1], lde, &f[is + js * f_dim1], ldf, &scaloc, &dsum, & dscale, &iwork[q + 2], &ppqq, &linfo);
                if (linfo > 0)
                {
                    *info = linfo;
                }
                if (scaloc != 1.f)
                {
                    i__3 = js - 1;
                    for (k = 1;
                            k <= i__3;
                            ++k)
                    {
                        sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
                        sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                        /* L160: */
                    }
                    i__3 = je;
                    for (k = js;
                            k <= i__3;
                            ++k)
                    {
                        i__4 = is - 1;
                        sscal_(&i__4, &scaloc, &c__[k * c_dim1 + 1], &c__1);
                        i__4 = is - 1;
                        sscal_(&i__4, &scaloc, &f[k * f_dim1 + 1], &c__1);
                        /* L170: */
                    }
                    i__3 = je;
                    for (k = js;
                            k <= i__3;
                            ++k)
                    {
                        i__4 = *m - ie;
                        sscal_(&i__4, &scaloc, &c__[ie + 1 + k * c_dim1], & c__1);
                        i__4 = *m - ie;
                        sscal_(&i__4, &scaloc, &f[ie + 1 + k * f_dim1], &c__1) ;
                        /* L180: */
                    }
                    i__3 = *n;
                    for (k = je + 1;
                            k <= i__3;
                            ++k)
                    {
                        sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
                        sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                        /* L190: */
                    }
                    *scale *= scaloc;
                }
                /* Substitute R(I, J) and L(I, J) into remaining equation. */
                if (j > p + 2)
                {
                    i__3 = js - 1;
                    sgemm_("N", "T", &mb, &i__3, &nb, &c_b52, &c__[is + js * c_dim1], ldc, &b[js * b_dim1 + 1], ldb, &c_b52, & f[is + f_dim1], ldf);
                    i__3 = js - 1;
                    sgemm_("N", "T", &mb, &i__3, &nb, &c_b52, &f[is + js * f_dim1], ldf, &e[js * e_dim1 + 1], lde, &c_b52, & f[is + f_dim1], ldf);
                }
                if (i__ < p)
                {
                    i__3 = *m - ie;
                    sgemm_("T", "N", &i__3, &nb, &mb, &c_b51, &a[is + (ie + 1) * a_dim1], lda, &c__[is + js * c_dim1], ldc, & c_b52, &c__[ie + 1 + js * c_dim1], ldc);
                    i__3 = *m - ie;
                    sgemm_("T", "N", &i__3, &nb, &mb, &c_b51, &d__[is + (ie + 1) * d_dim1], ldd, &f[is + js * f_dim1], ldf, & c_b52, &c__[ie + 1 + js * c_dim1], ldc);
                }
                /* L200: */
            }
            /* L210: */
        }
    }
    work[1] = (real) lwmin;
    return 0;
    /* End of STGSYL */
}
/* stgsyl_ */
