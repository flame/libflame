/* ../netlib/stgsy2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__8 = 8;
static integer c__1 = 1;
static real c_b27 = -1.f;
static real c_b42 = 1.f;
static real c_b56 = 0.f;
/* > \brief \b STGSY2 solves the generalized Sylvester equation (unblocked algorithm). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STGSY2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsy2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsy2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsy2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, */
/* LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL, */
/* IWORK, PQ, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N, */
/* $ PQ */
/* REAL RDSCAL, RDSUM, SCALE */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL A( LDA, * ), B( LDB, * ), C( LDC, * ), */
/* $ D( LDD, * ), E( LDE, * ), F( LDF, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGSY2 solves the generalized Sylvester equation: */
/* > */
/* > A * R - L * B = scale * C (1) */
/* > D * R - L * E = scale * F, */
/* > */
/* > using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices, */
/* > (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M, */
/* > N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E) */
/* > must be in generalized Schur canonical form, i.e. A, B are upper */
/* > quasi triangular and D, E are upper triangular. The solution (R, L) */
/* > overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor */
/* > chosen to avoid overflow. */
/* > */
/* > In matrix notation solving equation (1) corresponds to solve */
/* > Z*x = scale*b, where Z is defined as */
/* > */
/* > Z = [ kron(In, A) -kron(B**T, Im) ] (2) */
/* > [ kron(In, D) -kron(E**T, Im) ], */
/* > */
/* > Ik is the identity matrix of size k and X**T is the transpose of X. */
/* > kron(X, Y) is the Kronecker product between the matrices X and Y. */
/* > In the process of solving (1), we solve a number of such systems */
/* > where Dim(In), Dim(In) = 1 or 2. */
/* > */
/* > If TRANS = 'T', solve the transposed system Z**T*y = scale*b for y, */
/* > which is equivalent to solve for R and L in */
/* > */
/* > A**T * R + D**T * L = scale * C (3) */
/* > R * B**T + L * E**T = scale * -F */
/* > */
/* > This case is used to compute an estimate of Dif[(A, D), (B, E)] = */
/* > sigma_min(Z) using reverse communicaton with SLACON. */
/* > */
/* > STGSY2 also (IJOB >= 1) contributes to the computation in STGSYL */
/* > of an upper bound on the separation between to matrix pairs. Then */
/* > the input (A, D), (B, E) are sub-pencils of the matrix pair in */
/* > STGSYL. See STGSYL for details. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N', solve the generalized Sylvester equation (1). */
/* > = 'T': solve the 'transposed' system (3). */
/* > \endverbatim */
/* > */
/* > \param[in] IJOB */
/* > \verbatim */
/* > IJOB is INTEGER */
/* > Specifies what kind of functionality to be performed. */
/* > = 0: solve (1) only. */
/* > = 1: A contribution from this subsystem to a Frobenius */
/* > norm-based estimate of the separation between two matrix */
/* > pairs is computed. (look ahead strategy is used). */
/* > = 2: A contribution from this subsystem to a Frobenius */
/* > norm-based estimate of the separation between two matrix */
/* > pairs is computed. (SGECON on sub-systems is used.) */
/* > Not referenced if TRANS = 'T'. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > On entry, M specifies the order of A and D, and the row */
/* > dimension of C, F, R and L. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > On entry, N specifies the order of B and E, and the column */
/* > dimension of C, F, R and L. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, M) */
/* > On entry, A contains an upper quasi triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the matrix A. LDA >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB, N) */
/* > On entry, B contains an upper quasi triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the matrix B. LDB >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is REAL array, dimension (LDC, N) */
/* > On entry, C contains the right-hand-side of the first matrix */
/* > equation in (1). */
/* > On exit, if IJOB = 0, C has been overwritten by the */
/* > solution R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the matrix C. LDC >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (LDD, M) */
/* > On entry, D contains an upper triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDD */
/* > \verbatim */
/* > LDD is INTEGER */
/* > The leading dimension of the matrix D. LDD >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is REAL array, dimension (LDE, N) */
/* > On entry, E contains an upper triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDE */
/* > \verbatim */
/* > LDE is INTEGER */
/* > The leading dimension of the matrix E. LDE >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] F */
/* > \verbatim */
/* > F is REAL array, dimension (LDF, N) */
/* > On entry, F contains the right-hand-side of the second matrix */
/* > equation in (1). */
/* > On exit, if IJOB = 0, F has been overwritten by the */
/* > solution L. */
/* > \endverbatim */
/* > */
/* > \param[in] LDF */
/* > \verbatim */
/* > LDF is INTEGER */
/* > The leading dimension of the matrix F. LDF >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is REAL */
/* > On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions */
/* > R and L (C and F on entry) will hold the solutions to a */
/* > slightly perturbed system but the input matrices A, B, D and */
/* > E have not been changed. If SCALE = 0, R and L will hold the */
/* > solutions to the homogeneous system with C = F = 0. Normally, */
/* > SCALE = 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSUM */
/* > \verbatim */
/* > RDSUM is REAL */
/* > On entry, the sum of squares of computed contributions to */
/* > the Dif-estimate under computation by STGSYL, where the */
/* > scaling factor RDSCAL (see below) has been factored out. */
/* > On exit, the corresponding sum of squares updated with the */
/* > contributions from the current sub-system. */
/* > If TRANS = 'T' RDSUM is not touched. */
/* > NOTE: RDSUM only makes sense when STGSY2 is called by STGSYL. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSCAL */
/* > \verbatim */
/* > RDSCAL is REAL */
/* > On entry, scaling factor used to prevent overflow in RDSUM. */
/* > On exit, RDSCAL is updated w.r.t. the current contributions */
/* > in RDSUM. */
/* > If TRANS = 'T', RDSCAL is not touched. */
/* > NOTE: RDSCAL only makes sense when STGSY2 is called by */
/* > STGSYL. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (M+N+2) */
/* > \endverbatim */
/* > */
/* > \param[out] PQ */
/* > \verbatim */
/* > PQ is INTEGER */
/* > On exit, the number of subsystems (of size 2-by-2, 4-by-4 and */
/* > 8-by-8) solved by this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > On exit, if INFO is set to */
/* > =0: Successful exit */
/* > <0: If INFO = -i, the i-th argument had an illegal value. */
/* > >0: The matrix pairs (A, D) and (B, E) have common or very */
/* > close eigenvalues. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realSYauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* > Umea University, S-901 87 Umea, Sweden. */
/* ===================================================================== */
/* Subroutine */
int stgsy2_(char *trans, integer *ijob, integer *m, integer * n, real *a, integer *lda, real *b, integer *ldb, real *c__, integer * ldc, real *d__, integer *ldd, real *e, integer *lde, real *f, integer *ldf, real *scale, real *rdsum, real *rdscal, integer *iwork, integer *pq, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, d_offset, e_dim1, e_offset, f_dim1, f_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, j, k, p, q;
    real z__[64] /* was [8][8] */
    ;
    integer ie, je, mb, nb, ii, jj, is, js;
    real rhs[8];
    integer isp1, jsp1;
    extern /* Subroutine */
    int sger_(integer *, integer *, real *, real *, integer *, real *, integer *, real *, integer *);
    integer ierr, zdim, ipiv[8], jpiv[8];
    real alpha;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *), sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *), sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *), scopy_(integer *, real *, integer *, real *, integer *), saxpy_(integer *, real *, real *, integer *, real *, integer *), sgesc2_(integer *, real *, integer *, real *, integer *, integer * , real *), sgetc2_(integer *, real *, integer *, integer *, integer *, integer *);
    real scaloc;
    extern /* Subroutine */
    int slatdf_(integer *, integer *, real *, integer *, real *, real *, real *, integer *, integer *), xerbla_(char *, integer *), slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
    logical notran;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* Replaced various illegal calls to SCOPY by calls to SLASET. */
    /* Sven Hammarling, 27/5/02. */
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
    --iwork;
    /* Function Body */
    *info = 0;
    ierr = 0;
    notran = lsame_(trans, "N");
    if (! notran && ! lsame_(trans, "T"))
    {
        *info = -1;
    }
    else if (notran)
    {
        if (*ijob < 0 || *ijob > 2)
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
            *info = -5;
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
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("STGSY2", &i__1);
        return 0;
    }
    /* Determine block structure of A */
    *pq = 0;
    p = 0;
    i__ = 1;
L10:
    if (i__ > *m)
    {
        goto L20;
    }
    ++p;
    iwork[p] = i__;
    if (i__ == *m)
    {
        goto L20;
    }
    if (a[i__ + 1 + i__ * a_dim1] != 0.f)
    {
        i__ += 2;
    }
    else
    {
        ++i__;
    }
    goto L10;
L20:
    iwork[p + 1] = *m + 1;
    /* Determine block structure of B */
    q = p + 1;
    j = 1;
L30:
    if (j > *n)
    {
        goto L40;
    }
    ++q;
    iwork[q] = j;
    if (j == *n)
    {
        goto L40;
    }
    if (b[j + 1 + j * b_dim1] != 0.f)
    {
        j += 2;
    }
    else
    {
        ++j;
    }
    goto L30;
L40:
    iwork[q + 1] = *n + 1;
    *pq = p * (q - p - 1);
    if (notran)
    {
        /* Solve (I, J) - subsystem */
        /* A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J) */
        /* D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J) */
        /* for I = P, P - 1, ..., 1;
        J = 1, 2, ..., Q */
        *scale = 1.f;
        scaloc = 1.f;
        i__1 = q;
        for (j = p + 2;
                j <= i__1;
                ++j)
        {
            js = iwork[j];
            jsp1 = js + 1;
            je = iwork[j + 1] - 1;
            nb = je - js + 1;
            for (i__ = p;
                    i__ >= 1;
                    --i__)
            {
                is = iwork[i__];
                isp1 = is + 1;
                ie = iwork[i__ + 1] - 1;
                mb = ie - is + 1;
                zdim = mb * nb << 1;
                if (mb == 1 && nb == 1)
                {
                    /* Build a 2-by-2 system Z * x = RHS */
                    z__[0] = a[is + is * a_dim1];
                    z__[1] = d__[is + is * d_dim1];
                    z__[8] = -b[js + js * b_dim1];
                    z__[9] = -e[js + js * e_dim1];
                    /* Set up right hand side(s) */
                    rhs[0] = c__[is + js * c_dim1];
                    rhs[1] = f[is + js * f_dim1];
                    /* Solve Z * x = RHS */
                    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
                    if (ierr > 0)
                    {
                        *info = ierr;
                    }
                    if (*ijob == 0)
                    {
                        sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
                        if (scaloc != 1.f)
                        {
                            i__2 = *n;
                            for (k = 1;
                                    k <= i__2;
                                    ++k)
                            {
                                sscal_(m, &scaloc, &c__[k * c_dim1 + 1], & c__1);
                                sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                                /* L50: */
                            }
                            *scale *= scaloc;
                        }
                    }
                    else
                    {
                        slatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, ipiv, jpiv);
                    }
                    /* Unpack solution vector(s) */
                    c__[is + js * c_dim1] = rhs[0];
                    f[is + js * f_dim1] = rhs[1];
                    /* Substitute R(I, J) and L(I, J) into remaining */
                    /* equation. */
                    if (i__ > 1)
                    {
                        alpha = -rhs[0];
                        i__2 = is - 1;
                        saxpy_(&i__2, &alpha, &a[is * a_dim1 + 1], &c__1, & c__[js * c_dim1 + 1], &c__1);
                        i__2 = is - 1;
                        saxpy_(&i__2, &alpha, &d__[is * d_dim1 + 1], &c__1, & f[js * f_dim1 + 1], &c__1);
                    }
                    if (j < q)
                    {
                        i__2 = *n - je;
                        saxpy_(&i__2, &rhs[1], &b[js + (je + 1) * b_dim1], ldb, &c__[is + (je + 1) * c_dim1], ldc);
                        i__2 = *n - je;
                        saxpy_(&i__2, &rhs[1], &e[js + (je + 1) * e_dim1], lde, &f[is + (je + 1) * f_dim1], ldf);
                    }
                }
                else if (mb == 1 && nb == 2)
                {
                    /* Build a 4-by-4 system Z * x = RHS */
                    z__[0] = a[is + is * a_dim1];
                    z__[1] = 0.f;
                    z__[2] = d__[is + is * d_dim1];
                    z__[3] = 0.f;
                    z__[8] = 0.f;
                    z__[9] = a[is + is * a_dim1];
                    z__[10] = 0.f;
                    z__[11] = d__[is + is * d_dim1];
                    z__[16] = -b[js + js * b_dim1];
                    z__[17] = -b[js + jsp1 * b_dim1];
                    z__[18] = -e[js + js * e_dim1];
                    z__[19] = -e[js + jsp1 * e_dim1];
                    z__[24] = -b[jsp1 + js * b_dim1];
                    z__[25] = -b[jsp1 + jsp1 * b_dim1];
                    z__[26] = 0.f;
                    z__[27] = -e[jsp1 + jsp1 * e_dim1];
                    /* Set up right hand side(s) */
                    rhs[0] = c__[is + js * c_dim1];
                    rhs[1] = c__[is + jsp1 * c_dim1];
                    rhs[2] = f[is + js * f_dim1];
                    rhs[3] = f[is + jsp1 * f_dim1];
                    /* Solve Z * x = RHS */
                    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
                    if (ierr > 0)
                    {
                        *info = ierr;
                    }
                    if (*ijob == 0)
                    {
                        sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
                        if (scaloc != 1.f)
                        {
                            i__2 = *n;
                            for (k = 1;
                                    k <= i__2;
                                    ++k)
                            {
                                sscal_(m, &scaloc, &c__[k * c_dim1 + 1], & c__1);
                                sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                                /* L60: */
                            }
                            *scale *= scaloc;
                        }
                    }
                    else
                    {
                        slatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, ipiv, jpiv);
                    }
                    /* Unpack solution vector(s) */
                    c__[is + js * c_dim1] = rhs[0];
                    c__[is + jsp1 * c_dim1] = rhs[1];
                    f[is + js * f_dim1] = rhs[2];
                    f[is + jsp1 * f_dim1] = rhs[3];
                    /* Substitute R(I, J) and L(I, J) into remaining */
                    /* equation. */
                    if (i__ > 1)
                    {
                        i__2 = is - 1;
                        sger_(&i__2, &nb, &c_b27, &a[is * a_dim1 + 1], &c__1, rhs, &c__1, &c__[js * c_dim1 + 1], ldc);
                        i__2 = is - 1;
                        sger_(&i__2, &nb, &c_b27, &d__[is * d_dim1 + 1], & c__1, rhs, &c__1, &f[js * f_dim1 + 1], ldf);
                    }
                    if (j < q)
                    {
                        i__2 = *n - je;
                        saxpy_(&i__2, &rhs[2], &b[js + (je + 1) * b_dim1], ldb, &c__[is + (je + 1) * c_dim1], ldc);
                        i__2 = *n - je;
                        saxpy_(&i__2, &rhs[2], &e[js + (je + 1) * e_dim1], lde, &f[is + (je + 1) * f_dim1], ldf);
                        i__2 = *n - je;
                        saxpy_(&i__2, &rhs[3], &b[jsp1 + (je + 1) * b_dim1], ldb, &c__[is + (je + 1) * c_dim1], ldc);
                        i__2 = *n - je;
                        saxpy_(&i__2, &rhs[3], &e[jsp1 + (je + 1) * e_dim1], lde, &f[is + (je + 1) * f_dim1], ldf);
                    }
                }
                else if (mb == 2 && nb == 1)
                {
                    /* Build a 4-by-4 system Z * x = RHS */
                    z__[0] = a[is + is * a_dim1];
                    z__[1] = a[isp1 + is * a_dim1];
                    z__[2] = d__[is + is * d_dim1];
                    z__[3] = 0.f;
                    z__[8] = a[is + isp1 * a_dim1];
                    z__[9] = a[isp1 + isp1 * a_dim1];
                    z__[10] = d__[is + isp1 * d_dim1];
                    z__[11] = d__[isp1 + isp1 * d_dim1];
                    z__[16] = -b[js + js * b_dim1];
                    z__[17] = 0.f;
                    z__[18] = -e[js + js * e_dim1];
                    z__[19] = 0.f;
                    z__[24] = 0.f;
                    z__[25] = -b[js + js * b_dim1];
                    z__[26] = 0.f;
                    z__[27] = -e[js + js * e_dim1];
                    /* Set up right hand side(s) */
                    rhs[0] = c__[is + js * c_dim1];
                    rhs[1] = c__[isp1 + js * c_dim1];
                    rhs[2] = f[is + js * f_dim1];
                    rhs[3] = f[isp1 + js * f_dim1];
                    /* Solve Z * x = RHS */
                    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
                    if (ierr > 0)
                    {
                        *info = ierr;
                    }
                    if (*ijob == 0)
                    {
                        sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
                        if (scaloc != 1.f)
                        {
                            i__2 = *n;
                            for (k = 1;
                                    k <= i__2;
                                    ++k)
                            {
                                sscal_(m, &scaloc, &c__[k * c_dim1 + 1], & c__1);
                                sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                                /* L70: */
                            }
                            *scale *= scaloc;
                        }
                    }
                    else
                    {
                        slatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, ipiv, jpiv);
                    }
                    /* Unpack solution vector(s) */
                    c__[is + js * c_dim1] = rhs[0];
                    c__[isp1 + js * c_dim1] = rhs[1];
                    f[is + js * f_dim1] = rhs[2];
                    f[isp1 + js * f_dim1] = rhs[3];
                    /* Substitute R(I, J) and L(I, J) into remaining */
                    /* equation. */
                    if (i__ > 1)
                    {
                        i__2 = is - 1;
                        sgemv_("N", &i__2, &mb, &c_b27, &a[is * a_dim1 + 1], lda, rhs, &c__1, &c_b42, &c__[js * c_dim1 + 1] , &c__1);
                        i__2 = is - 1;
                        sgemv_("N", &i__2, &mb, &c_b27, &d__[is * d_dim1 + 1], ldd, rhs, &c__1, &c_b42, &f[js * f_dim1 + 1], &c__1);
                    }
                    if (j < q)
                    {
                        i__2 = *n - je;
                        sger_(&mb, &i__2, &c_b42, &rhs[2], &c__1, &b[js + (je + 1) * b_dim1], ldb, &c__[is + (je + 1) * c_dim1], ldc);
                        i__2 = *n - je;
                        sger_(&mb, &i__2, &c_b42, &rhs[2], &c__1, &e[js + (je + 1) * e_dim1], lde, &f[is + (je + 1) * f_dim1], ldf);
                    }
                }
                else if (mb == 2 && nb == 2)
                {
                    /* Build an 8-by-8 system Z * x = RHS */
                    slaset_("F", &c__8, &c__8, &c_b56, &c_b56, z__, &c__8);
                    z__[0] = a[is + is * a_dim1];
                    z__[1] = a[isp1 + is * a_dim1];
                    z__[4] = d__[is + is * d_dim1];
                    z__[8] = a[is + isp1 * a_dim1];
                    z__[9] = a[isp1 + isp1 * a_dim1];
                    z__[12] = d__[is + isp1 * d_dim1];
                    z__[13] = d__[isp1 + isp1 * d_dim1];
                    z__[18] = a[is + is * a_dim1];
                    z__[19] = a[isp1 + is * a_dim1];
                    z__[22] = d__[is + is * d_dim1];
                    z__[26] = a[is + isp1 * a_dim1];
                    z__[27] = a[isp1 + isp1 * a_dim1];
                    z__[30] = d__[is + isp1 * d_dim1];
                    z__[31] = d__[isp1 + isp1 * d_dim1];
                    z__[32] = -b[js + js * b_dim1];
                    z__[34] = -b[js + jsp1 * b_dim1];
                    z__[36] = -e[js + js * e_dim1];
                    z__[38] = -e[js + jsp1 * e_dim1];
                    z__[41] = -b[js + js * b_dim1];
                    z__[43] = -b[js + jsp1 * b_dim1];
                    z__[45] = -e[js + js * e_dim1];
                    z__[47] = -e[js + jsp1 * e_dim1];
                    z__[48] = -b[jsp1 + js * b_dim1];
                    z__[50] = -b[jsp1 + jsp1 * b_dim1];
                    z__[54] = -e[jsp1 + jsp1 * e_dim1];
                    z__[57] = -b[jsp1 + js * b_dim1];
                    z__[59] = -b[jsp1 + jsp1 * b_dim1];
                    z__[63] = -e[jsp1 + jsp1 * e_dim1];
                    /* Set up right hand side(s) */
                    k = 1;
                    ii = mb * nb + 1;
                    i__2 = nb - 1;
                    for (jj = 0;
                            jj <= i__2;
                            ++jj)
                    {
                        scopy_(&mb, &c__[is + (js + jj) * c_dim1], &c__1, & rhs[k - 1], &c__1);
                        scopy_(&mb, &f[is + (js + jj) * f_dim1], &c__1, &rhs[ ii - 1], &c__1);
                        k += mb;
                        ii += mb;
                        /* L80: */
                    }
                    /* Solve Z * x = RHS */
                    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
                    if (ierr > 0)
                    {
                        *info = ierr;
                    }
                    if (*ijob == 0)
                    {
                        sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
                        if (scaloc != 1.f)
                        {
                            i__2 = *n;
                            for (k = 1;
                                    k <= i__2;
                                    ++k)
                            {
                                sscal_(m, &scaloc, &c__[k * c_dim1 + 1], & c__1);
                                sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                                /* L90: */
                            }
                            *scale *= scaloc;
                        }
                    }
                    else
                    {
                        slatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, ipiv, jpiv);
                    }
                    /* Unpack solution vector(s) */
                    k = 1;
                    ii = mb * nb + 1;
                    i__2 = nb - 1;
                    for (jj = 0;
                            jj <= i__2;
                            ++jj)
                    {
                        scopy_(&mb, &rhs[k - 1], &c__1, &c__[is + (js + jj) * c_dim1], &c__1);
                        scopy_(&mb, &rhs[ii - 1], &c__1, &f[is + (js + jj) * f_dim1], &c__1);
                        k += mb;
                        ii += mb;
                        /* L100: */
                    }
                    /* Substitute R(I, J) and L(I, J) into remaining */
                    /* equation. */
                    if (i__ > 1)
                    {
                        i__2 = is - 1;
                        sgemm_("N", "N", &i__2, &nb, &mb, &c_b27, &a[is * a_dim1 + 1], lda, rhs, &mb, &c_b42, &c__[js * c_dim1 + 1], ldc);
                        i__2 = is - 1;
                        sgemm_("N", "N", &i__2, &nb, &mb, &c_b27, &d__[is * d_dim1 + 1], ldd, rhs, &mb, &c_b42, &f[js * f_dim1 + 1], ldf);
                    }
                    if (j < q)
                    {
                        k = mb * nb + 1;
                        i__2 = *n - je;
                        sgemm_("N", "N", &mb, &i__2, &nb, &c_b42, &rhs[k - 1], &mb, &b[js + (je + 1) * b_dim1], ldb, &c_b42, &c__[is + (je + 1) * c_dim1], ldc);
                        i__2 = *n - je;
                        sgemm_("N", "N", &mb, &i__2, &nb, &c_b42, &rhs[k - 1], &mb, &e[js + (je + 1) * e_dim1], lde, &c_b42, &f[is + (je + 1) * f_dim1], ldf);
                    }
                }
                /* L110: */
            }
            /* L120: */
        }
    }
    else
    {
        /* Solve (I, J) - subsystem */
        /* A(I, I)**T * R(I, J) + D(I, I)**T * L(J, J) = C(I, J) */
        /* R(I, I) * B(J, J) + L(I, J) * E(J, J) = -F(I, J) */
        /* for I = 1, 2, ..., P, J = Q, Q - 1, ..., 1 */
        *scale = 1.f;
        scaloc = 1.f;
        i__1 = p;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            is = iwork[i__];
            isp1 = is + 1;
            ie = iwork[i__ + 1] - 1;
            mb = ie - is + 1;
            i__2 = p + 2;
            for (j = q;
                    j >= i__2;
                    --j)
            {
                js = iwork[j];
                jsp1 = js + 1;
                je = iwork[j + 1] - 1;
                nb = je - js + 1;
                zdim = mb * nb << 1;
                if (mb == 1 && nb == 1)
                {
                    /* Build a 2-by-2 system Z**T * x = RHS */
                    z__[0] = a[is + is * a_dim1];
                    z__[1] = -b[js + js * b_dim1];
                    z__[8] = d__[is + is * d_dim1];
                    z__[9] = -e[js + js * e_dim1];
                    /* Set up right hand side(s) */
                    rhs[0] = c__[is + js * c_dim1];
                    rhs[1] = f[is + js * f_dim1];
                    /* Solve Z**T * x = RHS */
                    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
                    if (ierr > 0)
                    {
                        *info = ierr;
                    }
                    sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
                    if (scaloc != 1.f)
                    {
                        i__3 = *n;
                        for (k = 1;
                                k <= i__3;
                                ++k)
                        {
                            sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
                            sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                            /* L130: */
                        }
                        *scale *= scaloc;
                    }
                    /* Unpack solution vector(s) */
                    c__[is + js * c_dim1] = rhs[0];
                    f[is + js * f_dim1] = rhs[1];
                    /* Substitute R(I, J) and L(I, J) into remaining */
                    /* equation. */
                    if (j > p + 2)
                    {
                        alpha = rhs[0];
                        i__3 = js - 1;
                        saxpy_(&i__3, &alpha, &b[js * b_dim1 + 1], &c__1, &f[ is + f_dim1], ldf);
                        alpha = rhs[1];
                        i__3 = js - 1;
                        saxpy_(&i__3, &alpha, &e[js * e_dim1 + 1], &c__1, &f[ is + f_dim1], ldf);
                    }
                    if (i__ < p)
                    {
                        alpha = -rhs[0];
                        i__3 = *m - ie;
                        saxpy_(&i__3, &alpha, &a[is + (ie + 1) * a_dim1], lda, &c__[ie + 1 + js * c_dim1], &c__1);
                        alpha = -rhs[1];
                        i__3 = *m - ie;
                        saxpy_(&i__3, &alpha, &d__[is + (ie + 1) * d_dim1], ldd, &c__[ie + 1 + js * c_dim1], &c__1);
                    }
                }
                else if (mb == 1 && nb == 2)
                {
                    /* Build a 4-by-4 system Z**T * x = RHS */
                    z__[0] = a[is + is * a_dim1];
                    z__[1] = 0.f;
                    z__[2] = -b[js + js * b_dim1];
                    z__[3] = -b[jsp1 + js * b_dim1];
                    z__[8] = 0.f;
                    z__[9] = a[is + is * a_dim1];
                    z__[10] = -b[js + jsp1 * b_dim1];
                    z__[11] = -b[jsp1 + jsp1 * b_dim1];
                    z__[16] = d__[is + is * d_dim1];
                    z__[17] = 0.f;
                    z__[18] = -e[js + js * e_dim1];
                    z__[19] = 0.f;
                    z__[24] = 0.f;
                    z__[25] = d__[is + is * d_dim1];
                    z__[26] = -e[js + jsp1 * e_dim1];
                    z__[27] = -e[jsp1 + jsp1 * e_dim1];
                    /* Set up right hand side(s) */
                    rhs[0] = c__[is + js * c_dim1];
                    rhs[1] = c__[is + jsp1 * c_dim1];
                    rhs[2] = f[is + js * f_dim1];
                    rhs[3] = f[is + jsp1 * f_dim1];
                    /* Solve Z**T * x = RHS */
                    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
                    if (ierr > 0)
                    {
                        *info = ierr;
                    }
                    sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
                    if (scaloc != 1.f)
                    {
                        i__3 = *n;
                        for (k = 1;
                                k <= i__3;
                                ++k)
                        {
                            sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
                            sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                            /* L140: */
                        }
                        *scale *= scaloc;
                    }
                    /* Unpack solution vector(s) */
                    c__[is + js * c_dim1] = rhs[0];
                    c__[is + jsp1 * c_dim1] = rhs[1];
                    f[is + js * f_dim1] = rhs[2];
                    f[is + jsp1 * f_dim1] = rhs[3];
                    /* Substitute R(I, J) and L(I, J) into remaining */
                    /* equation. */
                    if (j > p + 2)
                    {
                        i__3 = js - 1;
                        saxpy_(&i__3, rhs, &b[js * b_dim1 + 1], &c__1, &f[is + f_dim1], ldf);
                        i__3 = js - 1;
                        saxpy_(&i__3, &rhs[1], &b[jsp1 * b_dim1 + 1], &c__1, & f[is + f_dim1], ldf);
                        i__3 = js - 1;
                        saxpy_(&i__3, &rhs[2], &e[js * e_dim1 + 1], &c__1, &f[ is + f_dim1], ldf);
                        i__3 = js - 1;
                        saxpy_(&i__3, &rhs[3], &e[jsp1 * e_dim1 + 1], &c__1, & f[is + f_dim1], ldf);
                    }
                    if (i__ < p)
                    {
                        i__3 = *m - ie;
                        sger_(&i__3, &nb, &c_b27, &a[is + (ie + 1) * a_dim1], lda, rhs, &c__1, &c__[ie + 1 + js * c_dim1], ldc);
                        i__3 = *m - ie;
                        sger_(&i__3, &nb, &c_b27, &d__[is + (ie + 1) * d_dim1] , ldd, &rhs[2], &c__1, &c__[ie + 1 + js * c_dim1], ldc);
                    }
                }
                else if (mb == 2 && nb == 1)
                {
                    /* Build a 4-by-4 system Z**T * x = RHS */
                    z__[0] = a[is + is * a_dim1];
                    z__[1] = a[is + isp1 * a_dim1];
                    z__[2] = -b[js + js * b_dim1];
                    z__[3] = 0.f;
                    z__[8] = a[isp1 + is * a_dim1];
                    z__[9] = a[isp1 + isp1 * a_dim1];
                    z__[10] = 0.f;
                    z__[11] = -b[js + js * b_dim1];
                    z__[16] = d__[is + is * d_dim1];
                    z__[17] = d__[is + isp1 * d_dim1];
                    z__[18] = -e[js + js * e_dim1];
                    z__[19] = 0.f;
                    z__[24] = 0.f;
                    z__[25] = d__[isp1 + isp1 * d_dim1];
                    z__[26] = 0.f;
                    z__[27] = -e[js + js * e_dim1];
                    /* Set up right hand side(s) */
                    rhs[0] = c__[is + js * c_dim1];
                    rhs[1] = c__[isp1 + js * c_dim1];
                    rhs[2] = f[is + js * f_dim1];
                    rhs[3] = f[isp1 + js * f_dim1];
                    /* Solve Z**T * x = RHS */
                    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
                    if (ierr > 0)
                    {
                        *info = ierr;
                    }
                    sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
                    if (scaloc != 1.f)
                    {
                        i__3 = *n;
                        for (k = 1;
                                k <= i__3;
                                ++k)
                        {
                            sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
                            sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                            /* L150: */
                        }
                        *scale *= scaloc;
                    }
                    /* Unpack solution vector(s) */
                    c__[is + js * c_dim1] = rhs[0];
                    c__[isp1 + js * c_dim1] = rhs[1];
                    f[is + js * f_dim1] = rhs[2];
                    f[isp1 + js * f_dim1] = rhs[3];
                    /* Substitute R(I, J) and L(I, J) into remaining */
                    /* equation. */
                    if (j > p + 2)
                    {
                        i__3 = js - 1;
                        sger_(&mb, &i__3, &c_b42, rhs, &c__1, &b[js * b_dim1 + 1], &c__1, &f[is + f_dim1], ldf);
                        i__3 = js - 1;
                        sger_(&mb, &i__3, &c_b42, &rhs[2], &c__1, &e[js * e_dim1 + 1], &c__1, &f[is + f_dim1], ldf);
                    }
                    if (i__ < p)
                    {
                        i__3 = *m - ie;
                        sgemv_("T", &mb, &i__3, &c_b27, &a[is + (ie + 1) * a_dim1], lda, rhs, &c__1, &c_b42, &c__[ie + 1 + js * c_dim1], &c__1);
                        i__3 = *m - ie;
                        sgemv_("T", &mb, &i__3, &c_b27, &d__[is + (ie + 1) * d_dim1], ldd, &rhs[2], &c__1, &c_b42, &c__[ie + 1 + js * c_dim1], &c__1);
                    }
                }
                else if (mb == 2 && nb == 2)
                {
                    /* Build an 8-by-8 system Z**T * x = RHS */
                    slaset_("F", &c__8, &c__8, &c_b56, &c_b56, z__, &c__8);
                    z__[0] = a[is + is * a_dim1];
                    z__[1] = a[is + isp1 * a_dim1];
                    z__[4] = -b[js + js * b_dim1];
                    z__[6] = -b[jsp1 + js * b_dim1];
                    z__[8] = a[isp1 + is * a_dim1];
                    z__[9] = a[isp1 + isp1 * a_dim1];
                    z__[13] = -b[js + js * b_dim1];
                    z__[15] = -b[jsp1 + js * b_dim1];
                    z__[18] = a[is + is * a_dim1];
                    z__[19] = a[is + isp1 * a_dim1];
                    z__[20] = -b[js + jsp1 * b_dim1];
                    z__[22] = -b[jsp1 + jsp1 * b_dim1];
                    z__[26] = a[isp1 + is * a_dim1];
                    z__[27] = a[isp1 + isp1 * a_dim1];
                    z__[29] = -b[js + jsp1 * b_dim1];
                    z__[31] = -b[jsp1 + jsp1 * b_dim1];
                    z__[32] = d__[is + is * d_dim1];
                    z__[33] = d__[is + isp1 * d_dim1];
                    z__[36] = -e[js + js * e_dim1];
                    z__[41] = d__[isp1 + isp1 * d_dim1];
                    z__[45] = -e[js + js * e_dim1];
                    z__[50] = d__[is + is * d_dim1];
                    z__[51] = d__[is + isp1 * d_dim1];
                    z__[52] = -e[js + jsp1 * e_dim1];
                    z__[54] = -e[jsp1 + jsp1 * e_dim1];
                    z__[59] = d__[isp1 + isp1 * d_dim1];
                    z__[61] = -e[js + jsp1 * e_dim1];
                    z__[63] = -e[jsp1 + jsp1 * e_dim1];
                    /* Set up right hand side(s) */
                    k = 1;
                    ii = mb * nb + 1;
                    i__3 = nb - 1;
                    for (jj = 0;
                            jj <= i__3;
                            ++jj)
                    {
                        scopy_(&mb, &c__[is + (js + jj) * c_dim1], &c__1, & rhs[k - 1], &c__1);
                        scopy_(&mb, &f[is + (js + jj) * f_dim1], &c__1, &rhs[ ii - 1], &c__1);
                        k += mb;
                        ii += mb;
                        /* L160: */
                    }
                    /* Solve Z**T * x = RHS */
                    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
                    if (ierr > 0)
                    {
                        *info = ierr;
                    }
                    sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
                    if (scaloc != 1.f)
                    {
                        i__3 = *n;
                        for (k = 1;
                                k <= i__3;
                                ++k)
                        {
                            sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
                            sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
                            /* L170: */
                        }
                        *scale *= scaloc;
                    }
                    /* Unpack solution vector(s) */
                    k = 1;
                    ii = mb * nb + 1;
                    i__3 = nb - 1;
                    for (jj = 0;
                            jj <= i__3;
                            ++jj)
                    {
                        scopy_(&mb, &rhs[k - 1], &c__1, &c__[is + (js + jj) * c_dim1], &c__1);
                        scopy_(&mb, &rhs[ii - 1], &c__1, &f[is + (js + jj) * f_dim1], &c__1);
                        k += mb;
                        ii += mb;
                        /* L180: */
                    }
                    /* Substitute R(I, J) and L(I, J) into remaining */
                    /* equation. */
                    if (j > p + 2)
                    {
                        i__3 = js - 1;
                        sgemm_("N", "T", &mb, &i__3, &nb, &c_b42, &c__[is + js * c_dim1], ldc, &b[js * b_dim1 + 1], ldb, & c_b42, &f[is + f_dim1], ldf);
                        i__3 = js - 1;
                        sgemm_("N", "T", &mb, &i__3, &nb, &c_b42, &f[is + js * f_dim1], ldf, &e[js * e_dim1 + 1], lde, & c_b42, &f[is + f_dim1], ldf);
                    }
                    if (i__ < p)
                    {
                        i__3 = *m - ie;
                        sgemm_("T", "N", &i__3, &nb, &mb, &c_b27, &a[is + (ie + 1) * a_dim1], lda, &c__[is + js * c_dim1], ldc, &c_b42, &c__[ie + 1 + js * c_dim1], ldc);
                        i__3 = *m - ie;
                        sgemm_("T", "N", &i__3, &nb, &mb, &c_b27, &d__[is + ( ie + 1) * d_dim1], ldd, &f[is + js * f_dim1], ldf, &c_b42, &c__[ie + 1 + js * c_dim1], ldc);
                    }
                }
                /* L190: */
            }
            /* L200: */
        }
    }
    return 0;
    /* End of STGSY2 */
}
/* stgsy2_ */
