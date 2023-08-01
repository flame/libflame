/* dtrsyl3.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b19 = 2.;
static doublereal c_b31 = -1.;
static doublereal c_b32 = 1.;
/* > \brief \b DTRSYL3 */
/* Definition: */
/* =========== */
/* > \par Purpose */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTRSYL3 solves the real Sylvester matrix equation: */
/* > */
/* > op(A)*X + X*op(B) = scale*C or */
/* > op(A)*X - X*op(B) = scale*C, */
/* > */
/* > where op(A) = A or A**T, and A and B are both upper quasi- */
/* > triangular. A is M-by-M and B is N-by-N;
the right hand side C and */
/* > the solution X are M-by-N;
and scale is an output scale factor, set */
/* > <= 1 to avoid overflow in X. */
/* > */
/* > A and B must be in Schur canonical form (as returned by DHSEQR), that */
/* > is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
*/
/* > each 2-by-2 diagonal block has its diagonal elements equal and its */
/* > off-diagonal elements of opposite sign. */
/* > */
/* > This is the block version of the algorithm. */
/* > \endverbatim */
/* Arguments */
/* ========= */
/* > \param[in] TRANA */
/* > \verbatim */
/* > TRANA is CHARACTER*1 */
/* > Specifies the option op(A): */
/* > = 'N': op(A) = A (No transpose) */
/* > = 'T': op(A) = A**T (Transpose) */
/* > = 'C': op(A) = A**H (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] TRANB */
/* > \verbatim */
/* > TRANB is CHARACTER*1 */
/* > Specifies the option op(B): */
/* > = 'N': op(B) = B (No transpose) */
/* > = 'T': op(B) = B**T (Transpose) */
/* > = 'C': op(B) = B**H (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] ISGN */
/* > \verbatim */
/* > ISGN is INTEGER */
/* > Specifies the sign in the equation: */
/* > = +1: solve op(A)*X + X*op(B) = scale*C */
/* > = -1: solve op(A)*X - X*op(B) = scale*C */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The order of the matrix A, and the number of rows in the */
/* > matrices X and C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix B, and the number of columns in the */
/* > matrices X and C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,M) */
/* > The upper quasi-triangular matrix A, in Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,N) */
/* > The upper quasi-triangular matrix B, in Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (LDC,N) */
/* > On entry, the M-by-N right hand side matrix C. */
/* > On exit, C is overwritten by the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M) */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is DOUBLE PRECISION */
/* > The scale factor, scale, set <= 1 to avoid overflow in X. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (fla_max(1,LIWORK)) */
/* > On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > IWORK is INTEGER */
/* > The dimension of the array IWORK. LIWORK >= ((M + NB - 1) / NB + 1) */
/* > + ((N + NB - 1) / NB + 1), where NB is the optimal block size. */
/* > */
/* > If LIWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal dimension of the IWORK array, */
/* > returns this value as the first entry of the IWORK array, and */
/* > no error message related to LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] SWORK */
/* > \verbatim */
/* > SWORK is DOUBLE PRECISION array, dimension (fla_max(2, ROWS), */
/* > fla_max(1,COLS)). */
/* > On exit, if INFO = 0, SWORK(1) returns the optimal value ROWS */
/* > and SWORK(2) returns the optimal COLS. */
/* > \endverbatim */
/* > */
/* > \param[in] LDSWORK */
/* > \verbatim */
/* > LDSWORK is INTEGER */
/* > LDSWORK >= fla_max(2,ROWS), where ROWS = ((M + NB - 1) / NB + 1) */
/* > and NB is the optimal block size. */
/* > */
/* > If LDSWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal dimensions of the SWORK matrix, */
/* > returns these values as the first and second entry of the SWORK */
/* > matrix, and no error message related LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > = 1: A and B have common or very close eigenvalues;
perturbed */
/* > values were used to solve the equation (but the matrices */
/* > A and B are unchanged). */
/* > \endverbatim */
/* ===================================================================== */
/* References: */
/* E. S. Quintana-Orti and R. A. Van De Geijn (2003). Formal derivation of */
/* algorithms: The triangular Sylvester equation, ACM Transactions */
/* on Mathematical Software (TOMS), volume 29, pages 218--243. */
/* A. Schwarz and C. C. Kjelgaard Mikkelsen (2020). Robust Task-Parallel */
/* Solution of the Triangular Sylvester Equation. Lecture Notes in */
/* Computer Science, vol 12043, pages 82--92, Springer. */
/* Contributor: */
/* Angelika Schwarz, Umea University, Sweden. */
/* ===================================================================== */
/* Subroutine */
int dtrsyl3_(char *trana, char *tranb, integer *isgn, integer *m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, doublereal *scale, integer *iwork, integer *liwork, doublereal *swork, integer *ldswork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dtrsyl3 inputs: trana %c, tranb %c, isgn %" FLA_IS ", m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS ", ldc %" FLA_IS "",*trana, *tranb, *isgn, *m, *n, *lda, *ldb, *ldc);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, swork_dim1, swork_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3;
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    /* Local variables */
    integer i__, j, k, l, i1, i2, j1, j2, k1, k2, l1, l2, nb, pc, jj, ll, nba, nbb;
    doublereal buf, sgn, scal, anrm, bnrm, cnrm;
    integer awrk, bwrk;
    int temp;
    logical skip;
    doublereal *wnrm, xnrm;
    extern /* Subroutine */
    int dscal_(integer *, doublereal *, doublereal *, integer *), dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *);
    integer iinfo;
    extern doublereal dlamch_(char *), dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *);
    extern /* Subroutine */
    int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *);
    doublereal scaloc, scamin;
    extern doublereal dlarmm_(doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    doublereal bignum;
    logical notrna, notrnb;
    doublereal smlnum;
    logical lquery;
    extern /* Subroutine */
    int dtrsyl_(char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
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
    /* Decode and Test input parameters */
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
    --iwork;
    swork_dim1 = *ldswork;
    swork_offset = 1 + swork_dim1;
    swork -= swork_offset;
    wnrm = (doublereal *)malloc(fla_max(*m,*n) * sizeof(doublereal));
    /* Function Body */
    notrna = lsame_(trana, "N");
    notrnb = lsame_(tranb, "N");
    /* Use the same block size for all matrices. */
    /* Computing fla_max */
    i__1 = 8;
    i__2 = ilaenv_(&c__1, "DTRSYL", "", m, n, &c_n1, &c_n1); // , expr subst
    nb = fla_max(i__1,i__2);
    /* Compute number of blocks in A and B */
    /* Computing fla_max */
    i__1 = 1;
    i__2 = (*m + nb - 1) / nb; // , expr subst
    nba = fla_max(i__1,i__2);
    /* Computing fla_max */
    i__1 = 1;
    i__2 = (*n + nb - 1) / nb; // , expr subst
    nbb = fla_max(i__1,i__2);
    /* Compute workspace */
    *info = 0;
    lquery = *liwork == -1 || *ldswork == -1;
    iwork[1] = nba + nbb + 2;
    if (lquery)
    {
        *ldswork = 2;
        swork[swork_dim1 + 1] = (doublereal) fla_max(nba,nbb);
        swork[swork_dim1 + 2] = (doublereal) ((nbb << 1) + nba);
    }
    /* Test the input arguments */
    if (! notrna && ! lsame_(trana, "T") && ! lsame_( trana, "C"))
    {
        *info = -1;
    }
    else if (! notrnb && ! lsame_(tranb, "T") && ! lsame_(tranb, "C"))
    {
        *info = -2;
    }
    else if (*isgn != 1 && *isgn != -1)
    {
        *info = -3;
    }
    else if (*m < 0)
    {
        *info = -4;
    }
    else if (*n < 0)
    {
        *info = -5;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -7;
    }
    else if (*ldb < fla_max(1,*n))
    {
        *info = -9;
    }
    else if (*ldc < fla_max(1,*m))
    {
        *info = -11;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DTRSYL3", &i__1, (ftnlen)7);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (lquery)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    *scale = 1.;
    if (*m == 0 || *n == 0)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Use unblocked code for small problems or if insufficient */
    /* workspaces are provided */
    if (fla_min(nba,nbb) == 1 || *ldswork < fla_max(nba,nbb) || *liwork < iwork[1])
    {
        dtrsyl_(trana, tranb, isgn, m, n, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset], ldc, scale, info);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Set constants to control overflow */
    smlnum = dlamch_("S");
    bignum = 1. / smlnum;
    /* Partition A such that 2-by-2 blocks on the diagonal are not split */
    skip = FALSE_;
    i__1 = nba;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        iwork[i__] = (i__ - 1) * nb + 1;
    }
    iwork[nba + 1] = *m + 1;
    i__1 = nba;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        l1 = iwork[k];
        l2 = iwork[k + 1] - 1;
        i__2 = l2;
        for (l = l1;
                l <= i__2;
                ++l)
        {
            if (skip)
            {
                skip = FALSE_;
                continue;
            }
            if (l >= *m)
            {
                /* A( M, M ) is a 1-by-1 block */
                continue;
            }
            if (a[l + (l + 1) * a_dim1] != 0. && a[l + 1 + l * a_dim1] != 0.)
            {
                /* Check if 2-by-2 block is split */
                if (l + 1 == iwork[k + 1])
                {
                    ++iwork[k + 1];
                    continue;
                }
                skip = TRUE_;
            }
        }
    }
    iwork[nba + 1] = *m + 1;
    if (iwork[nba] >= iwork[nba + 1])
    {
        iwork[nba] = iwork[nba + 1];
        --nba;
    }
    /* Partition B such that 2-by-2 blocks on the diagonal are not split */
    pc = nba + 1;
    skip = FALSE_;
    i__1 = nbb;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        iwork[pc + i__] = (i__ - 1) * nb + 1;
    }
    iwork[pc + nbb + 1] = *n + 1;
    i__1 = nbb;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        l1 = iwork[pc + k];
        l2 = iwork[pc + k + 1] - 1;
        i__2 = l2;
        for (l = l1;
                l <= i__2;
                ++l)
        {
            if (skip)
            {
                skip = FALSE_;
                continue;
            }
            if (l >= *n)
            {
                /* B( N, N ) is a 1-by-1 block */
                continue;
            }
            if (b[l + (l + 1) * b_dim1] != 0. && b[l + 1 + l * b_dim1] != 0.)
            {
                /* Check if 2-by-2 block is split */
                if (l + 1 == iwork[pc + k + 1])
                {
                    ++iwork[pc + k + 1];
                    continue;
                }
                skip = TRUE_;
            }
        }
    }
    iwork[pc + nbb + 1] = *n + 1;
    if (iwork[pc + nbb] >= iwork[pc + nbb + 1])
    {
        iwork[pc + nbb] = iwork[pc + nbb + 1];
        --nbb;
    }
    /* Set local scaling factors - must never attain zero. */
    i__1 = nbb;
    for (l = 1;
            l <= i__1;
            ++l)
    {
        i__2 = nba;
        for (k = 1;
                k <= i__2;
                ++k)
        {
            swork[k + l * swork_dim1] = 1.;
        }
    }
    /* Fallback scaling factor to prevent flushing of SWORK( K, L ) to zero. */
    /* This scaling is to ensure compatibility with TRSYL and may get flushed. */
    buf = 1.;
    /* Compute upper bounds of blocks of A and B */
    awrk = nbb;
    i__1 = nba;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        k1 = iwork[k];
        k2 = iwork[k + 1];
        i__2 = nba;
        for (l = k;
                l <= i__2;
                ++l)
        {
            l1 = iwork[l];
            l2 = iwork[l + 1];
            if (notrna)
            {
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                swork[k + (awrk + l) * swork_dim1] = dlange_("I", &i__3, & i__4, &a[k1 + l1 * a_dim1], lda, wnrm);
            }
            else
            {
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                swork[l + (awrk + k) * swork_dim1] = dlange_("1", &i__3, & i__4, &a[k1 + l1 * a_dim1], lda, wnrm);
            }
        }
    }
    bwrk = nbb + nba;
    i__1 = nbb;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        k1 = iwork[pc + k];
        k2 = iwork[pc + k + 1];
        i__2 = nbb;
        for (l = k;
                l <= i__2;
                ++l)
        {
            l1 = iwork[pc + l];
            l2 = iwork[pc + l + 1];
            if (notrnb)
            {
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                swork[k + (bwrk + l) * swork_dim1] = dlange_("I", &i__3, & i__4, &b[k1 + l1 * b_dim1], ldb, wnrm);
            }
            else
            {
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                swork[l + (bwrk + k) * swork_dim1] = dlange_("1", &i__3, & i__4, &b[k1 + l1 * b_dim1], ldb, wnrm);
            }
        }
    }
    sgn = (doublereal) (*isgn);
    if (notrna && notrnb)
    {
        /* Solve A*X + ISGN*X*B = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* bottom-left corner column by column by */
        /* A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */
        /* Where */
        /* M L-1 */
        /* R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]. */
        /* I=K+1 J=1 */
        /* Start loop over block rows (index = K) and block columns (index = L) */
        for (k = nba;
                k >= 1;
                --k)
        {
            /* K1: row index of the first row in X( K, L ) */
            /* K2: row index of the first row in X( K+1, L ) */
            /* so the K2 - K1 is the column count of the block X( K, L ) */
            k1 = iwork[k];
            k2 = iwork[k + 1];
            i__1 = nbb;
            for (l = 1;
                    l <= i__1;
                    ++l)
            {
                /* L1: column index of the first column in X( K, L ) */
                /* L2: column index of the first column in X( K, L + 1) */
                /* so that L2 - L1 is the row count of the block X( K, L ) */
                l1 = iwork[pc + l];
                l2 = iwork[pc + l + 1];
                i__2 = k2 - k1;
                i__3 = l2 - l1;
                dtrsyl_(trana, tranb, isgn, &i__2, &i__3, &a[k1 + k1 * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, &c__[k1 + l1 * c_dim1], ldc, &scaloc, &iinfo);
                *info = fla_max(*info,iinfo);
                if (scaloc * swork[k + l * swork_dim1] == 0.)
                {
                    if (scaloc == 0.)
                    {
                        /* The magnitude of the largest entry of X(K1:K2-1, L1:L2-1) */
                        /* is larger than the product of BIGNUM**2 and cannot be */
                        /* represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1). */
                        /* Mark the computation as pointless. */
                        buf = 0.;
                    }
                    else
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__3 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                    }
                    i__2 = nbb;
                    for (jj = 1;
                            jj <= i__2;
                            ++jj)
                    {
                        i__3 = nba;
                        for (ll = 1;
                                ll <= i__3;
                                ++ll)
                        {
                            /* Bound by BIGNUM to not introduce Inf. The value */
                            /* is irrelevant;
                            corresponding entries of the */
                            /* solution will be flushed in consistency scaling. */
                            /* Computing fla_min */
                            frexp(scaloc, (int *) &temp); d__3 = temp;
                            d__1 = bignum;
                            d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                            swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                        }
                    }
                }
                swork[k + l * swork_dim1] = scaloc * swork[k + l * swork_dim1] ;
                i__2 = k2 - k1;
                i__3 = l2 - l1;
                xnrm = dlange_("I", &i__2, &i__3, &c__[k1 + l1 * c_dim1], ldc, wnrm);
                for (i__ = k - 1;
                        i__ >= 1;
                        --i__)
                {
                    /* C( I, L ) := C( I, L ) - A( I, K ) * C( K, L ) */
                    i1 = iwork[i__];
                    i2 = iwork[i__ + 1];
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__2 = i2 - i1;
                    i__3 = l2 - l1;
                    cnrm = dlange_("I", &i__2, &i__3, &c__[i1 + l1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[i__ + l * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1,d__2);
                    cnrm *= scamin / swork[i__ + l * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    anrm = swork[i__ + (awrk + k) * swork_dim1];
                    scaloc = dlarmm_(&anrm, &xnrm, &cnrm);
                    if (scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                        i__2 = nbb;
                        for (jj = 1;
                                jj <= i__2;
                                ++jj)
                        {
                            i__3 = nba;
                            for (ll = 1;
                                    ll <= i__3;
                                    ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *) &temp); d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                            }
                        }
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        
                        scamin /= pow_dd(&c_b19, &d__1);
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        
                        scaloc /= pow_dd(&c_b19, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( I, L ) and C( K, L ). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__2 = l2 - 1;
                        for (jj = l1;
                                jj <= i__2;
                                ++jj)
                        {
                            i__3 = k2 - k1;
                            dscal_(&i__3, &scal, &c__[k1 + jj * c_dim1], & c__1);
                        }
                    }
                    scal = scamin / swork[i__ + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__2 = l2 - 1;
                        for (ll = l1;
                                ll <= i__2;
                                ++ll)
                        {
                            i__3 = i2 - i1;
                            dscal_(&i__3, &scal, &c__[i1 + ll * c_dim1], & c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[i__ + l * swork_dim1] = scamin * scaloc;
                    i__2 = i2 - i1;
                    i__3 = l2 - l1;
                    i__4 = k2 - k1;
                    dgemm_("N", "N", &i__2, &i__3, &i__4, &c_b31, &a[i1 + k1 * a_dim1], lda, &c__[k1 + l1 * c_dim1], ldc, & c_b32, &c__[i1 + l1 * c_dim1], ldc);
                }
                i__2 = nbb;
                for (j = l + 1;
                        j <= i__2;
                        ++j)
                {
                    /* C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J ) */
                    j1 = iwork[pc + j];
                    j2 = iwork[pc + j + 1];
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__3 = k2 - k1;
                    i__4 = j2 - j1;
                    cnrm = dlange_("I", &i__3, &i__4, &c__[k1 + j1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[k + j * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1,d__2);
                    cnrm *= scamin / swork[k + j * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    bnrm = swork[l + (bwrk + j) * swork_dim1];
                    scaloc = dlarmm_(&bnrm, &xnrm, &cnrm);
                    if (scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                        i__3 = nbb;
                        for (jj = 1;
                                jj <= i__3;
                                ++jj)
                        {
                            i__4 = nba;
                            for (ll = 1;
                                    ll <= i__4;
                                    ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *) &temp); d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                            }
                        }
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        scamin /= pow_dd(&c_b19, &d__1);
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        scaloc /= pow_dd(&c_b19, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( K, J ) and C( K, L). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__3 = l2 - 1;
                        for (ll = l1;
                                ll <= i__3;
                                ++ll)
                        {
                            i__4 = k2 - k1;
                            dscal_(&i__4, &scal, &c__[k1 + ll * c_dim1], & c__1);
                        }
                    }
                    scal = scamin / swork[k + j * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__3 = j2 - 1;
                        for (jj = j1;
                                jj <= i__3;
                                ++jj)
                        {
                            i__4 = k2 - k1;
                            dscal_(&i__4, &scal, &c__[k1 + jj * c_dim1], & c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[k + j * swork_dim1] = scamin * scaloc;
                    i__3 = k2 - k1;
                    i__4 = j2 - j1;
                    i__5 = l2 - l1;
                    d__1 = -sgn;
                    dgemm_("N", "N", &i__3, &i__4, &i__5, &d__1, &c__[k1 + l1 * c_dim1], ldc, &b[l1 + j1 * b_dim1], ldb, &c_b32, &c__[k1 + j1 * c_dim1], ldc);
                }
            }
        }
    }
    else if (! notrna && notrnb)
    {
        /* Solve A**T*X + ISGN*X*B = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* upper-left corner column by column by */
        /* A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */
        /* Where */
        /* K-1 L-1 */
        /* R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)] */
        /* I=1 J=1 */
        /* Start loop over block rows (index = K) and block columns (index = L) */
        i__1 = nba;
        for (k = 1;
                k <= i__1;
                ++k)
        {
            /* K1: row index of the first row in X( K, L ) */
            /* K2: row index of the first row in X( K+1, L ) */
            /* so the K2 - K1 is the column count of the block X( K, L ) */
            k1 = iwork[k];
            k2 = iwork[k + 1];
            i__2 = nbb;
            for (l = 1;
                    l <= i__2;
                    ++l)
            {
                /* L1: column index of the first column in X( K, L ) */
                /* L2: column index of the first column in X( K, L + 1) */
                /* so that L2 - L1 is the row count of the block X( K, L ) */
                l1 = iwork[pc + l];
                l2 = iwork[pc + l + 1];
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                dtrsyl_(trana, tranb, isgn, &i__3, &i__4, &a[k1 + k1 * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, &c__[k1 + l1 * c_dim1], ldc, &scaloc, &iinfo);
                *info = fla_max(*info,iinfo);
                if (scaloc * swork[k + l * swork_dim1] == 0.)
                {
                    if (scaloc == 0.)
                    {
                        /* The magnitude of the largest entry of X(K1:K2-1, L1:L2-1) */
                        /* is larger than the product of BIGNUM**2 and cannot be */
                        /* represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1). */
                        /* Mark the computation as pointless. */
                        buf = 0.;
                    }
                    else
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        
                        buf *= pow_dd(&c_b19, &d__1);
                    }
                    i__3 = nbb;
                    for (jj = 1;
                            jj <= i__3;
                            ++jj)
                    {
                        i__4 = nba;
                        for (ll = 1;
                                ll <= i__4;
                                ++ll)
                        {
                            /* Bound by BIGNUM to not introduce Inf. The value */
                            /* is irrelevant;
                            corresponding entries of the */
                            /* solution will be flushed in consistency scaling. */
                            /* Computing fla_min */
                            frexp(scaloc, (int *) &temp); d__3 = temp;
                            d__1 = bignum;
                            d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                            swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                        }
                    }
                }
                swork[k + l * swork_dim1] = scaloc * swork[k + l * swork_dim1] ;
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                xnrm = dlange_("I", &i__3, &i__4, &c__[k1 + l1 * c_dim1], ldc, wnrm);
                i__3 = nba;
                for (i__ = k + 1;
                        i__ <= i__3;
                        ++i__)
                {
                    /* C( I, L ) := C( I, L ) - A( K, I )**T * C( K, L ) */
                    i1 = iwork[i__];
                    i2 = iwork[i__ + 1];
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__4 = i2 - i1;
                    i__5 = l2 - l1;
                    cnrm = dlange_("I", &i__4, &i__5, &c__[i1 + l1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[i__ + l * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1,d__2);
                    cnrm *= scamin / swork[i__ + l * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    anrm = swork[i__ + (awrk + k) * swork_dim1];
                    scaloc = dlarmm_(&anrm, &xnrm, &cnrm);
                    if (scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                        i__4 = nbb;
                        for (jj = 1;
                                jj <= i__4;
                                ++jj)
                        {
                            i__5 = nba;
                            for (ll = 1;
                                    ll <= i__5;
                                    ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *) &temp); d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                            }
                        }
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        
                        scamin /= pow_dd(&c_b19, &d__1);
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        scaloc /= pow_dd(&c_b19, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to to C( I, L ) and C( K, L ). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__4 = l2 - 1;
                        for (ll = l1;
                                ll <= i__4;
                                ++ll)
                        {
                            i__5 = k2 - k1;
                            dscal_(&i__5, &scal, &c__[k1 + ll * c_dim1], & c__1);
                        }
                    }
                    scal = scamin / swork[i__ + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__4 = l2 - 1;
                        for (ll = l1;
                                ll <= i__4;
                                ++ll)
                        {
                            i__5 = i2 - i1;
                            dscal_(&i__5, &scal, &c__[i1 + ll * c_dim1], & c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[i__ + l * swork_dim1] = scamin * scaloc;
                    i__4 = i2 - i1;
                    i__5 = l2 - l1;
                    i__6 = k2 - k1;
                    dgemm_("T", "N", &i__4, &i__5, &i__6, &c_b31, &a[k1 + i1 * a_dim1], lda, &c__[k1 + l1 * c_dim1], ldc, & c_b32, &c__[i1 + l1 * c_dim1], ldc);
                }
                i__3 = nbb;
                for (j = l + 1;
                        j <= i__3;
                        ++j)
                {
                    /* C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J ) */
                    j1 = iwork[pc + j];
                    j2 = iwork[pc + j + 1];
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__4 = k2 - k1;
                    i__5 = j2 - j1;
                    cnrm = dlange_("I", &i__4, &i__5, &c__[k1 + j1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[k + j * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1,d__2);
                    cnrm *= scamin / swork[k + j * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    bnrm = swork[l + (bwrk + j) * swork_dim1];
                    scaloc = dlarmm_(&bnrm, &xnrm, &cnrm);
                    if (scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                        i__4 = nbb;
                        for (jj = 1;
                                jj <= i__4;
                                ++jj)
                        {
                            i__5 = nba;
                            for (ll = 1;
                                    ll <= i__5;
                                    ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *) &temp); d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                            }
                        }
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        scamin /= pow_dd(&c_b19, &d__1);
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        scaloc /= pow_dd(&c_b19, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to to C( K, J ) and C( K, L ). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__4 = l2 - 1;
                        for (ll = l1;
                                ll <= i__4;
                                ++ll)
                        {
                            i__5 = k2 - k1;
                            dscal_(&i__5, &scal, &c__[k1 + ll * c_dim1], & c__1);
                        }
                    }
                    scal = scamin / swork[k + j * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__4 = j2 - 1;
                        for (jj = j1;
                                jj <= i__4;
                                ++jj)
                        {
                            i__5 = k2 - k1;
                            dscal_(&i__5, &scal, &c__[k1 + jj * c_dim1], & c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[k + j * swork_dim1] = scamin * scaloc;
                    i__4 = k2 - k1;
                    i__5 = j2 - j1;
                    i__6 = l2 - l1;
                    d__1 = -sgn;
                    dgemm_("N", "N", &i__4, &i__5, &i__6, &d__1, &c__[k1 + l1 * c_dim1], ldc, &b[l1 + j1 * b_dim1], ldb, &c_b32, &c__[k1 + j1 * c_dim1], ldc);
                }
            }
        }
    }
    else if (! notrna && ! notrnb)
    {
        /* Solve A**T*X + ISGN*X*B**T = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* top-right corner column by column by */
        /* A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L) */
        /* Where */
        /* K-1 N */
        /* R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T]. */
        /* I=1 J=L+1 */
        /* Start loop over block rows (index = K) and block columns (index = L) */
        i__1 = nba;
        for (k = 1;
                k <= i__1;
                ++k)
        {
            /* K1: row index of the first row in X( K, L ) */
            /* K2: row index of the first row in X( K+1, L ) */
            /* so the K2 - K1 is the column count of the block X( K, L ) */
            k1 = iwork[k];
            k2 = iwork[k + 1];
            for (l = nbb;
                    l >= 1;
                    --l)
            {
                /* L1: column index of the first column in X( K, L ) */
                /* L2: column index of the first column in X( K, L + 1) */
                /* so that L2 - L1 is the row count of the block X( K, L ) */
                l1 = iwork[pc + l];
                l2 = iwork[pc + l + 1];
                i__2 = k2 - k1;
                i__3 = l2 - l1;
                dtrsyl_(trana, tranb, isgn, &i__2, &i__3, &a[k1 + k1 * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, &c__[k1 + l1 * c_dim1], ldc, &scaloc, &iinfo);
                *info = fla_max(*info,iinfo);
                swork[k + l * swork_dim1] = scaloc * swork[k + l * swork_dim1] ;
                if (scaloc * swork[k + l * swork_dim1] == 0.)
                {
                    if (scaloc == 0.)
                    {
                        /* The magnitude of the largest entry of X(K1:K2-1, L1:L2-1) */
                        /* is larger than the product of BIGNUM**2 and cannot be */
                        /* represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1). */
                        /* Mark the computation as pointless. */
                        buf = 0.;
                    }
                    else
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                    }
                    i__2 = nbb;
                    for (jj = 1;
                            jj <= i__2;
                            ++jj)
                    {
                        i__3 = nba;
                        for (ll = 1;
                                ll <= i__3;
                                ++ll)
                        {
                            /* Bound by BIGNUM to not introduce Inf. The value */
                            /* is irrelevant;
                            corresponding entries of the */
                            /* solution will be flushed in consistency scaling. */
                            /* Computing fla_min */
                            frexp(scaloc, (int *) &temp); d__3 = temp;
                            d__1 = bignum;
                            d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                            swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                        }
                    }
                }
                i__2 = k2 - k1;
                i__3 = l2 - l1;
                xnrm = dlange_("I", &i__2, &i__3, &c__[k1 + l1 * c_dim1], ldc, wnrm);
                i__2 = nba;
                for (i__ = k + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    /* C( I, L ) := C( I, L ) - A( K, I )**T * C( K, L ) */
                    i1 = iwork[i__];
                    i2 = iwork[i__ + 1];
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__3 = i2 - i1;
                    i__4 = l2 - l1;
                    cnrm = dlange_("I", &i__3, &i__4, &c__[i1 + l1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[i__ + l * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1,d__2);
                    cnrm *= scamin / swork[i__ + l * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    anrm = swork[i__ + (awrk + k) * swork_dim1];
                    scaloc = dlarmm_(&anrm, &xnrm, &cnrm);
                    if (scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                        i__3 = nbb;
                        for (jj = 1;
                                jj <= i__3;
                                ++jj)
                        {
                            i__4 = nba;
                            for (ll = 1;
                                    ll <= i__4;
                                    ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *) &temp); d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                            }
                        }
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        
                        
                        scamin /= pow_dd(&c_b19, &d__1);
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        scaloc /= pow_dd(&c_b19, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( I, L ) and C( K, L ). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__3 = l2 - 1;
                        for (ll = l1;
                                ll <= i__3;
                                ++ll)
                        {
                            i__4 = k2 - k1;
                            dscal_(&i__4, &scal, &c__[k1 + ll * c_dim1], & c__1);
                        }
                    }
                    scal = scamin / swork[i__ + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__3 = l2 - 1;
                        for (ll = l1;
                                ll <= i__3;
                                ++ll)
                        {
                            i__4 = i2 - i1;
                            dscal_(&i__4, &scal, &c__[i1 + ll * c_dim1], & c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[i__ + l * swork_dim1] = scamin * scaloc;
                    i__3 = i2 - i1;
                    i__4 = l2 - l1;
                    i__5 = k2 - k1;
                    dgemm_("T", "N", &i__3, &i__4, &i__5, &c_b31, &a[k1 + i1 * a_dim1], lda, &c__[k1 + l1 * c_dim1], ldc, & c_b32, &c__[i1 + l1 * c_dim1], ldc);
                }
                i__2 = l - 1;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    /* C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**T */
                    j1 = iwork[pc + j];
                    j2 = iwork[pc + j + 1];
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__3 = k2 - k1;
                    i__4 = j2 - j1;
                    cnrm = dlange_("I", &i__3, &i__4, &c__[k1 + j1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[k + j * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1,d__2);
                    cnrm *= scamin / swork[k + j * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    bnrm = swork[l + (bwrk + j) * swork_dim1];
                    scaloc = dlarmm_(&bnrm, &xnrm, &cnrm);
                    if (scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                        i__3 = nbb;
                        for (jj = 1;
                                jj <= i__3;
                                ++jj)
                        {
                            i__4 = nba;
                            for (ll = 1;
                                    ll <= i__4;
                                    ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *) &temp); d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                            }
                        }
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        
                        
                        scamin /= pow_dd(&c_b19, &d__1);
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        scaloc /= pow_dd(&c_b19, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( K, J ) and C( K, L ). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__3 = l2 - 1;
                        for (ll = l1;
                                ll <= i__3;
                                ++ll)
                        {
                            i__4 = k2 - k1;
                            dscal_(&i__4, &scal, &c__[k1 + ll * c_dim1], & c__1);
                        }
                    }
                    scal = scamin / swork[k + j * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__3 = j2 - 1;
                        for (jj = j1;
                                jj <= i__3;
                                ++jj)
                        {
                            i__4 = k2 - k1;
                            dscal_(&i__4, &scal, &c__[k1 + jj * c_dim1], & c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[k + j * swork_dim1] = scamin * scaloc;
                    i__3 = k2 - k1;
                    i__4 = j2 - j1;
                    i__5 = l2 - l1;
                    d__1 = -sgn;
                    dgemm_("N", "T", &i__3, &i__4, &i__5, &d__1, &c__[k1 + l1 * c_dim1], ldc, &b[j1 + l1 * b_dim1], ldb, &c_b32, &c__[k1 + j1 * c_dim1], ldc);
                }
            }
        }
    }
    else if (notrna && ! notrnb)
    {
        /* Solve A*X + ISGN*X*B**T = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* bottom-right corner column by column by */
        /* A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L) */
        /* Where */
        /* M N */
        /* R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T]. */
        /* I=K+1 J=L+1 */
        /* Start loop over block rows (index = K) and block columns (index = L) */
        for (k = nba;
                k >= 1;
                --k)
        {
            /* K1: row index of the first row in X( K, L ) */
            /* K2: row index of the first row in X( K+1, L ) */
            /* so the K2 - K1 is the column count of the block X( K, L ) */
            k1 = iwork[k];
            k2 = iwork[k + 1];
            for (l = nbb;
                    l >= 1;
                    --l)
            {
                /* L1: column index of the first column in X( K, L ) */
                /* L2: column index of the first column in X( K, L + 1) */
                /* so that L2 - L1 is the row count of the block X( K, L ) */
                l1 = iwork[pc + l];
                l2 = iwork[pc + l + 1];
                i__1 = k2 - k1;
                i__2 = l2 - l1;
                dtrsyl_(trana, tranb, isgn, &i__1, &i__2, &a[k1 + k1 * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, &c__[k1 + l1 * c_dim1], ldc, &scaloc, &iinfo);
                *info = fla_max(*info,iinfo);
                if (scaloc * swork[k + l * swork_dim1] == 0.)
                {
                    if (scaloc == 0.)
                    {
                        /* The magnitude of the largest entry of X(K1:K2-1, L1:L2-1) */
                        /* is larger than the product of BIGNUM**2 and cannot be */
                        /* represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1). */
                        /* Mark the computation as pointless. */
                        buf = 0.;
                    }
                    else
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                    }
                    i__1 = nbb;
                    for (jj = 1;
                            jj <= i__1;
                            ++jj)
                    {
                        i__2 = nba;
                        for (ll = 1;
                                ll <= i__2;
                                ++ll)
                        {
                            /* Bound by BIGNUM to not introduce Inf. The value */
                            /* is irrelevant;
                            corresponding entries of the */
                            /* solution will be flushed in consistency scaling. */
                            /* Computing fla_min */
                            frexp(scaloc, (int *) &temp); d__3 = temp;
                            d__1 = bignum;
                            d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                            swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                        }
                    }
                }
                swork[k + l * swork_dim1] = scaloc * swork[k + l * swork_dim1] ;
                i__1 = k2 - k1;
                i__2 = l2 - l1;
                xnrm = dlange_("I", &i__1, &i__2, &c__[k1 + l1 * c_dim1], ldc, wnrm);
                i__1 = k - 1;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    /* C( I, L ) := C( I, L ) - A( I, K ) * C( K, L ) */
                    i1 = iwork[i__];
                    i2 = iwork[i__ + 1];
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__2 = i2 - i1;
                    i__3 = l2 - l1;
                    cnrm = dlange_("I", &i__2, &i__3, &c__[i1 + l1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[i__ + l * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1,d__2);
                    cnrm *= scamin / swork[i__ + l * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    anrm = swork[i__ + (awrk + k) * swork_dim1];
                    scaloc = dlarmm_(&anrm, &xnrm, &cnrm);
                    if (scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                        i__2 = nbb;
                        for (jj = 1;
                                jj <= i__2;
                                ++jj)
                        {
                            i__3 = nba;
                            for (ll = 1;
                                    ll <= i__3;
                                    ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *) &temp); d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                            }
                        }
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        scamin /= pow_dd(&c_b19, &d__1);
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        scaloc /= pow_dd(&c_b19, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( I, L ) and C( K, L ). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__2 = l2 - 1;
                        for (ll = l1;
                                ll <= i__2;
                                ++ll)
                        {
                            i__3 = k2 - k1;
                            dscal_(&i__3, &scal, &c__[k1 + ll * c_dim1], & c__1);
                        }
                    }
                    scal = scamin / swork[i__ + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__2 = l2 - 1;
                        for (ll = l1;
                                ll <= i__2;
                                ++ll)
                        {
                            i__3 = i2 - i1;
                            dscal_(&i__3, &scal, &c__[i1 + ll * c_dim1], & c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[i__ + l * swork_dim1] = scamin * scaloc;
                    i__2 = i2 - i1;
                    i__3 = l2 - l1;
                    i__4 = k2 - k1;
                    dgemm_("N", "N", &i__2, &i__3, &i__4, &c_b31, &a[i1 + k1 * a_dim1], lda, &c__[k1 + l1 * c_dim1], ldc, & c_b32, &c__[i1 + l1 * c_dim1], ldc);
                }
                i__1 = l - 1;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    /* C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**T */
                    j1 = iwork[pc + j];
                    j2 = iwork[pc + j + 1];
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__2 = k2 - k1;
                    i__3 = j2 - j1;
                    cnrm = dlange_("I", &i__2, &i__3, &c__[k1 + j1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[k + j * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1,d__2);
                    cnrm *= scamin / swork[k + j * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    bnrm = swork[l + (bwrk + j) * swork_dim1];
                    scaloc = dlarmm_(&bnrm, &xnrm, &cnrm);
                    if (scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        buf *= pow_dd(&c_b19, &d__1);
                        i__2 = nbb;
                        for (jj = 1;
                                jj <= i__2;
                                ++jj)
                        {
                            i__3 = nba;
                            for (ll = 1;
                                    ll <= i__3;
                                    ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *) &temp); d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1] / pow_dd(&c_b19, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1,d__2);
                            }
                        }
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        
                        
                        scamin /= pow_dd(&c_b19, &d__1);
                        frexp(scaloc, (int *) &temp); d__1 = temp;
                        scaloc /= pow_dd(&c_b19, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( K, J ) and C( K, L ). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__2 = l2 - 1;
                        for (jj = l1;
                                jj <= i__2;
                                ++jj)
                        {
                            i__3 = k2 - k1;
                            dscal_(&i__3, &scal, &c__[k1 + jj * c_dim1], & c__1);
                        }
                    }
                    scal = scamin / swork[k + j * swork_dim1] * scaloc;
                    if (scal != 1.)
                    {
                        i__2 = j2 - 1;
                        for (jj = j1;
                                jj <= i__2;
                                ++jj)
                        {
                            i__3 = k2 - k1;
                            dscal_(&i__3, &scal, &c__[k1 + jj * c_dim1], & c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[k + j * swork_dim1] = scamin * scaloc;
                    i__2 = k2 - k1;
                    i__3 = j2 - j1;
                    i__4 = l2 - l1;
                    d__1 = -sgn;
                    dgemm_("N", "T", &i__2, &i__3, &i__4, &d__1, &c__[k1 + l1 * c_dim1], ldc, &b[j1 + l1 * b_dim1], ldb, &c_b32, &c__[k1 + j1 * c_dim1], ldc);
                }
            }
        }
    }
    /* Reduce local scaling factors */
    *scale = swork[swork_dim1 + 1];
    i__1 = nba;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        i__2 = nbb;
        for (l = 1;
                l <= i__2;
                ++l)
        {
            /* Computing fla_min */
            d__1 = *scale;
            d__2 = swork[k + l * swork_dim1]; // , expr subst
            *scale = fla_min(d__1,d__2);
        }
    }
    if (*scale == 0.)
    {
        /* The magnitude of the largest entry of the solution is larger */
        /* than the product of BIGNUM**2 and cannot be represented in the */
        /* form (1/SCALE)*X if SCALE is DOUBLE PRECISION. Set SCALE to */
        /* zero and give up. */
        iwork[1] = nba + nbb + 2;
        swork[swork_dim1 + 1] = (doublereal) fla_max(nba,nbb);
        swork[swork_dim1 + 2] = (doublereal) ((nbb << 1) + nba);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Realize consistent scaling */
    i__1 = nba;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        k1 = iwork[k];
        k2 = iwork[k + 1];
        i__2 = nbb;
        for (l = 1;
                l <= i__2;
                ++l)
        {
            l1 = iwork[pc + l];
            l2 = iwork[pc + l + 1];
            scal = *scale / swork[k + l * swork_dim1];
            if (scal != 1.)
            {
                i__3 = l2 - 1;
                for (ll = l1;
                        ll <= i__3;
                        ++ll)
                {
                    i__4 = k2 - k1;
                    dscal_(&i__4, &scal, &c__[k1 + ll * c_dim1], &c__1);
                }
            }
        }
    }
    if (buf != 1. && buf > 0.)
    {
        /* Decrease SCALE as much as possible. */
        /* Computing fla_min */
        d__1 = *scale / smlnum;
        d__2 = 1. / buf; // , expr subst
        scaloc = fla_min(d__1,d__2);
        buf *= scaloc;
        *scale /= scaloc;
    }
    if (buf != 1. && buf > 0.)
    {
        /* In case of overly aggressive scaling during the computation, */
        /* flushing of the global scale factor may be prevented by */
        /* undoing some of the scaling. This step is to ensure that */
        /* this routine flushes only scale factors that TRSYL also */
        /* flushes and be usable as a drop-in replacement. */
        /* How much can the normwise largest entry be upscaled? */
        scal = c__[c_dim1 + 1];
        i__1 = *m;
        for (k = 1;
                k <= i__1;
                ++k)
        {
            i__2 = *n;
            for (l = 1;
                    l <= i__2;
                    ++l)
            {
                /* Computing fla_max */
                d__2 = scal;
                d__3 = (d__1 = c__[k + l * c_dim1], f2c_abs(d__1)); // , expr subst
                scal = fla_max(d__2,d__3);
            }
        }
        /* Increase BUF as close to 1 as possible and apply scaling. */
        /* Computing fla_min */
        d__1 = bignum / scal;
        d__2 = 1. / buf; // , expr subst
        scaloc = fla_min(d__1,d__2);
        buf *= scaloc;
        dlascl_("G", &c_n1, &c_n1, &c_b32, &scaloc, m, n, &c__[c_offset], ldc, &iwork[1]);
    }
    /* Combine with buffer scaling factor. SCALE will be flushed if */
    /* BUF is less than one here. */
    *scale *= buf;
    /* Restore workspace dimensions */
    iwork[1] = nba + nbb + 2;
    swork[swork_dim1 + 1] = (doublereal) fla_max(nba,nbb);
    swork[swork_dim1 + 2] = (doublereal) ((nbb << 1) + nba);
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of DTRSYL3 */
}
/* dtrsyl3_ */
