/* ../netlib/sgesvj.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b17 = 0.f;
static real c_b18 = 1.f;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;
/* > \brief \b SGESVJ */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGESVJ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesvj. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesvj. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesvj. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V, */
/* LDV, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDV, LWORK, M, MV, N */
/* CHARACTER*1 JOBA, JOBU, JOBV */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), SVA( N ), V( LDV, * ), */
/* $ WORK( LWORK ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGESVJ computes the singular value decomposition (SVD) of a real */
/* > M-by-N matrix A, where M >= N. The SVD of A is written as */
/* > [++] [xx] [x0] [xx] */
/* > A = U * SIGMA * V^t, [++] = [xx] * [ox] * [xx] */
/* > [++] [xx] */
/* > where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal */
/* > matrix, and V is an N-by-N orthogonal matrix. The diagonal elements */
/* > of SIGMA are the singular values of A. The columns of U and V are the */
/* > left and the right singular vectors of A, respectively. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBA */
/* > \verbatim */
/* > JOBA is CHARACTER* 1 */
/* > Specifies the structure of A. */
/* > = 'L': The input matrix A is lower triangular;
*/
/* > = 'U': The input matrix A is upper triangular;
*/
/* > = 'G': The input matrix A is general M-by-N matrix, M >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBU */
/* > \verbatim */
/* > JOBU is CHARACTER*1 */
/* > Specifies whether to compute the left singular vectors */
/* > (columns of U): */
/* > = 'U': The left singular vectors corresponding to the nonzero */
/* > singular values are computed and returned in the leading */
/* > columns of A. See more details in the description of A. */
/* > The default numerical orthogonality threshold is set to */
/* > approximately TOL=CTOL*EPS, CTOL=SQRT(M), EPS=SLAMCH('E'). */
/* > = 'C': Analogous to JOBU='U', except that user can control the */
/* > level of numerical orthogonality of the computed left */
/* > singular vectors. TOL can be set to TOL = CTOL*EPS, where */
/* > CTOL is given on input in the array WORK. */
/* > No CTOL smaller than ONE is allowed. CTOL greater */
/* > than 1 / EPS is meaningless. The option 'C' */
/* > can be used if M*EPS is satisfactory orthogonality */
/* > of the computed left singular vectors, so CTOL=M could */
/* > save few sweeps of Jacobi rotations. */
/* > See the descriptions of A and WORK(1). */
/* > = 'N': The matrix U is not computed. However, see the */
/* > description of A. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* > JOBV is CHARACTER*1 */
/* > Specifies whether to compute the right singular vectors, that */
/* > is, the matrix V: */
/* > = 'V' : the matrix V is computed and returned in the array V */
/* > = 'A' : the Jacobi rotations are applied to the MV-by-N */
/* > array V. In other words, the right singular vector */
/* > matrix V is not computed explicitly;
instead it is */
/* > applied to an MV-by-N matrix initially stored in the */
/* > first MV rows of V. */
/* > = 'N' : the matrix V is not computed and the array V is not */
/* > referenced */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the input matrix A. 1/SLAMCH('E') > M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the input matrix A. */
/* > M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, */
/* > If JOBU .EQ. 'U' .OR. JOBU .EQ. 'C': */
/* > If INFO .EQ. 0 : */
/* > RANKA orthonormal columns of U are returned in the */
/* > leading RANKA columns of the array A. Here RANKA <= N */
/* > is the number of computed singular values of A that are */
/* > above the underflow threshold SLAMCH('S'). The singular */
/* > vectors corresponding to underflowed or zero singular */
/* > values are not computed. The value of RANKA is returned */
/* > in the array WORK as RANKA=NINT(WORK(2)). Also see the */
/* > descriptions of SVA and WORK. The computed columns of U */
/* > are mutually numerically orthogonal up to approximately */
/* > TOL=SQRT(M)*EPS (default);
or TOL=CTOL*EPS (JOBU.EQ.'C'), */
/* > see the description of JOBU. */
/* > If INFO .GT. 0, */
/* > the procedure SGESVJ did not converge in the given number */
/* > of iterations (sweeps). In that case, the computed */
/* > columns of U may not be orthogonal up to TOL. The output */
/* > U (stored in A), SIGMA (given by the computed singular */
/* > values in SVA(1:N)) and V is still a decomposition of the */
/* > input matrix A in the sense that the residual */
/* > ||A-SCALE*U*SIGMA*V^T||_2 / ||A||_2 is small. */
/* > If JOBU .EQ. 'N': */
/* > If INFO .EQ. 0 : */
/* > Note that the left singular vectors are 'for free' in the */
/* > one-sided Jacobi SVD algorithm. However, if only the */
/* > singular values are needed, the level of numerical */
/* > orthogonality of U is not an issue and iterations are */
/* > stopped when the columns of the iterated matrix are */
/* > numerically orthogonal up to approximately M*EPS. Thus, */
/* > on exit, A contains the columns of U scaled with the */
/* > corresponding singular values. */
/* > If INFO .GT. 0 : */
/* > the procedure SGESVJ did not converge in the given number */
/* > of iterations (sweeps). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] SVA */
/* > \verbatim */
/* > SVA is REAL array, dimension (N) */
/* > On exit, */
/* > If INFO .EQ. 0 : */
/* > depending on the value SCALE = WORK(1), we have: */
/* > If SCALE .EQ. ONE: */
/* > SVA(1:N) contains the computed singular values of A. */
/* > During the computation SVA contains the Euclidean column */
/* > norms of the iterated matrices in the array A. */
/* > If SCALE .NE. ONE: */
/* > The singular values of A are SCALE*SVA(1:N), and this */
/* > factored representation is due to the fact that some of the */
/* > singular values of A might underflow or overflow. */
/* > */
/* > If INFO .GT. 0 : */
/* > the procedure SGESVJ did not converge in the given number of */
/* > iterations (sweeps) and SCALE*SVA(1:N) may not be accurate. */
/* > \endverbatim */
/* > */
/* > \param[in] MV */
/* > \verbatim */
/* > MV is INTEGER */
/* > If JOBV .EQ. 'A', then the product of Jacobi rotations in SGESVJ */
/* > is applied to the first MV rows of V. See the description of JOBV. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* > V is REAL array, dimension (LDV,N) */
/* > If JOBV = 'V', then V contains on exit the N-by-N matrix of */
/* > the right singular vectors;
*/
/* > If JOBV = 'A', then V contains the product of the computed right */
/* > singular vector matrix and the initial matrix in */
/* > the array V. */
/* > If JOBV = 'N', then V is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V, LDV .GE. 1. */
/* > If JOBV .EQ. 'V', then LDV .GE. max(1,N). */
/* > If JOBV .EQ. 'A', then LDV .GE. max(1,MV) . */
/* > \endverbatim */
/* > */
/* > \param[in,out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension max(4,M+N). */
/* > On entry, */
/* > If JOBU .EQ. 'C' : */
/* > WORK(1) = CTOL, where CTOL defines the threshold for convergence. */
/* > The process stops if all columns of A are mutually */
/* > orthogonal up to CTOL*EPS, EPS=SLAMCH('E'). */
/* > It is required that CTOL >= ONE, i.e. it is not */
/* > allowed to force the routine to obtain orthogonality */
/* > below EPSILON. */
/* > On exit, */
/* > WORK(1) = SCALE is the scaling factor such that SCALE*SVA(1:N) */
/* > are the computed singular vcalues of A. */
/* > (See description of SVA().) */
/* > WORK(2) = NINT(WORK(2)) is the number of the computed nonzero */
/* > singular values. */
/* > WORK(3) = NINT(WORK(3)) is the number of the computed singular */
/* > values that are larger than the underflow threshold. */
/* > WORK(4) = NINT(WORK(4)) is the number of sweeps of Jacobi */
/* > rotations needed for numerical convergence. */
/* > WORK(5) = max_{
i.NE.j}
|COS(A(:,i),A(:,j))| in the last sweep. */
/* > This is useful information in cases when SGESVJ did */
/* > not converge, as it can be used to estimate whether */
/* > the output is stil useful and for post festum analysis. */
/* > WORK(6) = the largest absolute value over all sines of the */
/* > Jacobi rotation angles in the last sweep. It can be */
/* > useful for a post festum analysis. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > length of WORK, WORK >= MAX(6,M+N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0 : successful exit. */
/* > < 0 : if INFO = -i, then the i-th argument had an illegal value */
/* > > 0 : SGESVJ did not converge in the maximal allowed number (30) */
/* > of sweeps. The output may still be useful. See the */
/* > description of WORK. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > The orthogonal N-by-N matrix V is obtained as a product of Jacobi plane */
/* > rotations. The rotations are implemented as fast scaled rotations of */
/* > Anda and Park [1]. In the case of underflow of the Jacobi angle, a */
/* > modified Jacobi transformation of Drmac [4] is used. Pivot strategy uses */
/* > column interchanges of de Rijk [2]. The relative accuracy of the computed */
/* > singular values and the accuracy of the computed singular vectors (in */
/* > angle metric) is as guaranteed by the theory of Demmel and Veselic [3]. */
/* > The condition number that determines the accuracy in the full rank case */
/* > is essentially min_{
D=diag}
kappa(A*D), where kappa(.) is the */
/* > spectral condition number. The best performance of this Jacobi SVD */
/* > procedure is achieved if used in an accelerated version of Drmac and */
/* > Veselic [5,6], and it is the kernel routine in the SIGMA library [7]. */
/* > Some tunning parameters (marked with [TP]) are available for the */
/* > implementer. \n */
/* > The computational range for the nonzero singular values is the machine */
/* > number interval ( UNDERFLOW , OVERFLOW ). In extreme cases, even */
/* > denormalized singular values can be computed with the corresponding */
/* > gradual loss of accurate digits. */
/* > */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany) */
/* > */
/* > \par References: */
/* ================ */
/* > */
/* > [1] A. A. Anda and H. Park: Fast plane rotations with dynamic scaling. \n */
/* > SIAM J. matrix Anal. Appl., Vol. 15 (1994), pp. 162-174. \n\n */
/* > [2] P. P. M. De Rijk: A one-sided Jacobi algorithm for computing the */
/* > singular value decomposition on a vector computer. \n */
/* > SIAM J. Sci. Stat. Comp., Vol. 10 (1998), pp. 359-371. \n\n */
/* > [3] J. Demmel and K. Veselic: Jacobi method is more accurate than QR. \n */
/* > [4] Z. Drmac: Implementation of Jacobi rotations for accurate singular */
/* > value computation in floating point arithmetic. \n */
/* > SIAM J. Sci. Comp., Vol. 18 (1997), pp. 1200-1222. \n\n */
/* > [5] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I. \n */
/* > SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342. \n */
/* > LAPACK Working note 169. \n\n */
/* > [6] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II. \n */
/* > SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362. \n */
/* > LAPACK Working note 170. \n\n */
/* > [7] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV, */
/* > QSVD, (H,K)-SVD computations.\n */
/* > Department of Mathematics, University of Zagreb, 2008. */
/* > */
/* > \par Bugs, Examples and Comments: */
/* ================================= */
/* > */
/* > Please report all bugs and send interesting test examples and comments to */
/* > drmac@math.hr. Thank you. */
/* ===================================================================== */
/* Subroutine */
int sgesvj_(char *joba, char *jobu, char *jobv, integer *m, integer *n, real *a, integer *lda, real *sva, integer *mv, real *v, integer *ldv, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);
    /* Local variables */
    real bigtheta;
    integer pskipped, i__, p, q;
    real t;
    integer n2, n4;
    real rootsfmin;
    integer n34;
    real cs, sn;
    integer ir1, jbc;
    real big;
    integer kbl, igl, ibr, jgl, nbl;
    real skl, tol;
    integer mvl;
    real aapp, aapq, aaqq, ctol;
    integer ierr;
    extern real sdot_(integer *, real *, integer *, real *, integer *);
    real aapp0, temp1;
    extern real snrm2_(integer *, real *, integer *);
    real large, apoaq, aqoap;
    extern logical lsame_(char *, char *);
    real theta;
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    real small, sfmin;
    logical lsvec;
    real fastr[5], epsln;
    logical applv, rsvec, uctol, lower, upper;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *);
    logical rotok;
    extern /* Subroutine */
    int sswap_(integer *, real *, integer *, real *, integer *), saxpy_(integer *, real *, real *, integer *, real *, integer *), srotm_(integer *, real *, integer *, real *, integer * , real *), sgsvj0_(char *, integer *, integer *, real *, integer * , real *, real *, integer *, real *, integer *, real *, real *, real *, integer *, real *, integer *, integer *), sgsvj1_( char *, integer *, integer *, integer *, real *, integer *, real * , real *, integer *, real *, integer *, real *, real *, real *, integer *, real *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    integer ijblsk, swband;
    extern /* Subroutine */
    int slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *, integer *, integer *);
    extern integer isamax_(integer *, real *, integer *);
    integer blskip;
    real mxaapq;
    extern /* Subroutine */
    int slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
    real thsign;
    extern /* Subroutine */
    int slassq_(integer *, real *, integer *, real *, real *);
    real mxsinj;
    integer emptsw, notrot, iswrot, lkahead;
    logical goscale, noscale;
    real rootbig, rooteps;
    integer rowskip;
    real roottol;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Local Arrays .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* from BLAS */
    /* from LAPACK */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* from BLAS */
    /* from LAPACK */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    --sva;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --work;
    /* Function Body */
    lsvec = lsame_(jobu, "U");
    uctol = lsame_(jobu, "C");
    rsvec = lsame_(jobv, "V");
    applv = lsame_(jobv, "A");
    upper = lsame_(joba, "U");
    lower = lsame_(joba, "L");
    if (! (upper || lower || lsame_(joba, "G")))
    {
        *info = -1;
    }
    else if (! (lsvec || uctol || lsame_(jobu, "N")))
    {
        *info = -2;
    }
    else if (! (rsvec || applv || lsame_(jobv, "N")))
    {
        *info = -3;
    }
    else if (*m < 0)
    {
        *info = -4;
    }
    else if (*n < 0 || *n > *m)
    {
        *info = -5;
    }
    else if (*lda < *m)
    {
        *info = -7;
    }
    else if (*mv < 0)
    {
        *info = -9;
    }
    else if (rsvec && *ldv < *n || applv && *ldv < *mv)
    {
        *info = -11;
    }
    else if (uctol && work[1] <= 1.f)
    {
        *info = -12;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = *m + *n;
        if (*lwork < max(i__1,6))
        {
            *info = -13;
        }
        else
        {
            *info = 0;
        }
    }
    /* #:( */
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGESVJ", &i__1);
        return 0;
    }
    /* #:) Quick return for void matrix */
    if (*m == 0 || *n == 0)
    {
        return 0;
    }
    /* Set numerical parameters */
    /* The stopping criterion for Jacobi rotations is */
    /* max_{
    i<>j}
    |A(:,i)^T * A(:,j)|/(||A(:,i)||*||A(:,j)||) < CTOL*EPS */
    /* where EPS is the round-off and CTOL is defined as follows: */
    if (uctol)
    {
        /* ... user controlled */
        ctol = work[1];
    }
    else
    {
        /* ... default */
        if (lsvec || rsvec || applv)
        {
            ctol = sqrt((real) (*m));
        }
        else
        {
            ctol = (real) (*m);
        }
    }
    /* ... and the machine dependent parameters are */
    /* [!] (Make sure that SLAMCH() works properly on the target machine.) */
    epsln = slamch_("Epsilon");
    rooteps = sqrt(epsln);
    sfmin = slamch_("SafeMinimum");
    rootsfmin = sqrt(sfmin);
    small = sfmin / epsln;
    big = slamch_("Overflow");
    /* BIG = ONE / SFMIN */
    rootbig = 1.f / rootsfmin;
    large = big / sqrt((real) (*m * *n));
    bigtheta = 1.f / rooteps;
    tol = ctol * epsln;
    roottol = sqrt(tol);
    if ((real) (*m) * epsln >= 1.f)
    {
        *info = -4;
        i__1 = -(*info);
        xerbla_("SGESVJ", &i__1);
        return 0;
    }
    /* Initialize the right singular vector matrix. */
    if (rsvec)
    {
        mvl = *n;
        slaset_("A", &mvl, n, &c_b17, &c_b18, &v[v_offset], ldv);
    }
    else if (applv)
    {
        mvl = *mv;
    }
    rsvec = rsvec || applv;
    /* Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N ) */
    /* (!) If necessary, scale A to protect the largest singular value */
    /* from overflow. It is possible that saving the largest singular */
    /* value destroys the information about the small ones. */
    /* This initial scaling is almost minimal in the sense that the */
    /* goal is to make sure that no column norm overflows, and that */
    /* SQRT(N)*max_i SVA(i) does not overflow. If INFinite entries */
    /* in A are detected, the procedure returns with INFO=-6. */
    skl = 1.f / sqrt((real) (*m) * (real) (*n));
    noscale = TRUE_;
    goscale = TRUE_;
    if (lower)
    {
        /* the input matrix is M-by-N lower triangular (trapezoidal) */
        i__1 = *n;
        for (p = 1;
                p <= i__1;
                ++p)
        {
            aapp = 0.f;
            aaqq = 1.f;
            i__2 = *m - p + 1;
            slassq_(&i__2, &a[p + p * a_dim1], &c__1, &aapp, &aaqq);
            if (aapp > big)
            {
                *info = -6;
                i__2 = -(*info);
                xerbla_("SGESVJ", &i__2);
                return 0;
            }
            aaqq = sqrt(aaqq);
            if (aapp < big / aaqq && noscale)
            {
                sva[p] = aapp * aaqq;
            }
            else
            {
                noscale = FALSE_;
                sva[p] = aapp * (aaqq * skl);
                if (goscale)
                {
                    goscale = FALSE_;
                    i__2 = p - 1;
                    for (q = 1;
                            q <= i__2;
                            ++q)
                    {
                        sva[q] *= skl;
                        /* L1873: */
                    }
                }
            }
            /* L1874: */
        }
    }
    else if (upper)
    {
        /* the input matrix is M-by-N upper triangular (trapezoidal) */
        i__1 = *n;
        for (p = 1;
                p <= i__1;
                ++p)
        {
            aapp = 0.f;
            aaqq = 1.f;
            slassq_(&p, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
            if (aapp > big)
            {
                *info = -6;
                i__2 = -(*info);
                xerbla_("SGESVJ", &i__2);
                return 0;
            }
            aaqq = sqrt(aaqq);
            if (aapp < big / aaqq && noscale)
            {
                sva[p] = aapp * aaqq;
            }
            else
            {
                noscale = FALSE_;
                sva[p] = aapp * (aaqq * skl);
                if (goscale)
                {
                    goscale = FALSE_;
                    i__2 = p - 1;
                    for (q = 1;
                            q <= i__2;
                            ++q)
                    {
                        sva[q] *= skl;
                        /* L2873: */
                    }
                }
            }
            /* L2874: */
        }
    }
    else
    {
        /* the input matrix is M-by-N general dense */
        i__1 = *n;
        for (p = 1;
                p <= i__1;
                ++p)
        {
            aapp = 0.f;
            aaqq = 1.f;
            slassq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
            if (aapp > big)
            {
                *info = -6;
                i__2 = -(*info);
                xerbla_("SGESVJ", &i__2);
                return 0;
            }
            aaqq = sqrt(aaqq);
            if (aapp < big / aaqq && noscale)
            {
                sva[p] = aapp * aaqq;
            }
            else
            {
                noscale = FALSE_;
                sva[p] = aapp * (aaqq * skl);
                if (goscale)
                {
                    goscale = FALSE_;
                    i__2 = p - 1;
                    for (q = 1;
                            q <= i__2;
                            ++q)
                    {
                        sva[q] *= skl;
                        /* L3873: */
                    }
                }
            }
            /* L3874: */
        }
    }
    if (noscale)
    {
        skl = 1.f;
    }
    /* Move the smaller part of the spectrum from the underflow threshold */
    /* (!) Start by determining the position of the nonzero entries of the */
    /* array SVA() relative to ( SFMIN, BIG ). */
    aapp = 0.f;
    aaqq = big;
    i__1 = *n;
    for (p = 1;
            p <= i__1;
            ++p)
    {
        if (sva[p] != 0.f)
        {
            /* Computing MIN */
            r__1 = aaqq;
            r__2 = sva[p]; // , expr subst
            aaqq = min(r__1,r__2);
        }
        /* Computing MAX */
        r__1 = aapp;
        r__2 = sva[p]; // , expr subst
        aapp = max(r__1,r__2);
        /* L4781: */
    }
    /* #:) Quick return for zero matrix */
    if (aapp == 0.f)
    {
        if (lsvec)
        {
            slaset_("G", m, n, &c_b17, &c_b18, &a[a_offset], lda);
        }
        work[1] = 1.f;
        work[2] = 0.f;
        work[3] = 0.f;
        work[4] = 0.f;
        work[5] = 0.f;
        work[6] = 0.f;
        return 0;
    }
    /* #:) Quick return for one-column matrix */
    if (*n == 1)
    {
        if (lsvec)
        {
            slascl_("G", &c__0, &c__0, &sva[1], &skl, m, &c__1, &a[a_dim1 + 1] , lda, &ierr);
        }
        work[1] = 1.f / skl;
        if (sva[1] >= sfmin)
        {
            work[2] = 1.f;
        }
        else
        {
            work[2] = 0.f;
        }
        work[3] = 0.f;
        work[4] = 0.f;
        work[5] = 0.f;
        work[6] = 0.f;
        return 0;
    }
    /* Protect small singular values from underflow, and try to */
    /* avoid underflows/overflows in computing Jacobi rotations. */
    sn = sqrt(sfmin / epsln);
    temp1 = sqrt(big / (real) (*n));
    if (aapp <= sn || aaqq >= temp1 || sn <= aaqq && aapp <= temp1)
    {
        /* Computing MIN */
        r__1 = big;
        r__2 = temp1 / aapp; // , expr subst
        temp1 = min(r__1,r__2);
        /* AAQQ = AAQQ*TEMP1 */
        /* AAPP = AAPP*TEMP1 */
    }
    else if (aaqq <= sn && aapp <= temp1)
    {
        /* Computing MIN */
        r__1 = sn / aaqq;
        r__2 = big / (aapp * sqrt((real) (*n))); // , expr subst
        temp1 = min(r__1,r__2);
        /* AAQQ = AAQQ*TEMP1 */
        /* AAPP = AAPP*TEMP1 */
    }
    else if (aaqq >= sn && aapp >= temp1)
    {
        /* Computing MAX */
        r__1 = sn / aaqq;
        r__2 = temp1 / aapp; // , expr subst
        temp1 = max(r__1,r__2);
        /* AAQQ = AAQQ*TEMP1 */
        /* AAPP = AAPP*TEMP1 */
    }
    else if (aaqq <= sn && aapp >= temp1)
    {
        /* Computing MIN */
        r__1 = sn / aaqq;
        r__2 = big / (sqrt((real) (*n)) * aapp); // , expr subst
        temp1 = min(r__1,r__2);
        /* AAQQ = AAQQ*TEMP1 */
        /* AAPP = AAPP*TEMP1 */
    }
    else
    {
        temp1 = 1.f;
    }
    /* Scale, if necessary */
    if (temp1 != 1.f)
    {
        slascl_("G", &c__0, &c__0, &c_b18, &temp1, n, &c__1, &sva[1], n, & ierr);
    }
    skl = temp1 * skl;
    if (skl != 1.f)
    {
        slascl_(joba, &c__0, &c__0, &c_b18, &skl, m, n, &a[a_offset], lda, & ierr);
        skl = 1.f / skl;
    }
    /* Row-cyclic Jacobi SVD algorithm with column pivoting */
    emptsw = *n * (*n - 1) / 2;
    notrot = 0;
    fastr[0] = 0.f;
    /* A is represented in factored form A = A * diag(WORK), where diag(WORK) */
    /* is initialized to identity. WORK is updated during fast scaled */
    /* rotations. */
    i__1 = *n;
    for (q = 1;
            q <= i__1;
            ++q)
    {
        work[q] = 1.f;
        /* L1868: */
    }
    swband = 3;
    /* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
    /* if SGESVJ is used as a computational routine in the preconditioned */
    /* Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure */
    /* works on pivots inside a band-like region around the diagonal. */
    /* The boundaries are determined dynamically, based on the number of */
    /* pivots above a threshold. */
    kbl = min(8,*n);
    /* [TP] KBL is a tuning parameter that defines the tile size in the */
    /* tiling of the p-q loops of pivot pairs. In general, an optimal */
    /* value of KBL depends on the matrix dimensions and on the */
    /* parameters of the computer's memory. */
    nbl = *n / kbl;
    if (nbl * kbl != *n)
    {
        ++nbl;
    }
    /* Computing 2nd power */
    i__1 = kbl;
    blskip = i__1 * i__1;
    /* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */
    rowskip = min(5,kbl);
    /* [TP] ROWSKIP is a tuning parameter. */
    lkahead = 1;
    /* [TP] LKAHEAD is a tuning parameter. */
    /* Quasi block transformations, using the lower (upper) triangular */
    /* structure of the input matrix. The quasi-block-cycling usually */
    /* invokes cubic convergence. Big part of this cycle is done inside */
    /* canonical subspaces of dimensions less than M. */
    /* Computing MAX */
    i__1 = 64;
    i__2 = kbl << 2; // , expr subst
    if ((lower || upper) && *n > max(i__1,i__2))
    {
        /* [TP] The number of partition levels and the actual partition are */
        /* tuning parameters. */
        n4 = *n / 4;
        n2 = *n / 2;
        n34 = n4 * 3;
        if (applv)
        {
            q = 0;
        }
        else
        {
            q = 1;
        }
        if (lower)
        {
            /* This works very well on lower triangular matrices, in particular */
            /* in the framework of the preconditioned Jacobi SVD (xGEJSV). */
            /* The idea is simple: */
            /* [+ 0 0 0] Note that Jacobi transformations of [0 0] */
            /* [+ + 0 0] [0 0] */
            /* [+ + x 0] actually work on [x 0] [x 0] */
            /* [+ + x x] [x x]. [x x] */
            i__1 = *m - n34;
            i__2 = *n - n34;
            i__3 = *lwork - *n;
            sgsvj0_(jobv, &i__1, &i__2, &a[n34 + 1 + (n34 + 1) * a_dim1], lda, &work[n34 + 1], &sva[n34 + 1], &mvl, &v[n34 * q + 1 + ( n34 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, & work[*n + 1], &i__3, &ierr);
            i__1 = *m - n2;
            i__2 = n34 - n2;
            i__3 = *lwork - *n;
            sgsvj0_(jobv, &i__1, &i__2, &a[n2 + 1 + (n2 + 1) * a_dim1], lda, & work[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &work[*n + 1], &i__3, &ierr);
            i__1 = *m - n2;
            i__2 = *n - n2;
            i__3 = *lwork - *n;
            sgsvj1_(jobv, &i__1, &i__2, &n4, &a[n2 + 1 + (n2 + 1) * a_dim1], lda, &work[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + ( n2 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, & work[*n + 1], &i__3, &ierr);
            i__1 = *m - n4;
            i__2 = n2 - n4;
            i__3 = *lwork - *n;
            sgsvj0_(jobv, &i__1, &i__2, &a[n4 + 1 + (n4 + 1) * a_dim1], lda, & work[n4 + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 1], &i__3, &ierr);
            i__1 = *lwork - *n;
            sgsvj0_(jobv, m, &n4, &a[a_offset], lda, &work[1], &sva[1], &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 1], &i__1, &ierr);
            i__1 = *lwork - *n;
            sgsvj1_(jobv, m, &n2, &n4, &a[a_offset], lda, &work[1], &sva[1], & mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, & work[*n + 1], &i__1, &ierr);
        }
        else if (upper)
        {
            i__1 = *lwork - *n;
            sgsvj0_(jobv, &n4, &n4, &a[a_offset], lda, &work[1], &sva[1], & mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__2, & work[*n + 1], &i__1, &ierr);
            i__1 = *lwork - *n;
            sgsvj0_(jobv, &n2, &n4, &a[(n4 + 1) * a_dim1 + 1], lda, &work[n4 + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 1], &i__1, &ierr);
            i__1 = *lwork - *n;
            sgsvj1_(jobv, &n2, &n2, &n4, &a[a_offset], lda, &work[1], &sva[1], &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, & work[*n + 1], &i__1, &ierr);
            i__1 = n2 + n4;
            i__2 = *lwork - *n;
            sgsvj0_(jobv, &i__1, &n4, &a[(n2 + 1) * a_dim1 + 1], lda, &work[ n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 1], &i__2, &ierr);
        }
    }
    /* .. Row-cyclic pivot strategy with de Rijk's pivoting .. */
    for (i__ = 1;
            i__ <= 30;
            ++i__)
    {
        /* .. go go go ... */
        mxaapq = 0.f;
        mxsinj = 0.f;
        iswrot = 0;
        notrot = 0;
        pskipped = 0;
        /* Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
        /* 1 <= p < q <= N. This is the first step toward a blocked implementation */
        /* of the rotations. New implementation, based on block transformations, */
        /* is under development. */
        i__1 = nbl;
        for (ibr = 1;
                ibr <= i__1;
                ++ibr)
        {
            igl = (ibr - 1) * kbl + 1;
            /* Computing MIN */
            i__3 = lkahead;
            i__4 = nbl - ibr; // , expr subst
            i__2 = min(i__3,i__4);
            for (ir1 = 0;
                    ir1 <= i__2;
                    ++ir1)
            {
                igl += ir1 * kbl;
                /* Computing MIN */
                i__4 = igl + kbl - 1;
                i__5 = *n - 1; // , expr subst
                i__3 = min(i__4,i__5);
                for (p = igl;
                        p <= i__3;
                        ++p)
                {
                    /* .. de Rijk's pivoting */
                    i__4 = *n - p + 1;
                    q = isamax_(&i__4, &sva[p], &c__1) + p - 1;
                    if (p != q)
                    {
                        sswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
                        if (rsvec)
                        {
                            sswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &c__1);
                        }
                        temp1 = sva[p];
                        sva[p] = sva[q];
                        sva[q] = temp1;
                        temp1 = work[p];
                        work[p] = work[q];
                        work[q] = temp1;
                    }
                    if (ir1 == 0)
                    {
                        /* Column norms are periodically updated by explicit */
                        /* norm computation. */
                        /* Caveat: */
                        /* Unfortunately, some BLAS implementations compute SNRM2(M,A(1,p),1) */
                        /* as SQRT(SDOT(M,A(1,p),1,A(1,p),1)), which may cause the result to */
                        /* overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to */
                        /* underflow for ||A(:,p)||_2 < SQRT(underflow_threshold). */
                        /* Hence, SNRM2 cannot be trusted, not even in the case when */
                        /* the true norm is far from the under(over)flow boundaries. */
                        /* If properly implemented SNRM2 is available, the IF-THEN-ELSE */
                        /* below should read "AAPP = SNRM2( M, A(1,p), 1 ) * WORK(p)". */
                        if (sva[p] < rootbig && sva[p] > rootsfmin)
                        {
                            sva[p] = snrm2_(m, &a[p * a_dim1 + 1], &c__1) * work[p];
                        }
                        else
                        {
                            temp1 = 0.f;
                            aapp = 1.f;
                            slassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, & aapp);
                            sva[p] = temp1 * sqrt(aapp) * work[p];
                        }
                        aapp = sva[p];
                    }
                    else
                    {
                        aapp = sva[p];
                    }
                    if (aapp > 0.f)
                    {
                        pskipped = 0;
                        /* Computing MIN */
                        i__5 = igl + kbl - 1;
                        i__4 = min(i__5,*n);
                        for (q = p + 1;
                                q <= i__4;
                                ++q)
                        {
                            aaqq = sva[q];
                            if (aaqq > 0.f)
                            {
                                aapp0 = aapp;
                                if (aaqq >= 1.f)
                                {
                                    rotok = small * aapp <= aaqq;
                                    if (aapp < big / aaqq)
                                    {
                                        aapq = sdot_(m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], & c__1) * work[p] * work[q] / aaqq / aapp;
                                    }
                                    else
                                    {
                                        scopy_(m, &a[p * a_dim1 + 1], &c__1, & work[*n + 1], &c__1);
                                        slascl_("G", &c__0, &c__0, &aapp, & work[p], m, &c__1, &work[*n + 1], lda, &ierr);
                                        aapq = sdot_(m, &work[*n + 1], &c__1, &a[q * a_dim1 + 1], &c__1) * work[q] / aaqq;
                                    }
                                }
                                else
                                {
                                    rotok = aapp <= aaqq / small;
                                    if (aapp > small / aaqq)
                                    {
                                        aapq = sdot_(m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], & c__1) * work[p] * work[q] / aaqq / aapp;
                                    }
                                    else
                                    {
                                        scopy_(m, &a[q * a_dim1 + 1], &c__1, & work[*n + 1], &c__1);
                                        slascl_("G", &c__0, &c__0, &aaqq, & work[q], m, &c__1, &work[*n + 1], lda, &ierr);
                                        aapq = sdot_(m, &work[*n + 1], &c__1, &a[p * a_dim1 + 1], &c__1) * work[p] / aapp;
                                    }
                                }
                                /* Computing MAX */
                                r__1 = mxaapq;
                                r__2 = f2c_abs(aapq); // , expr subst
                                mxaapq = max(r__1,r__2);
                                /* TO rotate or NOT to rotate, THAT is the question ... */
                                if (f2c_abs(aapq) > tol)
                                {
                                    /* .. rotate */
                                    /* [RTD] ROTATED = ROTATED + ONE */
                                    if (ir1 == 0)
                                    {
                                        notrot = 0;
                                        pskipped = 0;
                                        ++iswrot;
                                    }
                                    if (rotok)
                                    {
                                        aqoap = aaqq / aapp;
                                        apoaq = aapp / aaqq;
                                        theta = (r__1 = aqoap - apoaq, f2c_abs( r__1)) * -.5f / aapq;
                                        if (f2c_abs(theta) > bigtheta)
                                        {
                                            t = .5f / theta;
                                            fastr[2] = t * work[p] / work[q];
                                            fastr[3] = -t * work[q] / work[p];
                                            srotm_(m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], &c__1, fastr);
                                            if (rsvec)
                                            {
                                                srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &c__1, fastr);
                                            }
                                            /* Computing MAX */
                                            r__1 = 0.f;
                                            r__2 = t * apoaq * aapq + 1.f; // , expr subst
                                            sva[q] = aaqq * sqrt((max(r__1, r__2)));
                                            /* Computing MAX */
                                            r__1 = 0.f;
                                            r__2 = 1.f - t * aqoap * aapq; // , expr subst
                                            aapp *= sqrt((max(r__1,r__2)));
                                            /* Computing MAX */
                                            r__1 = mxsinj;
                                            r__2 = f2c_abs(t); // , expr subst
                                            mxsinj = max(r__1,r__2);
                                        }
                                        else
                                        {
                                            /* .. choose correct signum for THETA and rotate */
                                            thsign = -r_sign(&c_b18, &aapq);
                                            t = 1.f / (theta + thsign * sqrt( theta * theta + 1.f));
                                            cs = sqrt(1.f / (t * t + 1.f));
                                            sn = t * cs;
                                            /* Computing MAX */
                                            r__1 = mxsinj;
                                            r__2 = f2c_abs(sn); // , expr subst
                                            mxsinj = max(r__1,r__2);
                                            /* Computing MAX */
                                            r__1 = 0.f;
                                            r__2 = t * apoaq * aapq + 1.f; // , expr subst
                                            sva[q] = aaqq * sqrt((max(r__1, r__2)));
                                            /* Computing MAX */
                                            r__1 = 0.f;
                                            r__2 = 1.f - t * aqoap * aapq; // , expr subst
                                            aapp *= sqrt((max(r__1,r__2)));
                                            apoaq = work[p] / work[q];
                                            aqoap = work[q] / work[p];
                                            if (work[p] >= 1.f)
                                            {
                                                if (work[q] >= 1.f)
                                                {
                                                    fastr[2] = t * apoaq;
                                                    fastr[3] = -t * aqoap;
                                                    work[p] *= cs;
                                                    work[q] *= cs;
                                                    srotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1, fastr);
                                                    if (rsvec)
                                                    {
                                                        srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[ q * v_dim1 + 1], &c__1, fastr);
                                                    }
                                                }
                                                else
                                                {
                                                    r__1 = -t * aqoap;
                                                    saxpy_(m, &r__1, &a[q * a_dim1 + 1], &c__1, &a[ p * a_dim1 + 1], &c__1);
                                                    r__1 = cs * sn * apoaq;
                                                    saxpy_(m, &r__1, &a[p * a_dim1 + 1], &c__1, &a[ q * a_dim1 + 1], &c__1);
                                                    work[p] *= cs;
                                                    work[q] /= cs;
                                                    if (rsvec)
                                                    {
                                                        r__1 = -t * aqoap;
                                                        saxpy_(&mvl, &r__1, &v[q * v_dim1 + 1], & c__1, &v[p * v_dim1 + 1], &c__1);
                                                        r__1 = cs * sn * apoaq;
                                                        saxpy_(&mvl, &r__1, &v[p * v_dim1 + 1], & c__1, &v[q * v_dim1 + 1], &c__1);
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (work[q] >= 1.f)
                                                {
                                                    r__1 = t * apoaq;
                                                    saxpy_(m, &r__1, &a[p * a_dim1 + 1], &c__1, &a[ q * a_dim1 + 1], &c__1);
                                                    r__1 = -cs * sn * aqoap;
                                                    saxpy_(m, &r__1, &a[q * a_dim1 + 1], &c__1, &a[ p * a_dim1 + 1], &c__1);
                                                    work[p] /= cs;
                                                    work[q] *= cs;
                                                    if (rsvec)
                                                    {
                                                        r__1 = t * apoaq;
                                                        saxpy_(&mvl, &r__1, &v[p * v_dim1 + 1], & c__1, &v[q * v_dim1 + 1], &c__1);
                                                        r__1 = -cs * sn * aqoap;
                                                        saxpy_(&mvl, &r__1, &v[q * v_dim1 + 1], & c__1, &v[p * v_dim1 + 1], &c__1);
                                                    }
                                                }
                                                else
                                                {
                                                    if (work[p] >= work[q])
                                                    {
                                                        r__1 = -t * aqoap;
                                                        saxpy_(m, &r__1, &a[q * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &c__1);
                                                        r__1 = cs * sn * apoaq;
                                                        saxpy_(m, &r__1, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
                                                        work[p] *= cs;
                                                        work[q] /= cs;
                                                        if (rsvec)
                                                        {
                                                            r__1 = -t * aqoap;
                                                            saxpy_(&mvl, &r__1, &v[q * v_dim1 + 1], &c__1, &v[p * v_dim1 + 1], & c__1);
                                                            r__1 = cs * sn * apoaq;
                                                            saxpy_(&mvl, &r__1, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], & c__1);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        r__1 = t * apoaq;
                                                        saxpy_(m, &r__1, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
                                                        r__1 = -cs * sn * aqoap;
                                                        saxpy_(m, &r__1, &a[q * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &c__1);
                                                        work[p] /= cs;
                                                        work[q] *= cs;
                                                        if (rsvec)
                                                        {
                                                            r__1 = t * apoaq;
                                                            saxpy_(&mvl, &r__1, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], & c__1);
                                                            r__1 = -cs * sn * aqoap;
                                                            saxpy_(&mvl, &r__1, &v[q * v_dim1 + 1], &c__1, &v[p * v_dim1 + 1], & c__1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        /* .. have to use modified Gram-Schmidt like transformation */
                                        scopy_(m, &a[p * a_dim1 + 1], &c__1, & work[*n + 1], &c__1);
                                        slascl_("G", &c__0, &c__0, &aapp, & c_b18, m, &c__1, &work[*n + 1] , lda, &ierr);
                                        slascl_("G", &c__0, &c__0, &aaqq, & c_b18, m, &c__1, &a[q * a_dim1 + 1], lda, &ierr);
                                        temp1 = -aapq * work[p] / work[q];
                                        saxpy_(m, &temp1, &work[*n + 1], & c__1, &a[q * a_dim1 + 1], & c__1);
                                        slascl_("G", &c__0, &c__0, &c_b18, & aaqq, m, &c__1, &a[q * a_dim1 + 1], lda, &ierr);
                                        /* Computing MAX */
                                        r__1 = 0.f;
                                        r__2 = 1.f - aapq * aapq; // , expr subst
                                        sva[q] = aaqq * sqrt((max(r__1,r__2))) ;
                                        mxsinj = max(mxsinj,sfmin);
                                    }
                                    /* END IF ROTOK THEN ... ELSE */
                                    /* In the case of cancellation in updating SVA(q), SVA(p) */
                                    /* recompute SVA(q), SVA(p). */
                                    /* Computing 2nd power */
                                    r__1 = sva[q] / aaqq;
                                    if (r__1 * r__1 <= rooteps)
                                    {
                                        if (aaqq < rootbig && aaqq > rootsfmin)
                                        {
                                            sva[q] = snrm2_(m, &a[q * a_dim1 + 1], &c__1) * work[q];
                                        }
                                        else
                                        {
                                            t = 0.f;
                                            aaqq = 1.f;
                                            slassq_(m, &a[q * a_dim1 + 1], & c__1, &t, &aaqq);
                                            sva[q] = t * sqrt(aaqq) * work[q];
                                        }
                                    }
                                    if (aapp / aapp0 <= rooteps)
                                    {
                                        if (aapp < rootbig && aapp > rootsfmin)
                                        {
                                            aapp = snrm2_(m, &a[p * a_dim1 + 1], &c__1) * work[p];
                                        }
                                        else
                                        {
                                            t = 0.f;
                                            aapp = 1.f;
                                            slassq_(m, &a[p * a_dim1 + 1], & c__1, &t, &aapp);
                                            aapp = t * sqrt(aapp) * work[p];
                                        }
                                        sva[p] = aapp;
                                    }
                                }
                                else
                                {
                                    /* A(:,p) and A(:,q) already numerically orthogonal */
                                    if (ir1 == 0)
                                    {
                                        ++notrot;
                                    }
                                    /* [RTD] SKIPPED = SKIPPED + 1 */
                                    ++pskipped;
                                }
                            }
                            else
                            {
                                /* A(:,q) is zero column */
                                if (ir1 == 0)
                                {
                                    ++notrot;
                                }
                                ++pskipped;
                            }
                            if (i__ <= swband && pskipped > rowskip)
                            {
                                if (ir1 == 0)
                                {
                                    aapp = -aapp;
                                }
                                notrot = 0;
                                goto L2103;
                            }
                            /* L2002: */
                        }
                        /* END q-LOOP */
L2103: /* bailed out of q-loop */
                        sva[p] = aapp;
                    }
                    else
                    {
                        sva[p] = aapp;
                        if (ir1 == 0 && aapp == 0.f)
                        {
                            /* Computing MIN */
                            i__4 = igl + kbl - 1;
                            notrot = notrot + min(i__4,*n) - p;
                        }
                    }
                    /* L2001: */
                }
                /* end of the p-loop */
                /* end of doing the block ( ibr, ibr ) */
                /* L1002: */
            }
            /* end of ir1-loop */
            /* ... go to the off diagonal blocks */
            igl = (ibr - 1) * kbl + 1;
            i__2 = nbl;
            for (jbc = ibr + 1;
                    jbc <= i__2;
                    ++jbc)
            {
                jgl = (jbc - 1) * kbl + 1;
                /* doing the block at ( ibr, jbc ) */
                ijblsk = 0;
                /* Computing MIN */
                i__4 = igl + kbl - 1;
                i__3 = min(i__4,*n);
                for (p = igl;
                        p <= i__3;
                        ++p)
                {
                    aapp = sva[p];
                    if (aapp > 0.f)
                    {
                        pskipped = 0;
                        /* Computing MIN */
                        i__5 = jgl + kbl - 1;
                        i__4 = min(i__5,*n);
                        for (q = jgl;
                                q <= i__4;
                                ++q)
                        {
                            aaqq = sva[q];
                            if (aaqq > 0.f)
                            {
                                aapp0 = aapp;
                                /* .. M x 2 Jacobi SVD .. */
                                /* Safe Gram matrix computation */
                                if (aaqq >= 1.f)
                                {
                                    if (aapp >= aaqq)
                                    {
                                        rotok = small * aapp <= aaqq;
                                    }
                                    else
                                    {
                                        rotok = small * aaqq <= aapp;
                                    }
                                    if (aapp < big / aaqq)
                                    {
                                        aapq = sdot_(m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], & c__1) * work[p] * work[q] / aaqq / aapp;
                                    }
                                    else
                                    {
                                        scopy_(m, &a[p * a_dim1 + 1], &c__1, & work[*n + 1], &c__1);
                                        slascl_("G", &c__0, &c__0, &aapp, & work[p], m, &c__1, &work[*n + 1], lda, &ierr);
                                        aapq = sdot_(m, &work[*n + 1], &c__1, &a[q * a_dim1 + 1], &c__1) * work[q] / aaqq;
                                    }
                                }
                                else
                                {
                                    if (aapp >= aaqq)
                                    {
                                        rotok = aapp <= aaqq / small;
                                    }
                                    else
                                    {
                                        rotok = aaqq <= aapp / small;
                                    }
                                    if (aapp > small / aaqq)
                                    {
                                        aapq = sdot_(m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], & c__1) * work[p] * work[q] / aaqq / aapp;
                                    }
                                    else
                                    {
                                        scopy_(m, &a[q * a_dim1 + 1], &c__1, & work[*n + 1], &c__1);
                                        slascl_("G", &c__0, &c__0, &aaqq, & work[q], m, &c__1, &work[*n + 1], lda, &ierr);
                                        aapq = sdot_(m, &work[*n + 1], &c__1, &a[p * a_dim1 + 1], &c__1) * work[p] / aapp;
                                    }
                                }
                                /* Computing MAX */
                                r__1 = mxaapq;
                                r__2 = f2c_abs(aapq); // , expr subst
                                mxaapq = max(r__1,r__2);
                                /* TO rotate or NOT to rotate, THAT is the question ... */
                                if (f2c_abs(aapq) > tol)
                                {
                                    notrot = 0;
                                    /* [RTD] ROTATED = ROTATED + 1 */
                                    pskipped = 0;
                                    ++iswrot;
                                    if (rotok)
                                    {
                                        aqoap = aaqq / aapp;
                                        apoaq = aapp / aaqq;
                                        theta = (r__1 = aqoap - apoaq, f2c_abs( r__1)) * -.5f / aapq;
                                        if (aaqq > aapp0)
                                        {
                                            theta = -theta;
                                        }
                                        if (f2c_abs(theta) > bigtheta)
                                        {
                                            t = .5f / theta;
                                            fastr[2] = t * work[p] / work[q];
                                            fastr[3] = -t * work[q] / work[p];
                                            srotm_(m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], &c__1, fastr);
                                            if (rsvec)
                                            {
                                                srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &c__1, fastr);
                                            }
                                            /* Computing MAX */
                                            r__1 = 0.f;
                                            r__2 = t * apoaq * aapq + 1.f; // , expr subst
                                            sva[q] = aaqq * sqrt((max(r__1, r__2)));
                                            /* Computing MAX */
                                            r__1 = 0.f;
                                            r__2 = 1.f - t * aqoap * aapq; // , expr subst
                                            aapp *= sqrt((max(r__1,r__2)));
                                            /* Computing MAX */
                                            r__1 = mxsinj;
                                            r__2 = f2c_abs(t); // , expr subst
                                            mxsinj = max(r__1,r__2);
                                        }
                                        else
                                        {
                                            /* .. choose correct signum for THETA and rotate */
                                            thsign = -r_sign(&c_b18, &aapq);
                                            if (aaqq > aapp0)
                                            {
                                                thsign = -thsign;
                                            }
                                            t = 1.f / (theta + thsign * sqrt( theta * theta + 1.f));
                                            cs = sqrt(1.f / (t * t + 1.f));
                                            sn = t * cs;
                                            /* Computing MAX */
                                            r__1 = mxsinj;
                                            r__2 = f2c_abs(sn); // , expr subst
                                            mxsinj = max(r__1,r__2);
                                            /* Computing MAX */
                                            r__1 = 0.f;
                                            r__2 = t * apoaq * aapq + 1.f; // , expr subst
                                            sva[q] = aaqq * sqrt((max(r__1, r__2)));
                                            /* Computing MAX */
                                            r__1 = 0.f;
                                            r__2 = 1.f - t * aqoap * aapq; // , expr subst
                                            aapp *= sqrt((max(r__1,r__2)));
                                            apoaq = work[p] / work[q];
                                            aqoap = work[q] / work[p];
                                            if (work[p] >= 1.f)
                                            {
                                                if (work[q] >= 1.f)
                                                {
                                                    fastr[2] = t * apoaq;
                                                    fastr[3] = -t * aqoap;
                                                    work[p] *= cs;
                                                    work[q] *= cs;
                                                    srotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1, fastr);
                                                    if (rsvec)
                                                    {
                                                        srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[ q * v_dim1 + 1], &c__1, fastr);
                                                    }
                                                }
                                                else
                                                {
                                                    r__1 = -t * aqoap;
                                                    saxpy_(m, &r__1, &a[q * a_dim1 + 1], &c__1, &a[ p * a_dim1 + 1], &c__1);
                                                    r__1 = cs * sn * apoaq;
                                                    saxpy_(m, &r__1, &a[p * a_dim1 + 1], &c__1, &a[ q * a_dim1 + 1], &c__1);
                                                    if (rsvec)
                                                    {
                                                        r__1 = -t * aqoap;
                                                        saxpy_(&mvl, &r__1, &v[q * v_dim1 + 1], & c__1, &v[p * v_dim1 + 1], &c__1);
                                                        r__1 = cs * sn * apoaq;
                                                        saxpy_(&mvl, &r__1, &v[p * v_dim1 + 1], & c__1, &v[q * v_dim1 + 1], &c__1);
                                                    }
                                                    work[p] *= cs;
                                                    work[q] /= cs;
                                                }
                                            }
                                            else
                                            {
                                                if (work[q] >= 1.f)
                                                {
                                                    r__1 = t * apoaq;
                                                    saxpy_(m, &r__1, &a[p * a_dim1 + 1], &c__1, &a[ q * a_dim1 + 1], &c__1);
                                                    r__1 = -cs * sn * aqoap;
                                                    saxpy_(m, &r__1, &a[q * a_dim1 + 1], &c__1, &a[ p * a_dim1 + 1], &c__1);
                                                    if (rsvec)
                                                    {
                                                        r__1 = t * apoaq;
                                                        saxpy_(&mvl, &r__1, &v[p * v_dim1 + 1], & c__1, &v[q * v_dim1 + 1], &c__1);
                                                        r__1 = -cs * sn * aqoap;
                                                        saxpy_(&mvl, &r__1, &v[q * v_dim1 + 1], & c__1, &v[p * v_dim1 + 1], &c__1);
                                                    }
                                                    work[p] /= cs;
                                                    work[q] *= cs;
                                                }
                                                else
                                                {
                                                    if (work[p] >= work[q])
                                                    {
                                                        r__1 = -t * aqoap;
                                                        saxpy_(m, &r__1, &a[q * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &c__1);
                                                        r__1 = cs * sn * apoaq;
                                                        saxpy_(m, &r__1, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
                                                        work[p] *= cs;
                                                        work[q] /= cs;
                                                        if (rsvec)
                                                        {
                                                            r__1 = -t * aqoap;
                                                            saxpy_(&mvl, &r__1, &v[q * v_dim1 + 1], &c__1, &v[p * v_dim1 + 1], & c__1);
                                                            r__1 = cs * sn * apoaq;
                                                            saxpy_(&mvl, &r__1, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], & c__1);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        r__1 = t * apoaq;
                                                        saxpy_(m, &r__1, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
                                                        r__1 = -cs * sn * aqoap;
                                                        saxpy_(m, &r__1, &a[q * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &c__1);
                                                        work[p] /= cs;
                                                        work[q] *= cs;
                                                        if (rsvec)
                                                        {
                                                            r__1 = t * apoaq;
                                                            saxpy_(&mvl, &r__1, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], & c__1);
                                                            r__1 = -cs * sn * aqoap;
                                                            saxpy_(&mvl, &r__1, &v[q * v_dim1 + 1], &c__1, &v[p * v_dim1 + 1], & c__1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (aapp > aaqq)
                                        {
                                            scopy_(m, &a[p * a_dim1 + 1], & c__1, &work[*n + 1], & c__1);
                                            slascl_("G", &c__0, &c__0, &aapp, &c_b18, m, &c__1, &work[* n + 1], lda, &ierr);
                                            slascl_("G", &c__0, &c__0, &aaqq, &c_b18, m, &c__1, &a[q * a_dim1 + 1], lda, &ierr);
                                            temp1 = -aapq * work[p] / work[q];
                                            saxpy_(m, &temp1, &work[*n + 1], & c__1, &a[q * a_dim1 + 1], &c__1);
                                            slascl_("G", &c__0, &c__0, &c_b18, &aaqq, m, &c__1, &a[q * a_dim1 + 1], lda, &ierr);
                                            /* Computing MAX */
                                            r__1 = 0.f;
                                            r__2 = 1.f - aapq * aapq; // , expr subst
                                            sva[q] = aaqq * sqrt((max(r__1, r__2)));
                                            mxsinj = max(mxsinj,sfmin);
                                        }
                                        else
                                        {
                                            scopy_(m, &a[q * a_dim1 + 1], & c__1, &work[*n + 1], & c__1);
                                            slascl_("G", &c__0, &c__0, &aaqq, &c_b18, m, &c__1, &work[* n + 1], lda, &ierr);
                                            slascl_("G", &c__0, &c__0, &aapp, &c_b18, m, &c__1, &a[p * a_dim1 + 1], lda, &ierr);
                                            temp1 = -aapq * work[q] / work[p];
                                            saxpy_(m, &temp1, &work[*n + 1], & c__1, &a[p * a_dim1 + 1], &c__1);
                                            slascl_("G", &c__0, &c__0, &c_b18, &aapp, m, &c__1, &a[p * a_dim1 + 1], lda, &ierr);
                                            /* Computing MAX */
                                            r__1 = 0.f;
                                            r__2 = 1.f - aapq * aapq; // , expr subst
                                            sva[p] = aapp * sqrt((max(r__1, r__2)));
                                            mxsinj = max(mxsinj,sfmin);
                                        }
                                    }
                                    /* END IF ROTOK THEN ... ELSE */
                                    /* In the case of cancellation in updating SVA(q) */
                                    /* .. recompute SVA(q) */
                                    /* Computing 2nd power */
                                    r__1 = sva[q] / aaqq;
                                    if (r__1 * r__1 <= rooteps)
                                    {
                                        if (aaqq < rootbig && aaqq > rootsfmin)
                                        {
                                            sva[q] = snrm2_(m, &a[q * a_dim1 + 1], &c__1) * work[q];
                                        }
                                        else
                                        {
                                            t = 0.f;
                                            aaqq = 1.f;
                                            slassq_(m, &a[q * a_dim1 + 1], & c__1, &t, &aaqq);
                                            sva[q] = t * sqrt(aaqq) * work[q];
                                        }
                                    }
                                    /* Computing 2nd power */
                                    r__1 = aapp / aapp0;
                                    if (r__1 * r__1 <= rooteps)
                                    {
                                        if (aapp < rootbig && aapp > rootsfmin)
                                        {
                                            aapp = snrm2_(m, &a[p * a_dim1 + 1], &c__1) * work[p];
                                        }
                                        else
                                        {
                                            t = 0.f;
                                            aapp = 1.f;
                                            slassq_(m, &a[p * a_dim1 + 1], & c__1, &t, &aapp);
                                            aapp = t * sqrt(aapp) * work[p];
                                        }
                                        sva[p] = aapp;
                                    }
                                    /* end of OK rotation */
                                }
                                else
                                {
                                    ++notrot;
                                    /* [RTD] SKIPPED = SKIPPED + 1 */
                                    ++pskipped;
                                    ++ijblsk;
                                }
                            }
                            else
                            {
                                ++notrot;
                                ++pskipped;
                                ++ijblsk;
                            }
                            if (i__ <= swband && ijblsk >= blskip)
                            {
                                sva[p] = aapp;
                                notrot = 0;
                                goto L2011;
                            }
                            if (i__ <= swband && pskipped > rowskip)
                            {
                                aapp = -aapp;
                                notrot = 0;
                                goto L2203;
                            }
                            /* L2200: */
                        }
                        /* end of the q-loop */
L2203:
                        sva[p] = aapp;
                    }
                    else
                    {
                        if (aapp == 0.f)
                        {
                            /* Computing MIN */
                            i__4 = jgl + kbl - 1;
                            notrot = notrot + min(i__4,*n) - jgl + 1;
                        }
                        if (aapp < 0.f)
                        {
                            notrot = 0;
                        }
                    }
                    /* L2100: */
                }
                /* end of the p-loop */
                /* L2010: */
            }
            /* end of the jbc-loop */
L2011: /* 2011 bailed out of the jbc-loop */
            /* Computing MIN */
            i__3 = igl + kbl - 1;
            i__2 = min(i__3,*n);
            for (p = igl;
                    p <= i__2;
                    ++p)
            {
                sva[p] = (r__1 = sva[p], f2c_abs(r__1));
                /* L2012: */
            }
            /* ** */
            /* L2000: */
        }
        /* 2000 :: end of the ibr-loop */
        /* .. update SVA(N) */
        if (sva[*n] < rootbig && sva[*n] > rootsfmin)
        {
            sva[*n] = snrm2_(m, &a[*n * a_dim1 + 1], &c__1) * work[*n];
        }
        else
        {
            t = 0.f;
            aapp = 1.f;
            slassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
            sva[*n] = t * sqrt(aapp) * work[*n];
        }
        /* Additional steering devices */
        if (i__ < swband && (mxaapq <= roottol || iswrot <= *n))
        {
            swband = i__;
        }
        if (i__ > swband + 1 && mxaapq < sqrt((real) (*n)) * tol && (real) (* n) * mxaapq * mxsinj < tol)
        {
            goto L1994;
        }
        if (notrot >= emptsw)
        {
            goto L1994;
        }
        /* L1993: */
    }
    /* end i=1:NSWEEP loop */
    /* #:( Reaching this point means that the procedure has not converged. */
    *info = 29;
    goto L1995;
L1994: /* #:) Reaching this point means numerical convergence after the i-th */
    /* sweep. */
    *info = 0;
    /* #:) INFO = 0 confirms successful iterations. */
L1995: /* Sort the singular values and find how many are above */
    /* the underflow threshold. */
    n2 = 0;
    n4 = 0;
    i__1 = *n - 1;
    for (p = 1;
            p <= i__1;
            ++p)
    {
        i__2 = *n - p + 1;
        q = isamax_(&i__2, &sva[p], &c__1) + p - 1;
        if (p != q)
        {
            temp1 = sva[p];
            sva[p] = sva[q];
            sva[q] = temp1;
            temp1 = work[p];
            work[p] = work[q];
            work[q] = temp1;
            sswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
            if (rsvec)
            {
                sswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], & c__1);
            }
        }
        if (sva[p] != 0.f)
        {
            ++n4;
            if (sva[p] * skl > sfmin)
            {
                ++n2;
            }
        }
        /* L5991: */
    }
    if (sva[*n] != 0.f)
    {
        ++n4;
        if (sva[*n] * skl > sfmin)
        {
            ++n2;
        }
    }
    /* Normalize the left singular vectors. */
    if (lsvec || uctol)
    {
        i__1 = n2;
        for (p = 1;
                p <= i__1;
                ++p)
        {
            r__1 = work[p] / sva[p];
            sscal_(m, &r__1, &a[p * a_dim1 + 1], &c__1);
            /* L1998: */
        }
    }
    /* Scale the product of Jacobi rotations (assemble the fast rotations). */
    if (rsvec)
    {
        if (applv)
        {
            i__1 = *n;
            for (p = 1;
                    p <= i__1;
                    ++p)
            {
                sscal_(&mvl, &work[p], &v[p * v_dim1 + 1], &c__1);
                /* L2398: */
            }
        }
        else
        {
            i__1 = *n;
            for (p = 1;
                    p <= i__1;
                    ++p)
            {
                temp1 = 1.f / snrm2_(&mvl, &v[p * v_dim1 + 1], &c__1);
                sscal_(&mvl, &temp1, &v[p * v_dim1 + 1], &c__1);
                /* L2399: */
            }
        }
    }
    /* Undo scaling, if necessary (and possible). */
    if (skl > 1.f && sva[1] < big / skl || skl < 1.f && sva[max(n2,1)] > sfmin / skl)
    {
        i__1 = *n;
        for (p = 1;
                p <= i__1;
                ++p)
        {
            sva[p] = skl * sva[p];
            /* L2400: */
        }
        skl = 1.f;
    }
    work[1] = skl;
    /* The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE */
    /* then some of the singular values may overflow or underflow and */
    /* the spectrum is given in this factored representation. */
    work[2] = (real) n4;
    /* N4 is the number of computed nonzero singular values of A. */
    work[3] = (real) n2;
    /* N2 is the number of singular values of A greater than SFMIN. */
    /* If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers */
    /* that may carry some information. */
    work[4] = (real) i__;
    /* i is the index of the last sweep before declaring convergence. */
    work[5] = mxaapq;
    /* MXAAPQ is the largest absolute value of scaled pivots in the */
    /* last sweep */
    work[6] = mxsinj;
    /* MXSINJ is the largest absolute value of the sines of Jacobi angles */
    /* in the last sweep */
    return 0;
    /* .. */
    /* .. END OF SGESVJ */
    /* .. */
}
/* sgesvj_ */
