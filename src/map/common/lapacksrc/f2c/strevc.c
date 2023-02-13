/* ../netlib/strevc.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static logical c_false = FALSE_;
static integer c__1 = 1;
static real c_b22 = 1.f;
static real c_b25 = 0.f;
static integer c__2 = 2;
static logical c_true = TRUE_;
/* > \brief \b STREVC */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STREVC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strevc. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strevc. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strevc. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, */
/* LDVR, MM, M, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER HOWMNY, SIDE */
/* INTEGER INFO, LDT, LDVL, LDVR, M, MM, N */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL SELECT( * ) */
/* REAL T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STREVC computes some or all of the right and/or left eigenvectors of */
/* > a real upper quasi-triangular matrix T. */
/* > Matrices of this type are produced by the Schur factorization of */
/* > a real general matrix: A = Q*T*Q**T, as computed by SHSEQR. */
/* > */
/* > The right eigenvector x and the left eigenvector y of T corresponding */
/* > to an eigenvalue w are defined by: */
/* > */
/* > T*x = w*x, (y**T)*T = w*(y**T) */
/* > */
/* > where y**T denotes the transpose of y. */
/* > The eigenvalues are not input to this routine, but are read directly */
/* > from the diagonal blocks of T. */
/* > */
/* > This routine returns the matrices X and/or Y of right and left */
/* > eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an */
/* > input matrix. If Q is the orthogonal factor that reduces a matrix */
/* > A to Schur form T, then Q*X and Q*Y are the matrices of right and */
/* > left eigenvectors of A. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'R': compute right eigenvectors only;
*/
/* > = 'L': compute left eigenvectors only;
*/
/* > = 'B': compute both right and left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] HOWMNY */
/* > \verbatim */
/* > HOWMNY is CHARACTER*1 */
/* > = 'A': compute all right and/or left eigenvectors;
*/
/* > = 'B': compute all right and/or left eigenvectors, */
/* > backtransformed by the matrices in VR and/or VL;
*/
/* > = 'S': compute selected right and/or left eigenvectors, */
/* > as indicated by the logical array SELECT. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SELECT */
/* > \verbatim */
/* > SELECT is LOGICAL array, dimension (N) */
/* > If HOWMNY = 'S', SELECT specifies the eigenvectors to be */
/* > computed. */
/* > If w(j) is a real eigenvalue, the corresponding real */
/* > eigenvector is computed if SELECT(j) is .TRUE.. */
/* > If w(j) and w(j+1) are the real and imaginary parts of a */
/* > complex eigenvalue, the corresponding complex eigenvector is */
/* > computed if either SELECT(j) or SELECT(j+1) is .TRUE., and */
/* > on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to */
/* > .FALSE.. */
/* > Not referenced if HOWMNY = 'A' or 'B'. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix T. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* > T is REAL array, dimension (LDT,N) */
/* > The upper quasi-triangular matrix T in Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* > VL is REAL array, dimension (LDVL,MM) */
/* > On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/* > contain an N-by-N matrix Q (usually the orthogonal matrix Q */
/* > of Schur vectors returned by SHSEQR). */
/* > On exit, if SIDE = 'L' or 'B', VL contains: */
/* > if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
*/
/* > if HOWMNY = 'B', the matrix Q*Y;
*/
/* > if HOWMNY = 'S', the left eigenvectors of T specified by */
/* > SELECT, stored consecutively in the columns */
/* > of VL, in the same order as their */
/* > eigenvalues. */
/* > A complex eigenvector corresponding to a complex eigenvalue */
/* > is stored in two consecutive columns, the first holding the */
/* > real part, and the second the imaginary part. */
/* > Not referenced if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* > LDVL is INTEGER */
/* > The leading dimension of the array VL. LDVL >= 1, and if */
/* > SIDE = 'L' or 'B', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* > VR is REAL array, dimension (LDVR,MM) */
/* > On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/* > contain an N-by-N matrix Q (usually the orthogonal matrix Q */
/* > of Schur vectors returned by SHSEQR). */
/* > On exit, if SIDE = 'R' or 'B', VR contains: */
/* > if HOWMNY = 'A', the matrix X of right eigenvectors of T;
*/
/* > if HOWMNY = 'B', the matrix Q*X;
*/
/* > if HOWMNY = 'S', the right eigenvectors of T specified by */
/* > SELECT, stored consecutively in the columns */
/* > of VR, in the same order as their */
/* > eigenvalues. */
/* > A complex eigenvector corresponding to a complex eigenvalue */
/* > is stored in two consecutive columns, the first holding the */
/* > real part and the second the imaginary part. */
/* > Not referenced if SIDE = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* > LDVR is INTEGER */
/* > The leading dimension of the array VR. LDVR >= 1, and if */
/* > SIDE = 'R' or 'B', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] MM */
/* > \verbatim */
/* > MM is INTEGER */
/* > The number of columns in the arrays VL and/or VR. MM >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of columns in the arrays VL and/or VR actually */
/* > used to store the eigenvectors. */
/* > If HOWMNY = 'A' or 'B', M is set to N. */
/* > Each selected real eigenvector occupies one column and each */
/* > selected complex eigenvector occupies two columns. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The algorithm used in this program is basically backward (forward) */
/* > substitution, with scaling to make the the code robust against */
/* > possible overflow. */
/* > */
/* > Each eigenvector is normalized so that the element of largest */
/* > magnitude has magnitude 1;
here the magnitude of a complex number */
/* > (x,y) is taken to be |x| + |y|. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int strevc_(char *side, char *howmny, logical *select, integer *n, real *t, integer *ldt, real *vl, integer *ldvl, real *vr, integer *ldvr, integer *mm, integer *m, real *work, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k;
    real x[4] /* was [2][2] */
    ;
    integer j1, j2, n2, ii, ki, ip, is;
    real wi, wr, rec, ulp, beta, emax;
    logical pair, allv;
    integer ierr;
    real unfl, ovfl, smin;
    extern real sdot_(integer *, real *, integer *, real *, integer *);
    logical over;
    real vmax;
    integer jnxt;
    real scale;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    real remax;
    logical leftv;
    extern /* Subroutine */
    int sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    logical bothv;
    real vcrit;
    logical somev;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *);
    real xnorm;
    extern /* Subroutine */
    int saxpy_(integer *, real *, real *, integer *, real *, integer *), slaln2_(logical *, integer *, integer *, real *, real *, real *, integer *, real *, real *, real *, integer *, real *, real *, real *, integer *, real *, real *, integer *), slabad_(real *, real *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    real bignum;
    extern integer isamax_(integer *, real *, integer *);
    logical rightv;
    real smlnum;
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Decode and test the input parameters */
    /* Parameter adjustments */
    --select;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    /* Function Body */
    bothv = lsame_(side, "B");
    rightv = lsame_(side, "R") || bothv;
    leftv = lsame_(side, "L") || bothv;
    allv = lsame_(howmny, "A");
    over = lsame_(howmny, "B");
    somev = lsame_(howmny, "S");
    *info = 0;
    if (! rightv && ! leftv)
    {
        *info = -1;
    }
    else if (! allv && ! over && ! somev)
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*ldt < max(1,*n))
    {
        *info = -6;
    }
    else if (*ldvl < 1 || leftv && *ldvl < *n)
    {
        *info = -8;
    }
    else if (*ldvr < 1 || rightv && *ldvr < *n)
    {
        *info = -10;
    }
    else
    {
        /* Set M to the number of columns required to store the selected */
        /* eigenvectors, standardize the array SELECT if necessary, and */
        /* test MM. */
        if (somev)
        {
            *m = 0;
            pair = FALSE_;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (pair)
                {
                    pair = FALSE_;
                    select[j] = FALSE_;
                }
                else
                {
                    if (j < *n)
                    {
                        if (t[j + 1 + j * t_dim1] == 0.f)
                        {
                            if (select[j])
                            {
                                ++(*m);
                            }
                        }
                        else
                        {
                            pair = TRUE_;
                            if (select[j] || select[j + 1])
                            {
                                select[j] = TRUE_;
                                *m += 2;
                            }
                        }
                    }
                    else
                    {
                        if (select[*n])
                        {
                            ++(*m);
                        }
                    }
                }
                /* L10: */
            }
        }
        else
        {
            *m = *n;
        }
        if (*mm < *m)
        {
            *info = -11;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("STREVC", &i__1);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0)
    {
        return 0;
    }
    /* Set the constants to control overflow. */
    unfl = slamch_("Safe minimum");
    ovfl = 1.f / unfl;
    slabad_(&unfl, &ovfl);
    ulp = slamch_("Precision");
    smlnum = unfl * (*n / ulp);
    bignum = (1.f - ulp) / smlnum;
    /* Compute 1-norm of each column of strictly upper triangular */
    /* part of T to control overflow in triangular solver. */
    work[1] = 0.f;
    i__1 = *n;
    for (j = 2;
            j <= i__1;
            ++j)
    {
        work[j] = 0.f;
        i__2 = j - 1;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            work[j] += (r__1 = t[i__ + j * t_dim1], f2c_abs(r__1));
            /* L20: */
        }
        /* L30: */
    }
    /* Index IP is used to specify the real or complex eigenvalue: */
    /* IP = 0, real eigenvalue, */
    /* 1, first of conjugate complex pair: (wr,wi) */
    /* -1, second of conjugate complex pair: (wr,wi) */
    n2 = *n << 1;
    if (rightv)
    {
        /* Compute right eigenvectors. */
        ip = 0;
        is = *m;
        for (ki = *n;
                ki >= 1;
                --ki)
        {
            if (ip == 1)
            {
                goto L130;
            }
            if (ki == 1)
            {
                goto L40;
            }
            if (t[ki + (ki - 1) * t_dim1] == 0.f)
            {
                goto L40;
            }
            ip = -1;
L40:
            if (somev)
            {
                if (ip == 0)
                {
                    if (! select[ki])
                    {
                        goto L130;
                    }
                }
                else
                {
                    if (! select[ki - 1])
                    {
                        goto L130;
                    }
                }
            }
            /* Compute the KI-th eigenvalue (WR,WI). */
            wr = t[ki + ki * t_dim1];
            wi = 0.f;
            if (ip != 0)
            {
                wi = sqrt((r__1 = t[ki + (ki - 1) * t_dim1], f2c_abs(r__1))) * sqrt((r__2 = t[ki - 1 + ki * t_dim1], f2c_abs(r__2)));
            }
            /* Computing MAX */
            r__1 = ulp * (f2c_abs(wr) + f2c_abs(wi));
            smin = max(r__1,smlnum);
            if (ip == 0)
            {
                /* Real right eigenvector */
                work[ki + *n] = 1.f;
                /* Form right-hand side */
                i__1 = ki - 1;
                for (k = 1;
                        k <= i__1;
                        ++k)
                {
                    work[k + *n] = -t[k + ki * t_dim1];
                    /* L50: */
                }
                /* Solve the upper quasi-triangular system: */
                /* (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK. */
                jnxt = ki - 1;
                for (j = ki - 1;
                        j >= 1;
                        --j)
                {
                    if (j > jnxt)
                    {
                        goto L60;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j - 1;
                    if (j > 1)
                    {
                        if (t[j + (j - 1) * t_dim1] != 0.f)
                        {
                            j1 = j - 1;
                            jnxt = j - 2;
                        }
                    }
                    if (j1 == j2)
                    {
                        /* 1-by-1 diagonal block */
                        slaln2_(&c_false, &c__1, &c__1, &smin, &c_b22, &t[j + j * t_dim1], ldt, &c_b22, &c_b22, &work[j + * n], n, &wr, &c_b25, x, &c__2, &scale, &xnorm, &ierr);
                        /* Scale X(1,1) to avoid overflow when updating */
                        /* the right-hand side. */
                        if (xnorm > 1.f)
                        {
                            if (work[j] > bignum / xnorm)
                            {
                                x[0] /= xnorm;
                                scale /= xnorm;
                            }
                        }
                        /* Scale if necessary */
                        if (scale != 1.f)
                        {
                            sscal_(&ki, &scale, &work[*n + 1], &c__1);
                        }
                        work[j + *n] = x[0];
                        /* Update right-hand side */
                        i__1 = j - 1;
                        r__1 = -x[0];
                        saxpy_(&i__1, &r__1, &t[j * t_dim1 + 1], &c__1, &work[ *n + 1], &c__1);
                    }
                    else
                    {
                        /* 2-by-2 diagonal block */
                        slaln2_(&c_false, &c__2, &c__1, &smin, &c_b22, &t[j - 1 + (j - 1) * t_dim1], ldt, &c_b22, &c_b22, & work[j - 1 + *n], n, &wr, &c_b25, x, &c__2, & scale, &xnorm, &ierr);
                        /* Scale X(1,1) and X(2,1) to avoid overflow when */
                        /* updating the right-hand side. */
                        if (xnorm > 1.f)
                        {
                            /* Computing MAX */
                            r__1 = work[j - 1];
                            r__2 = work[j]; // , expr subst
                            beta = max(r__1,r__2);
                            if (beta > bignum / xnorm)
                            {
                                x[0] /= xnorm;
                                x[1] /= xnorm;
                                scale /= xnorm;
                            }
                        }
                        /* Scale if necessary */
                        if (scale != 1.f)
                        {
                            sscal_(&ki, &scale, &work[*n + 1], &c__1);
                        }
                        work[j - 1 + *n] = x[0];
                        work[j + *n] = x[1];
                        /* Update right-hand side */
                        i__1 = j - 2;
                        r__1 = -x[0];
                        saxpy_(&i__1, &r__1, &t[(j - 1) * t_dim1 + 1], &c__1, &work[*n + 1], &c__1);
                        i__1 = j - 2;
                        r__1 = -x[1];
                        saxpy_(&i__1, &r__1, &t[j * t_dim1 + 1], &c__1, &work[ *n + 1], &c__1);
                    }
L60:
                    ;
                }
                /* Copy the vector x or Q*x to VR and normalize. */
                if (! over)
                {
                    scopy_(&ki, &work[*n + 1], &c__1, &vr[is * vr_dim1 + 1], & c__1);
                    ii = isamax_(&ki, &vr[is * vr_dim1 + 1], &c__1);
                    remax = 1.f / (r__1 = vr[ii + is * vr_dim1], f2c_abs(r__1));
                    sscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);
                    i__1 = *n;
                    for (k = ki + 1;
                            k <= i__1;
                            ++k)
                    {
                        vr[k + is * vr_dim1] = 0.f;
                        /* L70: */
                    }
                }
                else
                {
                    if (ki > 1)
                    {
                        i__1 = ki - 1;
                        sgemv_("N", n, &i__1, &c_b22, &vr[vr_offset], ldvr, & work[*n + 1], &c__1, &work[ki + *n], &vr[ki * vr_dim1 + 1], &c__1);
                    }
                    ii = isamax_(n, &vr[ki * vr_dim1 + 1], &c__1);
                    remax = 1.f / (r__1 = vr[ii + ki * vr_dim1], f2c_abs(r__1));
                    sscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
                }
            }
            else
            {
                /* Complex right eigenvector. */
                /* Initial solve */
                /* [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0. */
                /* [ (T(KI,KI-1) T(KI,KI) ) ] */
                if ((r__1 = t[ki - 1 + ki * t_dim1], f2c_abs(r__1)) >= (r__2 = t[ ki + (ki - 1) * t_dim1], f2c_abs(r__2)))
                {
                    work[ki - 1 + *n] = 1.f;
                    work[ki + n2] = wi / t[ki - 1 + ki * t_dim1];
                }
                else
                {
                    work[ki - 1 + *n] = -wi / t[ki + (ki - 1) * t_dim1];
                    work[ki + n2] = 1.f;
                }
                work[ki + *n] = 0.f;
                work[ki - 1 + n2] = 0.f;
                /* Form right-hand side */
                i__1 = ki - 2;
                for (k = 1;
                        k <= i__1;
                        ++k)
                {
                    work[k + *n] = -work[ki - 1 + *n] * t[k + (ki - 1) * t_dim1];
                    work[k + n2] = -work[ki + n2] * t[k + ki * t_dim1];
                    /* L80: */
                }
                /* Solve upper quasi-triangular system: */
                /* (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2) */
                jnxt = ki - 2;
                for (j = ki - 2;
                        j >= 1;
                        --j)
                {
                    if (j > jnxt)
                    {
                        goto L90;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j - 1;
                    if (j > 1)
                    {
                        if (t[j + (j - 1) * t_dim1] != 0.f)
                        {
                            j1 = j - 1;
                            jnxt = j - 2;
                        }
                    }
                    if (j1 == j2)
                    {
                        /* 1-by-1 diagonal block */
                        slaln2_(&c_false, &c__1, &c__2, &smin, &c_b22, &t[j + j * t_dim1], ldt, &c_b22, &c_b22, &work[j + * n], n, &wr, &wi, x, &c__2, &scale, &xnorm, & ierr);
                        /* Scale X(1,1) and X(1,2) to avoid overflow when */
                        /* updating the right-hand side. */
                        if (xnorm > 1.f)
                        {
                            if (work[j] > bignum / xnorm)
                            {
                                x[0] /= xnorm;
                                x[2] /= xnorm;
                                scale /= xnorm;
                            }
                        }
                        /* Scale if necessary */
                        if (scale != 1.f)
                        {
                            sscal_(&ki, &scale, &work[*n + 1], &c__1);
                            sscal_(&ki, &scale, &work[n2 + 1], &c__1);
                        }
                        work[j + *n] = x[0];
                        work[j + n2] = x[2];
                        /* Update the right-hand side */
                        i__1 = j - 1;
                        r__1 = -x[0];
                        saxpy_(&i__1, &r__1, &t[j * t_dim1 + 1], &c__1, &work[ *n + 1], &c__1);
                        i__1 = j - 1;
                        r__1 = -x[2];
                        saxpy_(&i__1, &r__1, &t[j * t_dim1 + 1], &c__1, &work[ n2 + 1], &c__1);
                    }
                    else
                    {
                        /* 2-by-2 diagonal block */
                        slaln2_(&c_false, &c__2, &c__2, &smin, &c_b22, &t[j - 1 + (j - 1) * t_dim1], ldt, &c_b22, &c_b22, & work[j - 1 + *n], n, &wr, &wi, x, &c__2, & scale, &xnorm, &ierr);
                        /* Scale X to avoid overflow when updating */
                        /* the right-hand side. */
                        if (xnorm > 1.f)
                        {
                            /* Computing MAX */
                            r__1 = work[j - 1];
                            r__2 = work[j]; // , expr subst
                            beta = max(r__1,r__2);
                            if (beta > bignum / xnorm)
                            {
                                rec = 1.f / xnorm;
                                x[0] *= rec;
                                x[2] *= rec;
                                x[1] *= rec;
                                x[3] *= rec;
                                scale *= rec;
                            }
                        }
                        /* Scale if necessary */
                        if (scale != 1.f)
                        {
                            sscal_(&ki, &scale, &work[*n + 1], &c__1);
                            sscal_(&ki, &scale, &work[n2 + 1], &c__1);
                        }
                        work[j - 1 + *n] = x[0];
                        work[j + *n] = x[1];
                        work[j - 1 + n2] = x[2];
                        work[j + n2] = x[3];
                        /* Update the right-hand side */
                        i__1 = j - 2;
                        r__1 = -x[0];
                        saxpy_(&i__1, &r__1, &t[(j - 1) * t_dim1 + 1], &c__1, &work[*n + 1], &c__1);
                        i__1 = j - 2;
                        r__1 = -x[1];
                        saxpy_(&i__1, &r__1, &t[j * t_dim1 + 1], &c__1, &work[ *n + 1], &c__1);
                        i__1 = j - 2;
                        r__1 = -x[2];
                        saxpy_(&i__1, &r__1, &t[(j - 1) * t_dim1 + 1], &c__1, &work[n2 + 1], &c__1);
                        i__1 = j - 2;
                        r__1 = -x[3];
                        saxpy_(&i__1, &r__1, &t[j * t_dim1 + 1], &c__1, &work[ n2 + 1], &c__1);
                    }
L90:
                    ;
                }
                /* Copy the vector x or Q*x to VR and normalize. */
                if (! over)
                {
                    scopy_(&ki, &work[*n + 1], &c__1, &vr[(is - 1) * vr_dim1 + 1], &c__1);
                    scopy_(&ki, &work[n2 + 1], &c__1, &vr[is * vr_dim1 + 1], & c__1);
                    emax = 0.f;
                    i__1 = ki;
                    for (k = 1;
                            k <= i__1;
                            ++k)
                    {
                        /* Computing MAX */
                        r__3 = emax;
                        r__4 = (r__1 = vr[k + (is - 1) * vr_dim1] , f2c_abs(r__1)) + (r__2 = vr[k + is * vr_dim1], f2c_abs(r__2)); // , expr subst
                        emax = max(r__3,r__4);
                        /* L100: */
                    }
                    remax = 1.f / emax;
                    sscal_(&ki, &remax, &vr[(is - 1) * vr_dim1 + 1], &c__1);
                    sscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);
                    i__1 = *n;
                    for (k = ki + 1;
                            k <= i__1;
                            ++k)
                    {
                        vr[k + (is - 1) * vr_dim1] = 0.f;
                        vr[k + is * vr_dim1] = 0.f;
                        /* L110: */
                    }
                }
                else
                {
                    if (ki > 2)
                    {
                        i__1 = ki - 2;
                        sgemv_("N", n, &i__1, &c_b22, &vr[vr_offset], ldvr, & work[*n + 1], &c__1, &work[ki - 1 + *n], &vr[( ki - 1) * vr_dim1 + 1], &c__1);
                        i__1 = ki - 2;
                        sgemv_("N", n, &i__1, &c_b22, &vr[vr_offset], ldvr, & work[n2 + 1], &c__1, &work[ki + n2], &vr[ki * vr_dim1 + 1], &c__1);
                    }
                    else
                    {
                        sscal_(n, &work[ki - 1 + *n], &vr[(ki - 1) * vr_dim1 + 1], &c__1);
                        sscal_(n, &work[ki + n2], &vr[ki * vr_dim1 + 1], & c__1);
                    }
                    emax = 0.f;
                    i__1 = *n;
                    for (k = 1;
                            k <= i__1;
                            ++k)
                    {
                        /* Computing MAX */
                        r__3 = emax;
                        r__4 = (r__1 = vr[k + (ki - 1) * vr_dim1] , f2c_abs(r__1)) + (r__2 = vr[k + ki * vr_dim1], f2c_abs(r__2)); // , expr subst
                        emax = max(r__3,r__4);
                        /* L120: */
                    }
                    remax = 1.f / emax;
                    sscal_(n, &remax, &vr[(ki - 1) * vr_dim1 + 1], &c__1);
                    sscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
                }
            }
            --is;
            if (ip != 0)
            {
                --is;
            }
L130:
            if (ip == 1)
            {
                ip = 0;
            }
            if (ip == -1)
            {
                ip = 1;
            }
            /* L140: */
        }
    }
    if (leftv)
    {
        /* Compute left eigenvectors. */
        ip = 0;
        is = 1;
        i__1 = *n;
        for (ki = 1;
                ki <= i__1;
                ++ki)
        {
            if (ip == -1)
            {
                goto L250;
            }
            if (ki == *n)
            {
                goto L150;
            }
            if (t[ki + 1 + ki * t_dim1] == 0.f)
            {
                goto L150;
            }
            ip = 1;
L150:
            if (somev)
            {
                if (! select[ki])
                {
                    goto L250;
                }
            }
            /* Compute the KI-th eigenvalue (WR,WI). */
            wr = t[ki + ki * t_dim1];
            wi = 0.f;
            if (ip != 0)
            {
                wi = sqrt((r__1 = t[ki + (ki + 1) * t_dim1], f2c_abs(r__1))) * sqrt((r__2 = t[ki + 1 + ki * t_dim1], f2c_abs(r__2)));
            }
            /* Computing MAX */
            r__1 = ulp * (f2c_abs(wr) + f2c_abs(wi));
            smin = max(r__1,smlnum);
            if (ip == 0)
            {
                /* Real left eigenvector. */
                work[ki + *n] = 1.f;
                /* Form right-hand side */
                i__2 = *n;
                for (k = ki + 1;
                        k <= i__2;
                        ++k)
                {
                    work[k + *n] = -t[ki + k * t_dim1];
                    /* L160: */
                }
                /* Solve the quasi-triangular system: */
                /* (T(KI+1:N,KI+1:N) - WR)**T*X = SCALE*WORK */
                vmax = 1.f;
                vcrit = bignum;
                jnxt = ki + 1;
                i__2 = *n;
                for (j = ki + 1;
                        j <= i__2;
                        ++j)
                {
                    if (j < jnxt)
                    {
                        goto L170;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j + 1;
                    if (j < *n)
                    {
                        if (t[j + 1 + j * t_dim1] != 0.f)
                        {
                            j2 = j + 1;
                            jnxt = j + 2;
                        }
                    }
                    if (j1 == j2)
                    {
                        /* 1-by-1 diagonal block */
                        /* Scale if necessary to avoid overflow when forming */
                        /* the right-hand side. */
                        if (work[j] > vcrit)
                        {
                            rec = 1.f / vmax;
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &rec, &work[ki + *n], &c__1);
                            vmax = 1.f;
                            vcrit = bignum;
                        }
                        i__3 = j - ki - 1;
                        work[j + *n] -= sdot_(&i__3, &t[ki + 1 + j * t_dim1], &c__1, &work[ki + 1 + *n], &c__1);
                        /* Solve (T(J,J)-WR)**T*X = WORK */
                        slaln2_(&c_false, &c__1, &c__1, &smin, &c_b22, &t[j + j * t_dim1], ldt, &c_b22, &c_b22, &work[j + * n], n, &wr, &c_b25, x, &c__2, &scale, &xnorm, &ierr);
                        /* Scale if necessary */
                        if (scale != 1.f)
                        {
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &scale, &work[ki + *n], &c__1);
                        }
                        work[j + *n] = x[0];
                        /* Computing MAX */
                        r__2 = (r__1 = work[j + *n], f2c_abs(r__1));
                        vmax = max(r__2,vmax);
                        vcrit = bignum / vmax;
                    }
                    else
                    {
                        /* 2-by-2 diagonal block */
                        /* Scale if necessary to avoid overflow when forming */
                        /* the right-hand side. */
                        /* Computing MAX */
                        r__1 = work[j];
                        r__2 = work[j + 1]; // , expr subst
                        beta = max(r__1,r__2);
                        if (beta > vcrit)
                        {
                            rec = 1.f / vmax;
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &rec, &work[ki + *n], &c__1);
                            vmax = 1.f;
                            vcrit = bignum;
                        }
                        i__3 = j - ki - 1;
                        work[j + *n] -= sdot_(&i__3, &t[ki + 1 + j * t_dim1], &c__1, &work[ki + 1 + *n], &c__1);
                        i__3 = j - ki - 1;
                        work[j + 1 + *n] -= sdot_(&i__3, &t[ki + 1 + (j + 1) * t_dim1], &c__1, &work[ki + 1 + *n], &c__1);
                        /* Solve */
                        /* [T(J,J)-WR T(J,J+1) ]**T* X = SCALE*( WORK1 ) */
                        /* [T(J+1,J) T(J+1,J+1)-WR] ( WORK2 ) */
                        slaln2_(&c_true, &c__2, &c__1, &smin, &c_b22, &t[j + j * t_dim1], ldt, &c_b22, &c_b22, &work[j + * n], n, &wr, &c_b25, x, &c__2, &scale, &xnorm, &ierr);
                        /* Scale if necessary */
                        if (scale != 1.f)
                        {
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &scale, &work[ki + *n], &c__1);
                        }
                        work[j + *n] = x[0];
                        work[j + 1 + *n] = x[1];
                        /* Computing MAX */
                        r__3 = (r__1 = work[j + *n], f2c_abs(r__1));
                        r__4 = (r__2 = work[j + 1 + *n], f2c_abs(r__2));
                        r__3 = max( r__3,r__4); // ; expr subst
                        vmax = max(r__3,vmax);
                        vcrit = bignum / vmax;
                    }
L170:
                    ;
                }
                /* Copy the vector x or Q*x to VL and normalize. */
                if (! over)
                {
                    i__2 = *n - ki + 1;
                    scopy_(&i__2, &work[ki + *n], &c__1, &vl[ki + is * vl_dim1], &c__1);
                    i__2 = *n - ki + 1;
                    ii = isamax_(&i__2, &vl[ki + is * vl_dim1], &c__1) + ki - 1;
                    remax = 1.f / (r__1 = vl[ii + is * vl_dim1], f2c_abs(r__1));
                    i__2 = *n - ki + 1;
                    sscal_(&i__2, &remax, &vl[ki + is * vl_dim1], &c__1);
                    i__2 = ki - 1;
                    for (k = 1;
                            k <= i__2;
                            ++k)
                    {
                        vl[k + is * vl_dim1] = 0.f;
                        /* L180: */
                    }
                }
                else
                {
                    if (ki < *n)
                    {
                        i__2 = *n - ki;
                        sgemv_("N", n, &i__2, &c_b22, &vl[(ki + 1) * vl_dim1 + 1], ldvl, &work[ki + 1 + *n], &c__1, &work[ ki + *n], &vl[ki * vl_dim1 + 1], &c__1);
                    }
                    ii = isamax_(n, &vl[ki * vl_dim1 + 1], &c__1);
                    remax = 1.f / (r__1 = vl[ii + ki * vl_dim1], f2c_abs(r__1));
                    sscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
                }
            }
            else
            {
                /* Complex left eigenvector. */
                /* Initial solve: */
                /* ((T(KI,KI) T(KI,KI+1) )**T - (WR - I* WI))*X = 0. */
                /* ((T(KI+1,KI) T(KI+1,KI+1)) ) */
                if ((r__1 = t[ki + (ki + 1) * t_dim1], f2c_abs(r__1)) >= (r__2 = t[ki + 1 + ki * t_dim1], f2c_abs(r__2)))
                {
                    work[ki + *n] = wi / t[ki + (ki + 1) * t_dim1];
                    work[ki + 1 + n2] = 1.f;
                }
                else
                {
                    work[ki + *n] = 1.f;
                    work[ki + 1 + n2] = -wi / t[ki + 1 + ki * t_dim1];
                }
                work[ki + 1 + *n] = 0.f;
                work[ki + n2] = 0.f;
                /* Form right-hand side */
                i__2 = *n;
                for (k = ki + 2;
                        k <= i__2;
                        ++k)
                {
                    work[k + *n] = -work[ki + *n] * t[ki + k * t_dim1];
                    work[k + n2] = -work[ki + 1 + n2] * t[ki + 1 + k * t_dim1] ;
                    /* L190: */
                }
                /* Solve complex quasi-triangular system: */
                /* ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2 */
                vmax = 1.f;
                vcrit = bignum;
                jnxt = ki + 2;
                i__2 = *n;
                for (j = ki + 2;
                        j <= i__2;
                        ++j)
                {
                    if (j < jnxt)
                    {
                        goto L200;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j + 1;
                    if (j < *n)
                    {
                        if (t[j + 1 + j * t_dim1] != 0.f)
                        {
                            j2 = j + 1;
                            jnxt = j + 2;
                        }
                    }
                    if (j1 == j2)
                    {
                        /* 1-by-1 diagonal block */
                        /* Scale if necessary to avoid overflow when */
                        /* forming the right-hand side elements. */
                        if (work[j] > vcrit)
                        {
                            rec = 1.f / vmax;
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &rec, &work[ki + *n], &c__1);
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &rec, &work[ki + n2], &c__1);
                            vmax = 1.f;
                            vcrit = bignum;
                        }
                        i__3 = j - ki - 2;
                        work[j + *n] -= sdot_(&i__3, &t[ki + 2 + j * t_dim1], &c__1, &work[ki + 2 + *n], &c__1);
                        i__3 = j - ki - 2;
                        work[j + n2] -= sdot_(&i__3, &t[ki + 2 + j * t_dim1], &c__1, &work[ki + 2 + n2], &c__1);
                        /* Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2 */
                        r__1 = -wi;
                        slaln2_(&c_false, &c__1, &c__2, &smin, &c_b22, &t[j + j * t_dim1], ldt, &c_b22, &c_b22, &work[j + * n], n, &wr, &r__1, x, &c__2, &scale, &xnorm, & ierr);
                        /* Scale if necessary */
                        if (scale != 1.f)
                        {
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &scale, &work[ki + *n], &c__1);
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &scale, &work[ki + n2], &c__1);
                        }
                        work[j + *n] = x[0];
                        work[j + n2] = x[2];
                        /* Computing MAX */
                        r__3 = (r__1 = work[j + *n], f2c_abs(r__1));
                        r__4 = (r__2 = work[j + n2], f2c_abs(r__2));
                        r__3 = max(r__3, r__4); // ; expr subst
                        vmax = max(r__3,vmax);
                        vcrit = bignum / vmax;
                    }
                    else
                    {
                        /* 2-by-2 diagonal block */
                        /* Scale if necessary to avoid overflow when forming */
                        /* the right-hand side elements. */
                        /* Computing MAX */
                        r__1 = work[j];
                        r__2 = work[j + 1]; // , expr subst
                        beta = max(r__1,r__2);
                        if (beta > vcrit)
                        {
                            rec = 1.f / vmax;
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &rec, &work[ki + *n], &c__1);
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &rec, &work[ki + n2], &c__1);
                            vmax = 1.f;
                            vcrit = bignum;
                        }
                        i__3 = j - ki - 2;
                        work[j + *n] -= sdot_(&i__3, &t[ki + 2 + j * t_dim1], &c__1, &work[ki + 2 + *n], &c__1);
                        i__3 = j - ki - 2;
                        work[j + n2] -= sdot_(&i__3, &t[ki + 2 + j * t_dim1], &c__1, &work[ki + 2 + n2], &c__1);
                        i__3 = j - ki - 2;
                        work[j + 1 + *n] -= sdot_(&i__3, &t[ki + 2 + (j + 1) * t_dim1], &c__1, &work[ki + 2 + *n], &c__1);
                        i__3 = j - ki - 2;
                        work[j + 1 + n2] -= sdot_(&i__3, &t[ki + 2 + (j + 1) * t_dim1], &c__1, &work[ki + 2 + n2], &c__1);
                        /* Solve 2-by-2 complex linear equation */
                        /* ([T(j,j) T(j,j+1) ]**T-(wr-i*wi)*I)*X = SCALE*B */
                        /* ([T(j+1,j) T(j+1,j+1)] ) */
                        r__1 = -wi;
                        slaln2_(&c_true, &c__2, &c__2, &smin, &c_b22, &t[j + j * t_dim1], ldt, &c_b22, &c_b22, &work[j + * n], n, &wr, &r__1, x, &c__2, &scale, &xnorm, & ierr);
                        /* Scale if necessary */
                        if (scale != 1.f)
                        {
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &scale, &work[ki + *n], &c__1);
                            i__3 = *n - ki + 1;
                            sscal_(&i__3, &scale, &work[ki + n2], &c__1);
                        }
                        work[j + *n] = x[0];
                        work[j + n2] = x[2];
                        work[j + 1 + *n] = x[1];
                        work[j + 1 + n2] = x[3];
                        /* Computing MAX */
                        r__1 = f2c_abs(x[0]), r__2 = f2c_abs(x[2]), r__1 = max(r__1, r__2), r__2 = f2c_abs(x[1]), r__1 = max(r__1,r__2) ;
                        r__2 = f2c_abs(x[3]);
                        r__1 = max(r__1,r__2); // ; expr subst
                        vmax = max(r__1,vmax);
                        vcrit = bignum / vmax;
                    }
L200:
                    ;
                }
                /* Copy the vector x or Q*x to VL and normalize. */
                if (! over)
                {
                    i__2 = *n - ki + 1;
                    scopy_(&i__2, &work[ki + *n], &c__1, &vl[ki + is * vl_dim1], &c__1);
                    i__2 = *n - ki + 1;
                    scopy_(&i__2, &work[ki + n2], &c__1, &vl[ki + (is + 1) * vl_dim1], &c__1);
                    emax = 0.f;
                    i__2 = *n;
                    for (k = ki;
                            k <= i__2;
                            ++k)
                    {
                        /* Computing MAX */
                        r__3 = emax;
                        r__4 = (r__1 = vl[k + is * vl_dim1], f2c_abs( r__1)) + (r__2 = vl[k + (is + 1) * vl_dim1], f2c_abs(r__2)); // , expr subst
                        emax = max(r__3,r__4);
                        /* L220: */
                    }
                    remax = 1.f / emax;
                    i__2 = *n - ki + 1;
                    sscal_(&i__2, &remax, &vl[ki + is * vl_dim1], &c__1);
                    i__2 = *n - ki + 1;
                    sscal_(&i__2, &remax, &vl[ki + (is + 1) * vl_dim1], &c__1) ;
                    i__2 = ki - 1;
                    for (k = 1;
                            k <= i__2;
                            ++k)
                    {
                        vl[k + is * vl_dim1] = 0.f;
                        vl[k + (is + 1) * vl_dim1] = 0.f;
                        /* L230: */
                    }
                }
                else
                {
                    if (ki < *n - 1)
                    {
                        i__2 = *n - ki - 1;
                        sgemv_("N", n, &i__2, &c_b22, &vl[(ki + 2) * vl_dim1 + 1], ldvl, &work[ki + 2 + *n], &c__1, &work[ ki + *n], &vl[ki * vl_dim1 + 1], &c__1);
                        i__2 = *n - ki - 1;
                        sgemv_("N", n, &i__2, &c_b22, &vl[(ki + 2) * vl_dim1 + 1], ldvl, &work[ki + 2 + n2], &c__1, &work[ ki + 1 + n2], &vl[(ki + 1) * vl_dim1 + 1], & c__1);
                    }
                    else
                    {
                        sscal_(n, &work[ki + *n], &vl[ki * vl_dim1 + 1], & c__1);
                        sscal_(n, &work[ki + 1 + n2], &vl[(ki + 1) * vl_dim1 + 1], &c__1);
                    }
                    emax = 0.f;
                    i__2 = *n;
                    for (k = 1;
                            k <= i__2;
                            ++k)
                    {
                        /* Computing MAX */
                        r__3 = emax;
                        r__4 = (r__1 = vl[k + ki * vl_dim1], f2c_abs( r__1)) + (r__2 = vl[k + (ki + 1) * vl_dim1], f2c_abs(r__2)); // , expr subst
                        emax = max(r__3,r__4);
                        /* L240: */
                    }
                    remax = 1.f / emax;
                    sscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
                    sscal_(n, &remax, &vl[(ki + 1) * vl_dim1 + 1], &c__1);
                }
            }
            ++is;
            if (ip != 0)
            {
                ++is;
            }
L250:
            if (ip == -1)
            {
                ip = 0;
            }
            if (ip == 1)
            {
                ip = -1;
            }
            /* L260: */
        }
    }
    return 0;
    /* End of STREVC */
}
/* strevc_ */
