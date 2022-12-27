/* ../netlib/slatbs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b36 = .5f;
/* > \brief \b SLATBS solves a triangular banded system of equations. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLATBS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatbs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatbs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatbs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X, */
/* SCALE, CNORM, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORMIN, TRANS, UPLO */
/* INTEGER INFO, KD, LDAB, N */
/* REAL SCALE */
/* .. */
/* .. Array Arguments .. */
/* REAL AB( LDAB, * ), CNORM( * ), X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLATBS solves one of the triangular systems */
/* > */
/* > A *x = s*b or A**T*x = s*b */
/* > */
/* > with scaling to prevent overflow, where A is an upper or lower */
/* > triangular band matrix. Here A**T denotes the transpose of A, x and b */
/* > are n-element vectors, and s is a scaling factor, usually less than */
/* > or equal to 1, chosen so that the components of x will be less than */
/* > the overflow threshold. If the unscaled problem will not cause */
/* > overflow, the Level 2 BLAS routine STBSV is called. If the matrix A */
/* > is singular (A(j,j) = 0 for some j), then s is set to 0 and a */
/* > non-trivial solution to A*x = 0 is returned. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the matrix A is upper or lower triangular. */
/* > = 'U': Upper triangular */
/* > = 'L': Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the operation applied to A. */
/* > = 'N': Solve A * x = s*b (No transpose) */
/* > = 'T': Solve A**T* x = s*b (Transpose) */
/* > = 'C': Solve A**T* x = s*b (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* > DIAG is CHARACTER*1 */
/* > Specifies whether or not the matrix A is unit triangular. */
/* > = 'N': Non-unit triangular */
/* > = 'U': Unit triangular */
/* > \endverbatim */
/* > */
/* > \param[in] NORMIN */
/* > \verbatim */
/* > NORMIN is CHARACTER*1 */
/* > Specifies whether CNORM has been set or not. */
/* > = 'Y': CNORM contains the column norms on entry */
/* > = 'N': CNORM is not set on entry. On exit, the norms will */
/* > be computed and stored in CNORM. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of subdiagonals or superdiagonals in the */
/* > triangular matrix A. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is REAL array, dimension (LDAB,N) */
/* > The upper or lower triangular band matrix A, stored in the */
/* > first KD+1 rows of the array. The j-th column of A is stored */
/* > in the j-th column of the array AB as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j;
*/
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=fla_min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is REAL array, dimension (N) */
/* > On entry, the right hand side b of the triangular system. */
/* > On exit, X is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is REAL */
/* > The scaling factor s for the triangular system */
/* > A * x = s*b or A**T* x = s*b. */
/* > If SCALE = 0, the matrix A is singular or badly scaled, and */
/* > the vector x is an exact or approximate solution to A*x = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CNORM */
/* > \verbatim */
/* > CNORM is REAL array, dimension (N) */
/* > */
/* > If NORMIN = 'Y', CNORM is an input argument and CNORM(j) */
/* > contains the norm of the off-diagonal part of the j-th column */
/* > of A. If TRANS = 'N', CNORM(j) must be greater than or equal */
/* > to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j) */
/* > must be greater than or equal to the 1-norm. */
/* > */
/* > If NORMIN = 'N', CNORM is an output argument and CNORM(j) */
/* > returns the 1-norm of the offdiagonal part of the j-th column */
/* > of A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -k, the k-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > A rough bound on x is computed;
if that is less than overflow, STBSV */
/* > is called, otherwise, specific code is used which checks for possible */
/* > overflow or divide-by-zero at every operation. */
/* > */
/* > A columnwise scheme is used for solving A*x = b. The basic algorithm */
/* > if A is lower triangular is */
/* > */
/* > x[1:n] := b[1:n] */
/* > for j = 1, ..., n */
/* > x(j) := x(j) / A(j,j) */
/* > x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j] */
/* > end */
/* > */
/* > Define bounds on the components of x after j iterations of the loop: */
/* > M(j) = bound on x[1:j] */
/* > G(j) = bound on x[j+1:n] */
/* > Initially, let M(0) = 0 and G(0) = max{
x(i), i=1,...,n}
. */
/* > */
/* > Then for iteration j+1 we have */
/* > M(j+1) <= G(j) / | A(j+1,j+1) | */
/* > G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] | */
/* > <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | ) */
/* > */
/* > where CNORM(j+1) is greater than or equal to the infinity-norm of */
/* > column j+1 of A, not counting the diagonal. Hence */
/* > */
/* > G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | ) */
/* > 1<=i<=j */
/* > and */
/* > */
/* > |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| ) */
/* > 1<=i< j */
/* > */
/* > Since |x(j)| <= M(j), we use the Level 2 BLAS routine STBSV if the */
/* > reciprocal of the largest M(j), j=1,..,n, is larger than */
/* > fla_max(underflow, 1/overflow). */
/* > */
/* > The bound on x(j) is also used to determine when a step in the */
/* > columnwise method can be performed without fear of overflow. If */
/* > the computed bound is greater than a large constant, x is scaled to */
/* > prevent overflow, but if the bound overflows, x is set to 0, x(j) to */
/* > 1, and scale to 0, and a non-trivial solution to A*x = 0 is found. */
/* > */
/* > Similarly, a row-wise scheme is used to solve A**T*x = b. The basic */
/* > algorithm for A upper triangular is */
/* > */
/* > for j = 1, ..., n */
/* > x(j) := ( b(j) - A[1:j-1,j]**T * x[1:j-1] ) / A(j,j) */
/* > end */
/* > */
/* > We simultaneously compute two bounds */
/* > G(j) = bound on ( b(i) - A[1:i-1,i]**T * x[1:i-1] ), 1<=i<=j */
/* > M(j) = bound on x(i), 1<=i<=j */
/* > */
/* > The initial values are G(0) = 0, M(0) = max{
b(i), i=1,..,n}
, and we */
/* > add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1. */
/* > Then the bound on x(j) is */
/* > */
/* > M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) | */
/* > */
/* > <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| ) */
/* > 1<=i<=j */
/* > */
/* > and we can safely call STBSV if 1/M(n) and 1/G(n) are both greater */
/* > than fla_max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int slatbs_(char *uplo, char *trans, char *diag, char * normin, integer *n, integer *kd, real *ab, integer *ldab, real *x, real *scale, real *cnorm, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3;
    /* Local variables */
    integer i__, j;
    real xj, rec, tjj;
    integer jinc, jlen;
    real xbnd;
    integer imax;
    real tmax, tjjs;
    extern real sdot_(integer *, real *, integer *, real *, integer *);
    real xmax, grow, sumj;
    integer maind;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    real tscal, uscal;
    integer jlast;
    extern real sasum_(integer *, real *, integer *);
    logical upper;
    extern /* Subroutine */
    int stbsv_(char *, char *, char *, integer *, integer *, real *, integer *, real *, integer *), saxpy_(integer *, real *, real *, integer *, real *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    real bignum;
    extern integer isamax_(integer *, real *, integer *);
    logical notran;
    integer jfirst;
    real smlnum;
    logical nounit;
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
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --x;
    --cnorm;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    notran = lsame_(trans, "N");
    nounit = lsame_(diag, "N");
    /* Test the input parameters. */
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, "C"))
    {
        *info = -2;
    }
    else if (! nounit && ! lsame_(diag, "U"))
    {
        *info = -3;
    }
    else if (! lsame_(normin, "Y") && ! lsame_(normin, "N"))
    {
        *info = -4;
    }
    else if (*n < 0)
    {
        *info = -5;
    }
    else if (*kd < 0)
    {
        *info = -6;
    }
    else if (*ldab < *kd + 1)
    {
        *info = -8;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLATBS", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Determine machine dependent parameters to control overflow. */
    smlnum = slamch_("Safe minimum") / slamch_("Precision");
    bignum = 1.f / smlnum;
    *scale = 1.f;
    if (lsame_(normin, "N"))
    {
        /* Compute the 1-norm of each column, not including the diagonal. */
        if (upper)
        {
            /* A is upper triangular. */
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Computing MIN */
                i__2 = *kd;
                i__3 = j - 1; // , expr subst
                jlen = fla_min(i__2,i__3);
                cnorm[j] = sasum_(&jlen, &ab[*kd + 1 - jlen + j * ab_dim1], & c__1);
                /* L10: */
            }
        }
        else
        {
            /* A is lower triangular. */
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Computing MIN */
                i__2 = *kd;
                i__3 = *n - j; // , expr subst
                jlen = fla_min(i__2,i__3);
                if (jlen > 0)
                {
                    cnorm[j] = sasum_(&jlen, &ab[j * ab_dim1 + 2], &c__1);
                }
                else
                {
                    cnorm[j] = 0.f;
                }
                /* L20: */
            }
        }
    }
    /* Scale the column norms by TSCAL if the maximum element in CNORM is */
    /* greater than BIGNUM. */
    imax = isamax_(n, &cnorm[1], &c__1);
    tmax = cnorm[imax];
    if (tmax <= bignum)
    {
        tscal = 1.f;
    }
    else
    {
        tscal = 1.f / (smlnum * tmax);
        sscal_(n, &tscal, &cnorm[1], &c__1);
    }
    /* Compute a bound on the computed solution vector to see if the */
    /* Level 2 BLAS routine STBSV can be used. */
    j = isamax_(n, &x[1], &c__1);
    xmax = (r__1 = x[j], f2c_abs(r__1));
    xbnd = xmax;
    if (notran)
    {
        /* Compute the growth in A * x = b. */
        if (upper)
        {
            jfirst = *n;
            jlast = 1;
            jinc = -1;
            maind = *kd + 1;
        }
        else
        {
            jfirst = 1;
            jlast = *n;
            jinc = 1;
            maind = 1;
        }
        if (tscal != 1.f)
        {
            grow = 0.f;
            goto L50;
        }
        if (nounit)
        {
            /* A is non-unit triangular. */
            /* Compute GROW = 1/G(j) and XBND = 1/M(j). */
            /* Initially, G(0) = max{
            x(i), i=1,...,n}
            . */
            grow = 1.f / fla_max(xbnd,smlnum);
            xbnd = grow;
            i__1 = jlast;
            i__2 = jinc;
            for (j = jfirst;
                    i__2 < 0 ? j >= i__1 : j <= i__1;
                    j += i__2)
            {
                /* Exit the loop if the growth factor is too small. */
                if (grow <= smlnum)
                {
                    goto L50;
                }
                /* M(j) = G(j-1) / f2c_abs(A(j,j)) */
                tjj = (r__1 = ab[maind + j * ab_dim1], f2c_abs(r__1));
                /* Computing MIN */
                r__1 = xbnd;
                r__2 = fla_min(1.f,tjj) * grow; // , expr subst
                xbnd = fla_min(r__1,r__2);
                if (tjj + cnorm[j] >= smlnum)
                {
                    /* G(j) = G(j-1)*( 1 + CNORM(j) / f2c_abs(A(j,j)) ) */
                    grow *= tjj / (tjj + cnorm[j]);
                }
                else
                {
                    /* G(j) could overflow, set GROW to 0. */
                    grow = 0.f;
                }
                /* L30: */
            }
            grow = xbnd;
        }
        else
        {
            /* A is unit triangular. */
            /* Compute GROW = 1/G(j), where G(0) = max{
            x(i), i=1,...,n}
            . */
            /* Computing MIN */
            r__1 = 1.f;
            r__2 = 1.f / fla_max(xbnd,smlnum); // , expr subst
            grow = fla_min(r__1,r__2);
            i__2 = jlast;
            i__1 = jinc;
            for (j = jfirst;
                    i__1 < 0 ? j >= i__2 : j <= i__2;
                    j += i__1)
            {
                /* Exit the loop if the growth factor is too small. */
                if (grow <= smlnum)
                {
                    goto L50;
                }
                /* G(j) = G(j-1)*( 1 + CNORM(j) ) */
                grow *= 1.f / (cnorm[j] + 1.f);
                /* L40: */
            }
        }
L50:
        ;
    }
    else
    {
        /* Compute the growth in A**T * x = b. */
        if (upper)
        {
            jfirst = 1;
            jlast = *n;
            jinc = 1;
            maind = *kd + 1;
        }
        else
        {
            jfirst = *n;
            jlast = 1;
            jinc = -1;
            maind = 1;
        }
        if (tscal != 1.f)
        {
            grow = 0.f;
            goto L80;
        }
        if (nounit)
        {
            /* A is non-unit triangular. */
            /* Compute GROW = 1/G(j) and XBND = 1/M(j). */
            /* Initially, M(0) = max{
            x(i), i=1,...,n}
            . */
            grow = 1.f / fla_max(xbnd,smlnum);
            xbnd = grow;
            i__1 = jlast;
            i__2 = jinc;
            for (j = jfirst;
                    i__2 < 0 ? j >= i__1 : j <= i__1;
                    j += i__2)
            {
                /* Exit the loop if the growth factor is too small. */
                if (grow <= smlnum)
                {
                    goto L80;
                }
                /* G(j) = fla_max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */
                xj = cnorm[j] + 1.f;
                /* Computing MIN */
                r__1 = grow;
                r__2 = xbnd / xj; // , expr subst
                grow = fla_min(r__1,r__2);
                /* M(j) = M(j-1)*( 1 + CNORM(j) ) / f2c_abs(A(j,j)) */
                tjj = (r__1 = ab[maind + j * ab_dim1], f2c_abs(r__1));
                if (xj > tjj)
                {
                    xbnd *= tjj / xj;
                }
                /* L60: */
            }
            grow = fla_min(grow,xbnd);
        }
        else
        {
            /* A is unit triangular. */
            /* Compute GROW = 1/G(j), where G(0) = max{
            x(i), i=1,...,n}
            . */
            /* Computing MIN */
            r__1 = 1.f;
            r__2 = 1.f / fla_max(xbnd,smlnum); // , expr subst
            grow = fla_min(r__1,r__2);
            i__2 = jlast;
            i__1 = jinc;
            for (j = jfirst;
                    i__1 < 0 ? j >= i__2 : j <= i__2;
                    j += i__1)
            {
                /* Exit the loop if the growth factor is too small. */
                if (grow <= smlnum)
                {
                    goto L80;
                }
                /* G(j) = ( 1 + CNORM(j) )*G(j-1) */
                xj = cnorm[j] + 1.f;
                grow /= xj;
                /* L70: */
            }
        }
L80:
        ;
    }
    if (grow * tscal > smlnum)
    {
        /* Use the Level 2 BLAS solve if the reciprocal of the bound on */
        /* elements of X is not too small. */
        stbsv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &x[1], &c__1);
    }
    else
    {
        /* Use a Level 1 BLAS solve, scaling intermediate results. */
        if (xmax > bignum)
        {
            /* Scale X so that its components are less than or equal to */
            /* BIGNUM in absolute value. */
            *scale = bignum / xmax;
            sscal_(n, scale, &x[1], &c__1);
            xmax = bignum;
        }
        if (notran)
        {
            /* Solve A * x = b */
            i__1 = jlast;
            i__2 = jinc;
            for (j = jfirst;
                    i__2 < 0 ? j >= i__1 : j <= i__1;
                    j += i__2)
            {
                /* Compute x(j) = b(j) / A(j,j), scaling x if necessary. */
                xj = (r__1 = x[j], f2c_abs(r__1));
                if (nounit)
                {
                    tjjs = ab[maind + j * ab_dim1] * tscal;
                }
                else
                {
                    tjjs = tscal;
                    if (tscal == 1.f)
                    {
                        goto L95;
                    }
                }
                tjj = f2c_abs(tjjs);
                if (tjj > smlnum)
                {
                    /* f2c_abs(A(j,j)) > SMLNUM: */
                    if (tjj < 1.f)
                    {
                        if (xj > tjj * bignum)
                        {
                            /* Scale x by 1/b(j). */
                            rec = 1.f / xj;
                            sscal_(n, &rec, &x[1], &c__1);
                            *scale *= rec;
                            xmax *= rec;
                        }
                    }
                    x[j] /= tjjs;
                    xj = (r__1 = x[j], f2c_abs(r__1));
                }
                else if (tjj > 0.f)
                {
                    /* 0 < f2c_abs(A(j,j)) <= SMLNUM: */
                    if (xj > tjj * bignum)
                    {
                        /* Scale x by (1/f2c_abs(x(j)))*f2c_abs(A(j,j))*BIGNUM */
                        /* to avoid overflow when dividing by A(j,j). */
                        rec = tjj * bignum / xj;
                        if (cnorm[j] > 1.f)
                        {
                            /* Scale by 1/CNORM(j) to avoid overflow when */
                            /* multiplying x(j) times column j. */
                            rec /= cnorm[j];
                        }
                        sscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                        xmax *= rec;
                    }
                    x[j] /= tjjs;
                    xj = (r__1 = x[j], f2c_abs(r__1));
                }
                else
                {
                    /* A(j,j) = 0: Set x(1:n) = 0, x(j) = 1, and */
                    /* scale = 0, and compute a solution to A*x = 0. */
                    i__3 = *n;
                    for (i__ = 1;
                            i__ <= i__3;
                            ++i__)
                    {
                        x[i__] = 0.f;
                        /* L90: */
                    }
                    x[j] = 1.f;
                    xj = 1.f;
                    *scale = 0.f;
                    xmax = 0.f;
                }
L95: /* Scale x if necessary to avoid overflow when adding a */
                /* multiple of column j of A. */
                if (xj > 1.f)
                {
                    rec = 1.f / xj;
                    if (cnorm[j] > (bignum - xmax) * rec)
                    {
                        /* Scale x by 1/(2*f2c_abs(x(j))). */
                        rec *= .5f;
                        sscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                    }
                }
                else if (xj * cnorm[j] > bignum - xmax)
                {
                    /* Scale x by 1/2. */
                    sscal_(n, &c_b36, &x[1], &c__1);
                    *scale *= .5f;
                }
                if (upper)
                {
                    if (j > 1)
                    {
                        /* Compute the update */
                        /* x(fla_max(1,j-kd):j-1) := x(fla_max(1,j-kd):j-1) - */
                        /* x(j)* A(fla_max(1,j-kd):j-1,j) */
                        /* Computing MIN */
                        i__3 = *kd;
                        i__4 = j - 1; // , expr subst
                        jlen = fla_min(i__3,i__4);
                        r__1 = -x[j] * tscal;
                        saxpy_(&jlen, &r__1, &ab[*kd + 1 - jlen + j * ab_dim1], &c__1, &x[j - jlen], &c__1);
                        i__3 = j - 1;
                        i__ = isamax_(&i__3, &x[1], &c__1);
                        xmax = (r__1 = x[i__], f2c_abs(r__1));
                    }
                }
                else if (j < *n)
                {
                    /* Compute the update */
                    /* x(j+1:fla_min(j+kd,n)) := x(j+1:fla_min(j+kd,n)) - */
                    /* x(j) * A(j+1:fla_min(j+kd,n),j) */
                    /* Computing MIN */
                    i__3 = *kd;
                    i__4 = *n - j; // , expr subst
                    jlen = fla_min(i__3,i__4);
                    if (jlen > 0)
                    {
                        r__1 = -x[j] * tscal;
                        saxpy_(&jlen, &r__1, &ab[j * ab_dim1 + 2], &c__1, &x[ j + 1], &c__1);
                    }
                    i__3 = *n - j;
                    i__ = j + isamax_(&i__3, &x[j + 1], &c__1);
                    xmax = (r__1 = x[i__], f2c_abs(r__1));
                }
                /* L100: */
            }
        }
        else
        {
            /* Solve A**T * x = b */
            i__2 = jlast;
            i__1 = jinc;
            for (j = jfirst;
                    i__1 < 0 ? j >= i__2 : j <= i__2;
                    j += i__1)
            {
                /* Compute x(j) = b(j) - sum A(k,j)*x(k). */
                /* k<>j */
                xj = (r__1 = x[j], f2c_abs(r__1));
                uscal = tscal;
                rec = 1.f / fla_max(xmax,1.f);
                if (cnorm[j] > (bignum - xj) * rec)
                {
                    /* If x(j) could overflow, scale x by 1/(2*XMAX). */
                    rec *= .5f;
                    if (nounit)
                    {
                        tjjs = ab[maind + j * ab_dim1] * tscal;
                    }
                    else
                    {
                        tjjs = tscal;
                    }
                    tjj = f2c_abs(tjjs);
                    if (tjj > 1.f)
                    {
                        /* Divide by A(j,j) when scaling x if A(j,j) > 1. */
                        /* Computing MIN */
                        r__1 = 1.f;
                        r__2 = rec * tjj; // , expr subst
                        rec = fla_min(r__1,r__2);
                        uscal /= tjjs;
                    }
                    if (rec < 1.f)
                    {
                        sscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                        xmax *= rec;
                    }
                }
                sumj = 0.f;
                if (uscal == 1.f)
                {
                    /* If the scaling needed for A in the dot product is 1, */
                    /* call SDOT to perform the dot product. */
                    if (upper)
                    {
                        /* Computing MIN */
                        i__3 = *kd;
                        i__4 = j - 1; // , expr subst
                        jlen = fla_min(i__3,i__4);
                        sumj = sdot_(&jlen, &ab[*kd + 1 - jlen + j * ab_dim1], &c__1, &x[j - jlen], &c__1);
                    }
                    else
                    {
                        /* Computing MIN */
                        i__3 = *kd;
                        i__4 = *n - j; // , expr subst
                        jlen = fla_min(i__3,i__4);
                        if (jlen > 0)
                        {
                            sumj = sdot_(&jlen, &ab[j * ab_dim1 + 2], &c__1, & x[j + 1], &c__1);
                        }
                    }
                }
                else
                {
                    /* Otherwise, use in-line code for the dot product. */
                    if (upper)
                    {
                        /* Computing MIN */
                        i__3 = *kd;
                        i__4 = j - 1; // , expr subst
                        jlen = fla_min(i__3,i__4);
                        i__3 = jlen;
                        for (i__ = 1;
                                i__ <= i__3;
                                ++i__)
                        {
                            sumj += ab[*kd + i__ - jlen + j * ab_dim1] * uscal * x[j - jlen - 1 + i__];
                            /* L110: */
                        }
                    }
                    else
                    {
                        /* Computing MIN */
                        i__3 = *kd;
                        i__4 = *n - j; // , expr subst
                        jlen = fla_min(i__3,i__4);
                        i__3 = jlen;
                        for (i__ = 1;
                                i__ <= i__3;
                                ++i__)
                        {
                            sumj += ab[i__ + 1 + j * ab_dim1] * uscal * x[j + i__];
                            /* L120: */
                        }
                    }
                }
                if (uscal == tscal)
                {
                    /* Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j) */
                    /* was not used to scale the dotproduct. */
                    x[j] -= sumj;
                    xj = (r__1 = x[j], f2c_abs(r__1));
                    if (nounit)
                    {
                        /* Compute x(j) = x(j) / A(j,j), scaling if necessary. */
                        tjjs = ab[maind + j * ab_dim1] * tscal;
                    }
                    else
                    {
                        tjjs = tscal;
                        if (tscal == 1.f)
                        {
                            goto L135;
                        }
                    }
                    tjj = f2c_abs(tjjs);
                    if (tjj > smlnum)
                    {
                        /* f2c_abs(A(j,j)) > SMLNUM: */
                        if (tjj < 1.f)
                        {
                            if (xj > tjj * bignum)
                            {
                                /* Scale X by 1/f2c_abs(x(j)). */
                                rec = 1.f / xj;
                                sscal_(n, &rec, &x[1], &c__1);
                                *scale *= rec;
                                xmax *= rec;
                            }
                        }
                        x[j] /= tjjs;
                    }
                    else if (tjj > 0.f)
                    {
                        /* 0 < f2c_abs(A(j,j)) <= SMLNUM: */
                        if (xj > tjj * bignum)
                        {
                            /* Scale x by (1/f2c_abs(x(j)))*f2c_abs(A(j,j))*BIGNUM. */
                            rec = tjj * bignum / xj;
                            sscal_(n, &rec, &x[1], &c__1);
                            *scale *= rec;
                            xmax *= rec;
                        }
                        x[j] /= tjjs;
                    }
                    else
                    {
                        /* A(j,j) = 0: Set x(1:n) = 0, x(j) = 1, and */
                        /* scale = 0, and compute a solution to A**T*x = 0. */
                        i__3 = *n;
                        for (i__ = 1;
                                i__ <= i__3;
                                ++i__)
                        {
                            x[i__] = 0.f;
                            /* L130: */
                        }
                        x[j] = 1.f;
                        *scale = 0.f;
                        xmax = 0.f;
                    }
L135:
                    ;
                }
                else
                {
                    /* Compute x(j) := x(j) / A(j,j) - sumj if the dot */
                    /* product has already been divided by 1/A(j,j). */
                    x[j] = x[j] / tjjs - sumj;
                }
                /* Computing MAX */
                r__2 = xmax;
                r__3 = (r__1 = x[j], f2c_abs(r__1)); // , expr subst
                xmax = fla_max(r__2,r__3);
                /* L140: */
            }
        }
        *scale /= tscal;
    }
    /* Scale the column norms by 1/TSCAL for return. */
    if (tscal != 1.f)
    {
        r__1 = 1.f / tscal;
        sscal_(n, &r__1, &cnorm[1], &c__1);
    }
    return 0;
    /* End of SLATBS */
}
/* slatbs_ */
