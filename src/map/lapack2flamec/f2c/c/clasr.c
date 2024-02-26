/* ../netlib/clasr.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLASR applies a sequence of plane rotations to a general rectangular matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLASR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clasr.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clasr.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clasr.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIRECT, PIVOT, SIDE */
/* INTEGER LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL C( * ), S( * ) */
/* COMPLEX A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLASR applies a sequence of real plane rotations to a complex matrix */
/* > A, from either the left or the right. */
/* > */
/* > When SIDE = 'L', the transformation takes the form */
/* > */
/* > A := P*A */
/* > */
/* > and when SIDE = 'R', the transformation takes the form */
/* > */
/* > A := A*P**T */
/* > */
/* > where P is an orthogonal matrix consisting of a sequence of z plane */
/* > rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R', */
/* > and P**T is the transpose of P. */
/* > */
/* > When DIRECT = 'F' (Forward sequence), then */
/* > */
/* > P = P(z-1) * ... * P(2) * P(1) */
/* > */
/* > and when DIRECT = 'B' (Backward sequence), then */
/* > */
/* > P = P(1) * P(2) * ... * P(z-1) */
/* > */
/* > where P(k) is a plane rotation matrix defined by the 2-by-2 rotation */
/* > */
/* > R(k) = ( c(k) s(k) ) */
/* > = ( -s(k) c(k) ). */
/* > */
/* > When PIVOT = 'V' (Variable pivot), the rotation is performed */
/* > for the plane (k,k+1), i.e., P(k) has the form */
/* > */
/* > P(k) = ( 1 ) */
/* > ( ... ) */
/* > ( 1 ) */
/* > ( c(k) s(k) ) */
/* > ( -s(k) c(k) ) */
/* > ( 1 ) */
/* > ( ... ) */
/* > ( 1 ) */
/* > */
/* > where R(k) appears as a rank-2 modification to the identity matrix in */
/* > rows and columns k and k+1. */
/* > */
/* > When PIVOT = 'T' (Top pivot), the rotation is performed for the */
/* > plane (1,k+1), so P(k) has the form */
/* > */
/* > P(k) = ( c(k) s(k) ) */
/* > ( 1 ) */
/* > ( ... ) */
/* > ( 1 ) */
/* > ( -s(k) c(k) ) */
/* > ( 1 ) */
/* > ( ... ) */
/* > ( 1 ) */
/* > */
/* > where R(k) appears in rows and columns 1 and k+1. */
/* > */
/* > Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is */
/* > performed for the plane (k,z), giving P(k) the form */
/* > */
/* > P(k) = ( 1 ) */
/* > ( ... ) */
/* > ( 1 ) */
/* > ( c(k) s(k) ) */
/* > ( 1 ) */
/* > ( ... ) */
/* > ( 1 ) */
/* > ( -s(k) c(k) ) */
/* > */
/* > where R(k) appears in rows and columns k and z. The rotations are */
/* > performed without ever forming P(k) explicitly. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > Specifies whether the plane rotation matrix P is applied to */
/* > A on the left or the right. */
/* > = 'L': Left, compute A := P*A */
/* > = 'R': Right, compute A:= A*P**T */
/* > \endverbatim */
/* > */
/* > \param[in] PIVOT */
/* > \verbatim */
/* > PIVOT is CHARACTER*1 */
/* > Specifies the plane for which P(k) is a plane rotation */
/* > matrix. */
/* > = 'V': Variable pivot, the plane (k,k+1) */
/* > = 'T': Top pivot, the plane (1,k+1) */
/* > = 'B': Bottom pivot, the plane (k,z) */
/* > \endverbatim */
/* > */
/* > \param[in] DIRECT */
/* > \verbatim */
/* > DIRECT is CHARACTER*1 */
/* > Specifies whether P is a forward or backward sequence of */
/* > plane rotations. */
/* > = 'F': Forward, P = P(z-1)*...*P(2)*P(1) */
/* > = 'B': Backward, P = P(1)*P(2)*...*P(z-1) */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. If m <= 1, an immediate */
/* > return is effected. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. If n <= 1, an */
/* > immediate return is effected. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is REAL array, dimension */
/* > (M-1) if SIDE = 'L' */
/* > (N-1) if SIDE = 'R' */
/* > The cosines c(k) of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is REAL array, dimension */
/* > (M-1) if SIDE = 'L' */
/* > (N-1) if SIDE = 'R' */
/* > The sines s(k) of the plane rotations. The 2-by-2 plane */
/* > rotation part of the matrix P(k), R(k), has the form */
/* > R(k) = ( c(k) s(k) ) */
/* > ( -s(k) c(k) ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > The M-by-N matrix A. On exit, A is overwritten by P*A if */
/* > SIDE = 'R' or by A*P**T if SIDE = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
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
int clasr_(char *side, char *pivot, char *direct, integer *m, integer *n, real *c__, real *s, complex *a, integer *lda)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"clasr inputs: side %c, pivot %c, direct %c, m %lld, n %lld, lda %lld",*side, *pivot, *direct, *m, *n, *lda);
#else
    snprintf(buffer, 256,"clasr inputs: side %c, pivot %c, direct %c, m %d, n %d, lda %d",*side, *pivot, *direct, *m, *n, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    complex q__1, q__2, q__3;
    /* Local variables */
    integer i__, j, info;
    complex temp;
    extern logical lsame_(char *, char *);
    real ctemp, stemp;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    --c__;
    --s;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    /* Function Body */
    info = 0;
    if (! (lsame_(side, "L") || lsame_(side, "R")))
    {
        info = 1;
    }
    else if (! (lsame_(pivot, "V") || lsame_(pivot, "T") || lsame_(pivot, "B")))
    {
        info = 2;
    }
    else if (! (lsame_(direct, "F") || lsame_(direct, "B")))
    {
        info = 3;
    }
    else if (*m < 0)
    {
        info = 4;
    }
    else if (*n < 0)
    {
        info = 5;
    }
    else if (*lda < fla_max(1,*m))
    {
        info = 9;
    }
    if (info != 0)
    {
        xerbla_("CLASR ", &info, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    if (lsame_(side, "L"))
    {
        /* Form P * A */
        if (lsame_(pivot, "V"))
        {
            if (lsame_(direct, "F"))
            {
                i__1 = *m - 1;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__2 = *n;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = j + 1 + i__ * a_dim1;
                            temp.r = a[i__3].r;
                            temp.i = a[i__3].i; // , expr subst
                            i__3 = j + 1 + i__ * a_dim1;
                            q__2.r = ctemp * temp.r;
                            q__2.i = ctemp * temp.i; // , expr subst
                            i__4 = j + i__ * a_dim1;
                            q__3.r = stemp * a[i__4].r;
                            q__3.i = stemp * a[ i__4].i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            i__3 = j + i__ * a_dim1;
                            q__2.r = stemp * temp.r;
                            q__2.i = stemp * temp.i; // , expr subst
                            i__4 = j + i__ * a_dim1;
                            q__3.r = ctemp * a[i__4].r;
                            q__3.i = ctemp * a[ i__4].i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            /* L10: */
                        }
                    }
                    /* L20: */
                }
            }
            else if (lsame_(direct, "B"))
            {
                for (j = *m - 1;
                        j >= 1;
                        --j)
                {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__1 = *n;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = j + 1 + i__ * a_dim1;
                            temp.r = a[i__2].r;
                            temp.i = a[i__2].i; // , expr subst
                            i__2 = j + 1 + i__ * a_dim1;
                            q__2.r = ctemp * temp.r;
                            q__2.i = ctemp * temp.i; // , expr subst
                            i__3 = j + i__ * a_dim1;
                            q__3.r = stemp * a[i__3].r;
                            q__3.i = stemp * a[ i__3].i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            i__2 = j + i__ * a_dim1;
                            q__2.r = stemp * temp.r;
                            q__2.i = stemp * temp.i; // , expr subst
                            i__3 = j + i__ * a_dim1;
                            q__3.r = ctemp * a[i__3].r;
                            q__3.i = ctemp * a[ i__3].i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            /* L30: */
                        }
                    }
                    /* L40: */
                }
            }
        }
        else if (lsame_(pivot, "T"))
        {
            if (lsame_(direct, "F"))
            {
                i__1 = *m;
                for (j = 2;
                        j <= i__1;
                        ++j)
                {
                    ctemp = c__[j - 1];
                    stemp = s[j - 1];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__2 = *n;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = j + i__ * a_dim1;
                            temp.r = a[i__3].r;
                            temp.i = a[i__3].i; // , expr subst
                            i__3 = j + i__ * a_dim1;
                            q__2.r = ctemp * temp.r;
                            q__2.i = ctemp * temp.i; // , expr subst
                            i__4 = i__ * a_dim1 + 1;
                            q__3.r = stemp * a[i__4].r;
                            q__3.i = stemp * a[ i__4].i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            i__3 = i__ * a_dim1 + 1;
                            q__2.r = stemp * temp.r;
                            q__2.i = stemp * temp.i; // , expr subst
                            i__4 = i__ * a_dim1 + 1;
                            q__3.r = ctemp * a[i__4].r;
                            q__3.i = ctemp * a[ i__4].i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            /* L50: */
                        }
                    }
                    /* L60: */
                }
            }
            else if (lsame_(direct, "B"))
            {
                for (j = *m;
                        j >= 2;
                        --j)
                {
                    ctemp = c__[j - 1];
                    stemp = s[j - 1];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__1 = *n;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = j + i__ * a_dim1;
                            temp.r = a[i__2].r;
                            temp.i = a[i__2].i; // , expr subst
                            i__2 = j + i__ * a_dim1;
                            q__2.r = ctemp * temp.r;
                            q__2.i = ctemp * temp.i; // , expr subst
                            i__3 = i__ * a_dim1 + 1;
                            q__3.r = stemp * a[i__3].r;
                            q__3.i = stemp * a[ i__3].i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            i__2 = i__ * a_dim1 + 1;
                            q__2.r = stemp * temp.r;
                            q__2.i = stemp * temp.i; // , expr subst
                            i__3 = i__ * a_dim1 + 1;
                            q__3.r = ctemp * a[i__3].r;
                            q__3.i = ctemp * a[ i__3].i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            /* L70: */
                        }
                    }
                    /* L80: */
                }
            }
        }
        else if (lsame_(pivot, "B"))
        {
            if (lsame_(direct, "F"))
            {
                i__1 = *m - 1;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__2 = *n;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = j + i__ * a_dim1;
                            temp.r = a[i__3].r;
                            temp.i = a[i__3].i; // , expr subst
                            i__3 = j + i__ * a_dim1;
                            i__4 = *m + i__ * a_dim1;
                            q__2.r = stemp * a[i__4].r;
                            q__2.i = stemp * a[ i__4].i; // , expr subst
                            q__3.r = ctemp * temp.r;
                            q__3.i = ctemp * temp.i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            i__3 = *m + i__ * a_dim1;
                            i__4 = *m + i__ * a_dim1;
                            q__2.r = ctemp * a[i__4].r;
                            q__2.i = ctemp * a[ i__4].i; // , expr subst
                            q__3.r = stemp * temp.r;
                            q__3.i = stemp * temp.i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            /* L90: */
                        }
                    }
                    /* L100: */
                }
            }
            else if (lsame_(direct, "B"))
            {
                for (j = *m - 1;
                        j >= 1;
                        --j)
                {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__1 = *n;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = j + i__ * a_dim1;
                            temp.r = a[i__2].r;
                            temp.i = a[i__2].i; // , expr subst
                            i__2 = j + i__ * a_dim1;
                            i__3 = *m + i__ * a_dim1;
                            q__2.r = stemp * a[i__3].r;
                            q__2.i = stemp * a[ i__3].i; // , expr subst
                            q__3.r = ctemp * temp.r;
                            q__3.i = ctemp * temp.i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            i__2 = *m + i__ * a_dim1;
                            i__3 = *m + i__ * a_dim1;
                            q__2.r = ctemp * a[i__3].r;
                            q__2.i = ctemp * a[ i__3].i; // , expr subst
                            q__3.r = stemp * temp.r;
                            q__3.i = stemp * temp.i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            /* L110: */
                        }
                    }
                    /* L120: */
                }
            }
        }
    }
    else if (lsame_(side, "R"))
    {
        /* Form A * P**T */
        if (lsame_(pivot, "V"))
        {
            if (lsame_(direct, "F"))
            {
                i__1 = *n - 1;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__2 = *m;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + (j + 1) * a_dim1;
                            temp.r = a[i__3].r;
                            temp.i = a[i__3].i; // , expr subst
                            i__3 = i__ + (j + 1) * a_dim1;
                            q__2.r = ctemp * temp.r;
                            q__2.i = ctemp * temp.i; // , expr subst
                            i__4 = i__ + j * a_dim1;
                            q__3.r = stemp * a[i__4].r;
                            q__3.i = stemp * a[ i__4].i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            i__3 = i__ + j * a_dim1;
                            q__2.r = stemp * temp.r;
                            q__2.i = stemp * temp.i; // , expr subst
                            i__4 = i__ + j * a_dim1;
                            q__3.r = ctemp * a[i__4].r;
                            q__3.i = ctemp * a[ i__4].i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            /* L130: */
                        }
                    }
                    /* L140: */
                }
            }
            else if (lsame_(direct, "B"))
            {
                for (j = *n - 1;
                        j >= 1;
                        --j)
                {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__1 = *m;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = i__ + (j + 1) * a_dim1;
                            temp.r = a[i__2].r;
                            temp.i = a[i__2].i; // , expr subst
                            i__2 = i__ + (j + 1) * a_dim1;
                            q__2.r = ctemp * temp.r;
                            q__2.i = ctemp * temp.i; // , expr subst
                            i__3 = i__ + j * a_dim1;
                            q__3.r = stemp * a[i__3].r;
                            q__3.i = stemp * a[ i__3].i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            i__2 = i__ + j * a_dim1;
                            q__2.r = stemp * temp.r;
                            q__2.i = stemp * temp.i; // , expr subst
                            i__3 = i__ + j * a_dim1;
                            q__3.r = ctemp * a[i__3].r;
                            q__3.i = ctemp * a[ i__3].i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            /* L150: */
                        }
                    }
                    /* L160: */
                }
            }
        }
        else if (lsame_(pivot, "T"))
        {
            if (lsame_(direct, "F"))
            {
                i__1 = *n;
                for (j = 2;
                        j <= i__1;
                        ++j)
                {
                    ctemp = c__[j - 1];
                    stemp = s[j - 1];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__2 = *m;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + j * a_dim1;
                            temp.r = a[i__3].r;
                            temp.i = a[i__3].i; // , expr subst
                            i__3 = i__ + j * a_dim1;
                            q__2.r = ctemp * temp.r;
                            q__2.i = ctemp * temp.i; // , expr subst
                            i__4 = i__ + a_dim1;
                            q__3.r = stemp * a[i__4].r;
                            q__3.i = stemp * a[ i__4].i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            i__3 = i__ + a_dim1;
                            q__2.r = stemp * temp.r;
                            q__2.i = stemp * temp.i; // , expr subst
                            i__4 = i__ + a_dim1;
                            q__3.r = ctemp * a[i__4].r;
                            q__3.i = ctemp * a[ i__4].i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            /* L170: */
                        }
                    }
                    /* L180: */
                }
            }
            else if (lsame_(direct, "B"))
            {
                for (j = *n;
                        j >= 2;
                        --j)
                {
                    ctemp = c__[j - 1];
                    stemp = s[j - 1];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__1 = *m;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = i__ + j * a_dim1;
                            temp.r = a[i__2].r;
                            temp.i = a[i__2].i; // , expr subst
                            i__2 = i__ + j * a_dim1;
                            q__2.r = ctemp * temp.r;
                            q__2.i = ctemp * temp.i; // , expr subst
                            i__3 = i__ + a_dim1;
                            q__3.r = stemp * a[i__3].r;
                            q__3.i = stemp * a[ i__3].i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            i__2 = i__ + a_dim1;
                            q__2.r = stemp * temp.r;
                            q__2.i = stemp * temp.i; // , expr subst
                            i__3 = i__ + a_dim1;
                            q__3.r = ctemp * a[i__3].r;
                            q__3.i = ctemp * a[ i__3].i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            /* L190: */
                        }
                    }
                    /* L200: */
                }
            }
        }
        else if (lsame_(pivot, "B"))
        {
            if (lsame_(direct, "F"))
            {
                i__1 = *n - 1;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__2 = *m;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + j * a_dim1;
                            temp.r = a[i__3].r;
                            temp.i = a[i__3].i; // , expr subst
                            i__3 = i__ + j * a_dim1;
                            i__4 = i__ + *n * a_dim1;
                            q__2.r = stemp * a[i__4].r;
                            q__2.i = stemp * a[ i__4].i; // , expr subst
                            q__3.r = ctemp * temp.r;
                            q__3.i = ctemp * temp.i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            i__3 = i__ + *n * a_dim1;
                            i__4 = i__ + *n * a_dim1;
                            q__2.r = ctemp * a[i__4].r;
                            q__2.i = ctemp * a[ i__4].i; // , expr subst
                            q__3.r = stemp * temp.r;
                            q__3.i = stemp * temp.i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__3].r = q__1.r;
                            a[i__3].i = q__1.i; // , expr subst
                            /* L210: */
                        }
                    }
                    /* L220: */
                }
            }
            else if (lsame_(direct, "B"))
            {
                for (j = *n - 1;
                        j >= 1;
                        --j)
                {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1.f || stemp != 0.f)
                    {
                        i__1 = *m;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = i__ + j * a_dim1;
                            temp.r = a[i__2].r;
                            temp.i = a[i__2].i; // , expr subst
                            i__2 = i__ + j * a_dim1;
                            i__3 = i__ + *n * a_dim1;
                            q__2.r = stemp * a[i__3].r;
                            q__2.i = stemp * a[ i__3].i; // , expr subst
                            q__3.r = ctemp * temp.r;
                            q__3.i = ctemp * temp.i; // , expr subst
                            q__1.r = q__2.r + q__3.r;
                            q__1.i = q__2.i + q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            i__2 = i__ + *n * a_dim1;
                            i__3 = i__ + *n * a_dim1;
                            q__2.r = ctemp * a[i__3].r;
                            q__2.i = ctemp * a[ i__3].i; // , expr subst
                            q__3.r = stemp * temp.r;
                            q__3.i = stemp * temp.i; // , expr subst
                            q__1.r = q__2.r - q__3.r;
                            q__1.i = q__2.i - q__3.i; // , expr subst
                            a[i__2].r = q__1.r;
                            a[i__2].i = q__1.i; // , expr subst
                            /* L230: */
                        }
                    }
                    /* L240: */
                }
            }
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CLASR */
}
/* clasr_ */
