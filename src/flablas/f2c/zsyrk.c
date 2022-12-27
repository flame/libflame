/* zsyrk.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zsyrk_(char *uplo, char *trans, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex * beta, doublecomplex *c__, integer *ldc)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1, z__2, z__3;
    /* Local variables */
    integer info;
    doublecomplex temp;
    integer i__, j, l;
    extern logical lsame_(char *, char *);
    integer nrowa;
    logical upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* ZSYRK performs one of the symmetric rank k operations */
    /* C := alpha*A*A' + beta*C, */
    /* or */
    /* C := alpha*A'*A + beta*C, */
    /* where alpha and beta are scalars, C is an n by n symmetric matrix */
    /* and A is an n by k matrix in the first case and a k by n matrix */
    /* in the second case. */
    /* Parameters */
    /* ========== */
    /* UPLO - CHARACTER*1. */
    /* On entry, UPLO specifies whether the upper or lower */
    /* triangular part of the array C is to be referenced as */
    /* follows: */
    /* UPLO = 'U' or 'u' Only the upper triangular part of C */
    /* is to be referenced. */
    /* UPLO = 'L' or 'l' Only the lower triangular part of C */
    /* is to be referenced. */
    /* Unchanged on exit. */
    /* TRANS - CHARACTER*1. */
    /* On entry, TRANS specifies the operation to be performed as */
    /* follows: */
    /* TRANS = 'N' or 'n' C := alpha*A*A' + beta*C. */
    /* TRANS = 'T' or 't' C := alpha*A'*A + beta*C. */
    /* Unchanged on exit. */
    /* N - INTEGER. */
    /* On entry, N specifies the order of the matrix C. N must be */
    /* at least zero. */
    /* Unchanged on exit. */
    /* K - INTEGER. */
    /* On entry with TRANS = 'N' or 'n', K specifies the number */
    /* of columns of the matrix A, and on entry with */
    /* TRANS = 'T' or 't', K specifies the number of rows of the */
    /* matrix A. K must be at least zero. */
    /* Unchanged on exit. */
    /* ALPHA - COMPLEX*16 . */
    /* On entry, ALPHA specifies the scalar alpha. */
    /* Unchanged on exit. */
    /* A - COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is */
    /* k when TRANS = 'N' or 'n', and is n otherwise. */
    /* Before entry with TRANS = 'N' or 'n', the leading n by k */
    /* part of the array A must contain the matrix A, otherwise */
    /* the leading k by n part of the array A must contain the */
    /* matrix A. */
    /* Unchanged on exit. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. When TRANS = 'N' or 'n' */
    /* then LDA must be at least fla_max( 1, n ), otherwise LDA must */
    /* be at least fla_max( 1, k ). */
    /* Unchanged on exit. */
    /* BETA - COMPLEX*16 . */
    /* On entry, BETA specifies the scalar beta. */
    /* Unchanged on exit. */
    /* C - COMPLEX*16 array of DIMENSION ( LDC, n ). */
    /* Before entry with UPLO = 'U' or 'u', the leading n by n */
    /* upper triangular part of the array C must contain the upper */
    /* triangular part of the symmetric matrix and the strictly */
    /* lower triangular part of C is not referenced. On exit, the */
    /* upper triangular part of the array C is overwritten by the */
    /* upper triangular part of the updated matrix. */
    /* Before entry with UPLO = 'L' or 'l', the leading n by n */
    /* lower triangular part of the array C must contain the lower */
    /* triangular part of the symmetric matrix and the strictly */
    /* upper triangular part of C is not referenced. On exit, the */
    /* lower triangular part of the array C is overwritten by the */
    /* lower triangular part of the updated matrix. */
    /* LDC - INTEGER. */
    /* On entry, LDC specifies the first dimension of C as declared */
    /* in the calling (sub) program. LDC must be at least */
    /* fla_max( 1, n ). */
    /* Unchanged on exit. */
    /* Level 3 Blas routine. */
    /* -- Written on 8-February-1989. */
    /* Jack Dongarra, Argonne National Laboratory. */
    /* Iain Duff, AERE Harwell. */
    /* Jeremy Du Croz, Numerical Algorithms Group Ltd. */
    /* Sven Hammarling, Numerical Algorithms Group Ltd. */
    /* .. External Functions .. */
    /* .. External Subroutines .. */
    /* .. Intrinsic Functions .. */
    /* .. Local Scalars .. */
    /* .. Parameters .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    /* Function Body */
    if (lsame_(trans, "N"))
    {
        nrowa = *n;
    }
    else
    {
        nrowa = *k;
    }
    upper = lsame_(uplo, "U");
    info = 0;
    if (! upper && ! lsame_(uplo, "L"))
    {
        info = 1;
    }
    else if (! lsame_(trans, "N") && ! lsame_(trans, "T"))
    {
        info = 2;
    }
    else if (*n < 0)
    {
        info = 3;
    }
    else if (*k < 0)
    {
        info = 4;
    }
    else if (*lda < fla_max(1,nrowa))
    {
        info = 7;
    }
    else if (*ldc < fla_max(1,*n))
    {
        info = 10;
    }
    if (info != 0)
    {
        xerbla_("ZSYRK ", &info);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0 || (alpha->r == 0. && alpha->i == 0. || *k == 0) && (beta->r == 1. && beta->i == 0.))
    {
        return 0;
    }
    /* And when alpha.eq.zero. */
    if (alpha->r == 0. && alpha->i == 0.)
    {
        if (upper)
        {
            if (beta->r == 0. && beta->i == 0.)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = j;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                        /* L10: */
                    }
                    /* L20: */
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = j;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4] .i, z__1.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                        /* L30: */
                    }
                    /* L40: */
                }
            }
        }
        else
        {
            if (beta->r == 0. && beta->i == 0.)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *n;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                        /* L50: */
                    }
                    /* L60: */
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *n;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4] .i, z__1.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                        /* L70: */
                    }
                    /* L80: */
                }
            }
        }
        return 0;
    }
    /* Start the operations. */
    if (lsame_(trans, "N"))
    {
        /* Form C := alpha*A*A' + beta*C. */
        if (upper)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (beta->r == 0. && beta->i == 0.)
                {
                    i__2 = j;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                        /* L90: */
                    }
                }
                else if (beta->r != 1. || beta->i != 0.)
                {
                    i__2 = j;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4] .i, z__1.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                        /* L100: */
                    }
                }
                i__2 = *k;
                for (l = 1;
                        l <= i__2;
                        ++l)
                {
                    i__3 = j + l * a_dim1;
                    if (a[i__3].r != 0. || a[i__3].i != 0.)
                    {
                        i__3 = j + l * a_dim1;
                        z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, z__1.i = alpha->r * a[i__3].i + alpha->i * a[ i__3].r;
                        temp.r = z__1.r, temp.i = z__1.i;
                        i__3 = j;
                        for (i__ = 1;
                                i__ <= i__3;
                                ++i__)
                        {
                            i__4 = i__ + j * c_dim1;
                            i__5 = i__ + j * c_dim1;
                            i__6 = i__ + l * a_dim1;
                            z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, z__2.i = temp.r * a[i__6].i + temp.i * a[ i__6].r;
                            z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5] .i + z__2.i;
                            c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
                            /* L110: */
                        }
                    }
                    /* L120: */
                }
                /* L130: */
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (beta->r == 0. && beta->i == 0.)
                {
                    i__2 = *n;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                        /* L140: */
                    }
                }
                else if (beta->r != 1. || beta->i != 0.)
                {
                    i__2 = *n;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4] .i, z__1.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                        /* L150: */
                    }
                }
                i__2 = *k;
                for (l = 1;
                        l <= i__2;
                        ++l)
                {
                    i__3 = j + l * a_dim1;
                    if (a[i__3].r != 0. || a[i__3].i != 0.)
                    {
                        i__3 = j + l * a_dim1;
                        z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, z__1.i = alpha->r * a[i__3].i + alpha->i * a[ i__3].r;
                        temp.r = z__1.r, temp.i = z__1.i;
                        i__3 = *n;
                        for (i__ = j;
                                i__ <= i__3;
                                ++i__)
                        {
                            i__4 = i__ + j * c_dim1;
                            i__5 = i__ + j * c_dim1;
                            i__6 = i__ + l * a_dim1;
                            z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, z__2.i = temp.r * a[i__6].i + temp.i * a[ i__6].r;
                            z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5] .i + z__2.i;
                            c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
                            /* L160: */
                        }
                    }
                    /* L170: */
                }
                /* L180: */
            }
        }
    }
    else
    {
        /* Form C := alpha*A'*A + beta*C. */
        if (upper)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    temp.r = 0., temp.i = 0.;
                    i__3 = *k;
                    for (l = 1;
                            l <= i__3;
                            ++l)
                    {
                        i__4 = l + i__ * a_dim1;
                        i__5 = l + j * a_dim1;
                        z__2.r = a[i__4].r * a[i__5].r - a[i__4].i * a[i__5] .i, z__2.i = a[i__4].r * a[i__5].i + a[i__4] .i * a[i__5].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                        /* L190: */
                    }
                    if (beta->r == 0. && beta->i == 0.)
                    {
                        i__3 = i__ + j * c_dim1;
                        z__1.r = alpha->r * temp.r - alpha->i * temp.i, z__1.i = alpha->r * temp.i + alpha->i * temp.r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                    else
                    {
                        i__3 = i__ + j * c_dim1;
                        z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = alpha->r * temp.i + alpha->i * temp.r;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4] .i, z__3.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                    /* L200: */
                }
                /* L210: */
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = *n;
                for (i__ = j;
                        i__ <= i__2;
                        ++i__)
                {
                    temp.r = 0., temp.i = 0.;
                    i__3 = *k;
                    for (l = 1;
                            l <= i__3;
                            ++l)
                    {
                        i__4 = l + i__ * a_dim1;
                        i__5 = l + j * a_dim1;
                        z__2.r = a[i__4].r * a[i__5].r - a[i__4].i * a[i__5] .i, z__2.i = a[i__4].r * a[i__5].i + a[i__4] .i * a[i__5].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                        /* L220: */
                    }
                    if (beta->r == 0. && beta->i == 0.)
                    {
                        i__3 = i__ + j * c_dim1;
                        z__1.r = alpha->r * temp.r - alpha->i * temp.i, z__1.i = alpha->r * temp.i + alpha->i * temp.r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                    else
                    {
                        i__3 = i__ + j * c_dim1;
                        z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = alpha->r * temp.i + alpha->i * temp.r;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4] .i, z__3.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                    /* L230: */
                }
                /* L240: */
            }
        }
    }
    return 0;
    /* End of ZSYRK . */
}
/* zsyrk_ */

