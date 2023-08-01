#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int ctrtri_check(char *uplo, char *diag, integer *n, scomplex *a, integer *lda, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    logical upper;
    logical nounit;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    nounit = lsame_(diag, "N");
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (! nounit && ! lsame_(diag, "U"))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*lda < fla_max(1,*n))
    {
        *info = -5;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CTRTRI", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    /* Check for singularity if non-unit. */
    if (nounit)
    {
        i__1 = *n;
        for (*info = 1;
                *info <= i__1;
                ++(*info))
        {
            i__2 = *info + *info * a_dim1;
            if (a[i__2].real == 0.f && a[i__2].imag == 0.f)
            {
                return LAPACK_FAILURE;
            }
            /* L10: */
        }
        *info = 0;
    }
    return LAPACK_SUCCESS;
}
