#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int dtrtri_check(char *uplo, char *diag, integer *n, double * a, integer *lda, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
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
        xerbla_("DTRTRI", &i__1);
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
            if (a[*info + *info * a_dim1] == 0.)
            {
                return LAPACK_FAILURE;
            }
            /* L10: */
        }
        *info = 0;
    }

    return LAPACK_SUCCESS;
}
