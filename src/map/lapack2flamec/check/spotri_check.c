#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int spotri_check(char *uplo, integer *n, float *a, integer *lda, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    /* Function Body */
    *info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*n))
    {
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SPOTRI", &i__1);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    /* Check for singularity */
    i__1 = *n;
    for (*info = 1;
            *info <= i__1;
            ++(*info))
    {
        if (a[*info + *info * a_dim1] == 0.f)
        {
            return LAPACK_FAILURE;
        }
    }
    *info = 0;
    return LAPACK_SUCCESS;
}

