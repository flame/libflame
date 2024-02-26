#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int zpotri_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "zpotri inputs: uplo %c, n %d, lda %d\n", *uplo, *n, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    /* Function Body */
    *info = 0;
    if (! lsame_(uplo, "U", 1, 1) && ! lsame_(uplo, "L", 1, 1))
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
        xerbla_("ZPOTRI", &i__1, (ftnlen)6);
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
        i__2 = *info + *info * a_dim1;
        if (a[i__2].real == 0. && a[i__2].imag == 0.)
        {
            return LAPACK_FAILURE;
        }
    }
    *info = 0;
    return LAPACK_SUCCESS;
}
