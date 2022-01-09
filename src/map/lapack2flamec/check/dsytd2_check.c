#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int dsytd2_check(char *uplo, integer *n, double *a, integer * lda, double *d__, double *e, double *tau, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    logical upper;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "hetd2-dsytd2 inputs: uplo %c, n %d, lda %d\n", *uplo, *n, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < max(1,*n))
    {
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DSYTD2", &i__1);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if (*n <= 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
