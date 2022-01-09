#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int ztrti2_check(char *uplo, char *diag, integer *n, dcomplex *a, integer *lda, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    logical upper;
    logical nounit;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "ztrti2 inputs: uplo %c, diag %c, n %d, lda %d\n", *uplo, *diag, *n, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

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
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZTRTI2", &i__1);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}
