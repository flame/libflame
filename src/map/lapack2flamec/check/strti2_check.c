#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" /* Table of constant values */

int strti2_check(char *uplo, char *diag, integer *n, float *a, integer *lda, integer *info)
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
    upper = lsame_(uplo, "U", 1, 1);
    nounit = lsame_(diag, "N", 1, 1);
    if (! upper && ! lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if (! nounit && ! lsame_(diag, "U", 1, 1))
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
        xerbla_("STRTI2", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}
