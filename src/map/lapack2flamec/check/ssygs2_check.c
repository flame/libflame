#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int ssygs2_check(int *itype, char *uplo, int *n, float *a, int *lda, float *b, int *ldb, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1;
    logical upper;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (*itype < 1 || *itype > 3)
    {
        *info = -1;
    }
    else if (! upper && ! lsame_(uplo, "L"))
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
    else if (*ldb < max(1,*n))
    {
        *info = -7;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSYGS2", &i__1);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}
