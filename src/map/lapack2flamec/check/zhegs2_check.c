#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int zhegs2_check(integer *itype, char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
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
    upper = lsame_(uplo, "U", 1, 1);
    if (*itype < 1 || *itype > 3)
    {
        *info = -1;
    }
    else if (! upper && ! lsame_(uplo, "L", 1, 1))
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
    else if (*ldb < fla_max(1,*n))
    {
        *info = -7;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHEGS2", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}

