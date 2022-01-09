#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int chegs2_check(integer *itype, char *uplo, integer *n, scomplex * a, integer *lda, scomplex *b, integer *ldb, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    logical upper;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "chegs2 inputs: itype %d, uplo %c, n %d, lda %d, ldb %d\n", *itype, *uplo, *n, *lda, *ldb);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

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
        xerbla_("CHEGS2", &i__1);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}
    
