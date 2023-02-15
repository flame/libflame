#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int zungl2_check(int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < *m)
    {
        *info = -2;
    }
    else if (*k < 0 || *k > *m)
    {
        *info = -3;
    }
    else if (*lda < max(1,*m))
    {
        *info = -5;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZUNGL2", &i__1);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if (*m <= 0)
    {
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
