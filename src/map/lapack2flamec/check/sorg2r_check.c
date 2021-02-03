#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int sorg2r_check(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "sorg2r inputs: m %d, n %d, k %d, lda %d\n", *m, *n, *k, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

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
    else if (*n < 0 || *n > *m)
    {
        *info = -2;
    }
    else if (*k < 0 || *k > *n)
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
        xerbla_("SORG2R", &i__1);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if (*n <= 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
