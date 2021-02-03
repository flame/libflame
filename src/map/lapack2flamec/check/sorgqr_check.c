#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static int c__1 = 1;
static int c_n1 = -1;

int sorgqr_check(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;
    /* Local variables */
    int nb;
    int lwkopt;
    logical lquery;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "sorgqr inputs: m %d, n %d, k %d, lda %d\n", *m, *n, *k, *lda);
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
    nb = ilaenv_(&c__1, "SORGQR", " ", m, n, k, &c_n1);
    lwkopt = max(1,*n) * nb;
    work[1] = (float) lwkopt;
    lquery = *lwork == -1;
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
    else if (*lwork < max(1,*n) && ! lquery)
    {
        *info = -8;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORGQR", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if (*n <= 0)
    {
        work[1] = 1.f;
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
