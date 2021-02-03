#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int dgebd2_check(int *m, int *n, double *a, int * lda, double *d__, double *e, double *tauq, double * taup, double *work, int *info)
{
    /* System generated locals */
  int a_dim1, a_offset, i__1;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "dgebd2 inputs: m %d, n %d, lda %d\n", *m, *n, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    --work;
    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < max(1,*m))
    {
        *info = -4;
    }
    if (*info < 0)
    {
        i__1 = -(*info);
        xerbla_("DGEBD2", &i__1);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}
