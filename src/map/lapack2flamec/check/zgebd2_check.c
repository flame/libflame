#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int zgebd2_check(int *m, int *n, dcomplex *a, int *lda, doublereal *d__, doublereal *e, dcomplex *tauq, dcomplex *taup, dcomplex *work, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;

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
        xerbla_("ZGEBD2", &i__1);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}
