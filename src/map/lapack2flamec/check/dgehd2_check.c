#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int dgehd2_check(int *n, int *ilo, int *ihi, double *a, int *lda, double *tau, double *work, int *info)
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
    if (*n < 0)
    {
        *info = -1;
    }
    else if (*ilo < 1 || *ilo > max(1,*n))
    {
        *info = -2;
    }
    else if (*ihi < min(*ilo,*n) || *ihi > *n)
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
        xerbla_("DGEHD2", &i__1);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}
