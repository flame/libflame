#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" 

int zgehd2_check(integer *n, integer *ilo, integer *ihi, dcomplex *a, integer *lda, dcomplex *tau, dcomplex * work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

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
    else if (*ilo < 1 || *ilo > fla_max(1,*n))
    {
        *info = -2;
    }
    else if (*ihi < fla_min(*ilo,*n) || *ihi > *n)
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
        xerbla_("ZGEHD2", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}
