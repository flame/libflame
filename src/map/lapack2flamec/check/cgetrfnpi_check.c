/*
 * Copyright (c) 2020 Advanced Micro Devices, Inc. All rights reserved.
 * */

#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int cgetrfnpi_check(int *m, int *n, int *nfact, scomplex *a, int *lda, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
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
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGETRFNPI", &i__1);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if(*nfact < 0 || *nfact > min(*m,*n))
    {
        return LAPACK_FAILURE;
    }
    if (*m == 0 || *n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
