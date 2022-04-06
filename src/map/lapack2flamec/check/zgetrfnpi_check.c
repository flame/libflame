/*
 * Copyright (c) 2021 Advanced Micro Devices, Inc. All rights reserved.
 * */

#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int zgetrfnpi_check(integer *m, integer *n, integer *nfact, dcomplex *a, integer *lda, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
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
    else if ((*nfact < 0) || (*nfact > min(*m,*n)))
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
        xerbla_("ZGETRFNPI", &i__1);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
     if (*m == 0 || *n == 0 || *nfact == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
