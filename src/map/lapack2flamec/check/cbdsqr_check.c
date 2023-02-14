#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int cbdsqr_check(char *uplo, integer *n, integer *ncvt, integer * nru, integer *ncc, float *d__, float *e, scomplex *vt, integer *ldvt, scomplex *u, integer *ldu, scomplex *c__, integer *ldc, float *rwork, integer *info)
{
    /* System generated locals */
  integer c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1;
    logical lower;

#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "cbdsqr inputs: uplo %c, n %d, ncvt %d, nru %d, ncc %d, ldvt %d, ldu %d, ldc %d\n", *uplo, *n, *ncvt, *nru, *ncc, *ldvt, *ldu, *ldc);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* Parameter adjustments */
    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --rwork;
    /* Function Body */
    *info = 0;
    lower = lsame_(uplo, "L");
    if (! lsame_(uplo, "U") && ! lower)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*ncvt < 0)
    {
        *info = -3;
    }
    else if (*nru < 0)
    {
        *info = -4;
    }
    else if (*ncc < 0)
    {
        *info = -5;
    }
    else if (*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < fla_max(1,*n))
    {
        *info = -9;
    }
    else if (*ldu < fla_max(1,*nru))
    {
        *info = -11;
    }
    else if (*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < fla_max(1,*n))
    {
        *info = -13;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CBDSQR", &i__1);
        return LAPACK_FAILURE;
    }
    if (*n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
