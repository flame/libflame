#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" /* Table of constant values */

int zgeqpf_check(integer *m, integer *n, dcomplex *a, integer *lda, integer *jpvt, dcomplex *tau, dcomplex *work, double *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "zgeqpf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "\n", *m, *n, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --tau;
    --work;
    --rwork;
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
    else if (*lda < fla_max(1,*m))
    {
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGEQPF", &i__1);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
      {
        return LAPACK_QUICK_RETURN;
      }
    return LAPACK_SUCCESS;
}
