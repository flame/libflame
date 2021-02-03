#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int cunm2r_check(char *side, char *trans, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, c_dim1, c_offset, i__1;
    /* Builtin functions */
    /* Local variables */
    int nq;
    logical left;
    logical notran;
    
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "cunm2r inputs: side %c, trans %c, m %d, n %d, k %d, lda %d, ldc %d\n", *side, *trans, *m, *n, *k, *lda, *ldc);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    *info = 0;
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");
    /* NQ is the order of Q */
    if (left)
    {
        nq = *m;
    }
    else
    {
        nq = *n;
    }
    if (! left && ! lsame_(side, "R"))
    {
        *info = -1;
    }
    else if (! notran && ! lsame_(trans, "C"))
    {
        *info = -2;
    }
    else if (*m < 0)
    {
        *info = -3;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*k < 0 || *k > nq)
    {
        *info = -5;
    }
    else if (*lda < max(1,nq))
    {
        *info = -7;
    }
    else if (*ldc < max(1,*m))
    {
        *info = -10;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CUNM2R", &i__1);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0 || *k == 0)
    {
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
