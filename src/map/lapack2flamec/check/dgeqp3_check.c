#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;

int dgeqp3_check(integer *m, integer *n, double *a, integer * lda, integer *jpvt, double *tau, double *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    integer nb, iws;
    integer minmn;
    integer lwkopt;
    logical lquery;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "dgeqp3 inputs: m %d, n %d, lda %d, jpvt %d\n", *m, *n, *lda, *jpvt);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
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
    if (*info == 0)
    {
        minmn = min(*m,*n);
        if (minmn == 0)
        {
            iws = 1;
            lwkopt = 1;
        }
        else
        {
            iws = *n * 3 + 1;
            nb = ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1);
            lwkopt = (*n << 1) + (*n + 1) * nb;
        }
        work[1] = (double) lwkopt;
        if (*lwork < iws && ! lquery)
        {
            *info = -8;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGEQP3", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }

    return LAPACK_SUCCESS;
}
