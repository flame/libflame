#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;

int sgebrd_check(integer *m, integer *n, float *a, integer *lda, float *d__, float *e, float *tauq, float *taup, float *work, integer * lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    integer nb, minmn;
    integer lwkopt;
    logical lquery;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "sgebrd inputs: m %d, n %d, lda %d\n", *m, *n, *lda);
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
    /* Computing MAX */
    i__1 = 1;
    i__2 = ilaenv_(&c__1, "SGEBRD", " ", m, n, &c_n1, &c_n1); // , expr subst
    nb = max(i__1,i__2);
    lwkopt = (*m + *n) * nb;
    work[1] = (float) lwkopt;
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
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = max(1,*m);
        if (*lwork < max(i__1,*n) && ! lquery)
        {
            *info = -10;
        }
    }
    if (*info < 0)
    {
        i__1 = -(*info);
        xerbla_("SGEBRD", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    minmn = min(*m,*n);
    if (minmn == 0)
    {
        work[1] = 1.f;
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
