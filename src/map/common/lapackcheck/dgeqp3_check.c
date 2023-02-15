#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" /* Table of constant values */
static int c__1 = 1;
static int c_n1 = -1;

int dgeqp3_check(int *m, int *n, double *a, int * lda, int *jpvt, double *tau, double *work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;
    /* Local variables */
    int nb, iws;
    int minmn;
    int lwkopt;
    logical lquery;

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
    /* Quick return if possible. */
    if (minmn == 0)
    {
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
