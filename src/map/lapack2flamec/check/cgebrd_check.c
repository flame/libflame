#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

static int c__1 = 1;
static int c_n1 = -1;

int cgebrd_check(int *m, int *n, scomplex *a, int *lda, real *d__, real *e, scomplex *tauq, scomplex *taup, scomplex *work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    int nb;
    int minmn;
    int lwkopt;
    logical lquery;

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
    i__2 = ilaenv_(&c__1, "CGEBRD", " ", m, n, &c_n1, &c_n1); // , expr subst
    nb = max(i__1,i__2);
    lwkopt = (*m + *n) * nb;
    r__1 = (real) lwkopt;
    work[1].real = r__1;
    work[1].imag = 0.f; // , expr subst
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
        xerbla_("CGEBRD", &i__1);
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
        work[1].real = 1.f;
        work[1].imag = 0.f; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
