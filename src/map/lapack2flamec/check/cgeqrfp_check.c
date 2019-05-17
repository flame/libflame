#include "FLA_lapack2flame_return_defs.h"


#include "FLA_f2c.h" /* Table of constant values */
static int c__1 = 1;
static int c_n1 = -1;
int cgeqrfp_check(int *m, int *n, scomplex *a, int * lda, scomplex *tau, scomplex *work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;
    /* Local variables */
    int k, nb;
    int lwkopt;
    logical lquery;
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "CGEQRF", " ", m, n, &c_n1, &c_n1);
    lwkopt = *n * nb;
    work[1].real = (float) lwkopt;
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
    else if (*lwork < max(1,*n) && ! lquery)
    {
        *info = -7;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGEQRFP", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    k = min(*m,*n);
    if (k == 0)
    {
        work[1].real = 1.f;
        work[1].imag = 0.f; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
