#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;

int zgelqf_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    integer k, nb;
    integer lwkopt;
    logical lquery;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "ZGELQF", " ", m, n, &c_n1, &c_n1);
    lwkopt = *m * nb;
    work[1].real = (double) lwkopt;
    work[1].imag = 0.; // , expr subst
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
    else if (*lwork < max(1,*m) && ! lquery)
    {
        *info = -7;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGELQF", &i__1);
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
        work[1].real = 1.;
        work[1].imag = 0.; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
