#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;

int zunglq_check(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex * work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    integer nb;
    logical lquery;
    integer lwkopt;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "zunglq inputs: m %d, n %d, k %d, lda %d\n", *m, *n, *k, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "ZUNGLQ", " ", m, n, k, &c_n1);
    lwkopt = max(1,*m) * nb;
    work[1].real = (double) lwkopt;
    work[1].imag = 0.; // , expr subst
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < *m)
    {
        *info = -2;
    }
    else if (*k < 0 || *k > *m)
    {
        *info = -3;
    }
    else if (*lda < max(1,*m))
    {
        *info = -5;
    }
    else if (*lwork < max(1,*m) && ! lquery)
    {
        *info = -8;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZUNGLQ", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if (*m <= 0)
    {
        work[1].real = 1.;
        work[1].imag = 0.; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
