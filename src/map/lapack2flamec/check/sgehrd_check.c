#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" 
static int c__1 = 1;
static int c_n1 = -1;

int sgehrd_check(int *n, int *ilo, int *ihi, float *a, int *lda, float *tau, float *work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    int i__;
    int nb, nh;
    int lwkopt;
    logical lquery;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "sgehrd inputs: n %d, ilo %d, ihi %d, lda %d\n", *n, *ilo, *ihi, *lda);
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
    /* Computing MIN */
    i__1 = 64;
    i__2 = ilaenv_(&c__1, "SGEHRD", " ", n, ilo, ihi, &c_n1); // , expr subst
    nb = min(i__1,i__2);
    lwkopt = *n * nb;
    work[1] = (float) lwkopt;
    lquery = *lwork == -1;
    if (*n < 0)
    {
        *info = -1;
    }
    else if (*ilo < 1 || *ilo > max(1,*n))
    {
        *info = -2;
    }
    else if (*ihi < min(*ilo,*n) || *ihi > *n)
    {
        *info = -3;
    }
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    else if (*lwork < max(1,*n) && ! lquery)
    {
        *info = -8;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGEHRD", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Set elements 1:ILO-1 and IHI:N-1 of TAU to zero */
    i__1 = *ilo - 1;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        tau[i__] = 0.f;
        /* L10: */
    }
    i__1 = *n - 1;
    for (i__ = max(1,*ihi);
            i__ <= i__1;
            ++i__)
    {
        tau[i__] = 0.f;
        /* L20: */
    }
    /* Quick return if possible */
    nh = *ihi - *ilo + 1;
    if (nh <= 1)
    {
        work[1] = 1.f;
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
