#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;

int zunmqr_check(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    char ch__1[2];

    /* Local variables */
    integer nb, nq, nw;
    logical left;
    logical notran;
    integer lwkopt;
    logical lquery;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "zunmqr inputs: side %c, trans %c, m %d, n %d, k %d, lda %d, ldc %d\n", *side, *trans, *m, *n, *k, *lda, *ldc);
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
    lquery = *lwork == -1;
    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left)
    {
        nq = *m;
        nw = *n;
    }
    else
    {
        nq = *n;
        nw = *m;
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
    else if (*lwork < max(1,nw) && ! lquery)
    {
        *info = -12;
    }
    if (*info == 0)
    {
        /* Determine the block size. NB may be at most NBMAX, where NBMAX */
        /* is used to define the local array T. */
        /* Computing MIN */
        i__1 = 64;
        i__2 = ilaenv_(&c__1, "ZUNMQR", ch__1, m, n, k, &c_n1); // , expr subst
        nb = min(i__1,i__2);
        lwkopt = max(1,nw) * nb;
        work[1].real = (double) lwkopt;
        work[1].imag = 0.; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZUNMQR", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0 || *k == 0)
    {
        work[1].real = 1.;
        work[1].imag = 0.; // , expr subst
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
