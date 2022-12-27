#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" 
static integer c__1 = 1;
static integer c_n1 = -1;

int dsytrd_check(char *uplo, integer *n, double *a, integer * lda, double *d__, double *e, double *tau, double * work, integer *lwork, integer *info)
{
    /* System generated locals */
  integer a_dim1, a_offset, i__1;
    /* Local variables */
    integer nb;
    logical upper;
    integer lwkopt;
    logical lquery;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    lquery = *lwork == -1;
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*n))
    {
        *info = -4;
    }
    else if (*lwork < 1 && ! lquery)
    {
        *info = -9;
    }
    if (*info == 0)
    {
        /* Determine the block size. */
        nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1);
        lwkopt = *n * nb;
        work[1] = (double) lwkopt;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DSYTRD", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        work[1] = 1.;
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
