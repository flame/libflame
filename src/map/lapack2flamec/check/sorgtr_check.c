#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;

int sorgtr_check(char *uplo, integer *n, float *a, integer *lda, float *tau, float *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    integer nb;
    logical upper;
    logical lquery;
    integer lwkopt;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    upper = lsame_(uplo, "U");
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
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n - 1; // , expr subst
        if (*lwork < fla_max(i__1,i__2) && ! lquery)
        {
            *info = -7;
        }
    }
    if (*info == 0)
    {
        if (upper)
        {
            i__1 = *n - 1;
            i__2 = *n - 1;
            i__3 = *n - 1;
            nb = ilaenv_(&c__1, "SORGQL", " ", &i__1, &i__2, &i__3, &c_n1);
        }
        else
        {
            i__1 = *n - 1;
            i__2 = *n - 1;
            i__3 = *n - 1;
            nb = ilaenv_(&c__1, "SORGQR", " ", &i__1, &i__2, &i__3, &c_n1);
        }
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n - 1; // , expr subst
        lwkopt = fla_max(i__1,i__2) * nb;
        work[1] = (float) lwkopt;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORGTR", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        work[1] = 1.f;
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
