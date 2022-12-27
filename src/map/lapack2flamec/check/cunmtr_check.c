#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" 
static integer c__1 = 1;
static integer c_n1 = -1;

int cunmtr_check(char *side, char *uplo, char *trans, integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *c__, integer *ldc, scomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__2, i__3;
    char ch__1[2];
    /* Builtin functions */
    /* Subroutine */

    /* Local variables */
    integer nb, nq, nw;
    logical left;
    logical upper;
    integer lwkopt;
    logical lquery;

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
    upper = lsame_(uplo, "U");
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
    else if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -2;
    }
    else if (! lsame_(trans, "N") && ! lsame_(trans, "C"))
    {
        *info = -3;
    }
    else if (*m < 0)
    {
        *info = -4;
    }
    else if (*n < 0)
    {
        *info = -5;
    }
    else if (*lda < fla_max(1,nq))
    {
        *info = -7;
    }
    else if (*ldc < fla_max(1,*m))
    {
        *info = -10;
    }
    else if (*lwork < fla_max(1,nw) && ! lquery)
    {
        *info = -12;
    }
    if (*info == 0)
    {
        if (upper)
        {
            if (left)
            {
                i__2 = *m - 1;
                i__3 = *m - 1;
                nb = ilaenv_(&c__1, "CUNMQL", ch__1, &i__2, n, &i__3, &c_n1);
            }
            else
            {
                i__2 = *n - 1;
                i__3 = *n - 1;
                nb = ilaenv_(&c__1, "CUNMQL", ch__1, m, &i__2, &i__3, &c_n1);
            }
        }
        else
        {
            if (left)
            {
                i__2 = *m - 1;
                i__3 = *m - 1;
                nb = ilaenv_(&c__1, "CUNMQR", ch__1, &i__2, n, &i__3, &c_n1);
            }
            else
            {
                i__2 = *n - 1;
                i__3 = *n - 1;
                nb = ilaenv_(&c__1, "CUNMQR", ch__1, m, &i__2, &i__3, &c_n1);
            }
        }
        lwkopt = fla_max(1,nw) * nb;
        work[1].real = (float) lwkopt;
        work[1].imag = 0.f; // , expr subst
    }
    if (*info != 0)
    {
        i__2 = -(*info);
        xerbla_("CUNMTR", &i__2);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0 || nq == 1)
    {
        work[1].real = 1.f;
        work[1].imag = 0.f; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}

