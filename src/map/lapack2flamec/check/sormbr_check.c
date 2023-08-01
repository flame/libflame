#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;

int sormbr_check(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    char ch__1[2];

    /* Local variables */
    integer nb, nq, nw;
    logical left;
    logical notran, applyq;
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
    applyq = lsame_(vect, "Q");
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");
    lquery = *lwork == -1;
    /* NQ is the order of Q or P and NW is the minimum dimension of WORK */
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
    if (! applyq && ! lsame_(vect, "P"))
    {
        *info = -1;
    }
    else if (! left && ! lsame_(side, "R"))
    {
        *info = -2;
    }
    else if (! notran && ! lsame_(trans, "T"))
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
    else if (*k < 0)
    {
        *info = -6;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = fla_min(nq,*k); // , expr subst
        if (applyq && *lda < fla_max(1,nq) || ! applyq && *lda < fla_max(i__1,i__2))
        {
            *info = -8;
        }
        else if (*ldc < fla_max(1,*m))
        {
            *info = -11;
        }
        else if (*lwork < fla_max(1,nw) && ! lquery)
        {
            *info = -13;
        }
    }
    if (*info == 0)
    {
        if (applyq)
        {
            if (left)
            {
                i__1 = *m - 1;
                i__2 = *m - 1;
                nb = ilaenv_(&c__1, "SORMQR", ch__1, &i__1, n, &i__2, &c_n1);
            }
            else
            {
                i__1 = *n - 1;
                i__2 = *n - 1;
                nb = ilaenv_(&c__1, "SORMQR", ch__1, m, &i__1, &i__2, &c_n1);
            }
        }
        else
        {
            if (left)
            {
                i__1 = *m - 1;
                i__2 = *m - 1;
                nb = ilaenv_(&c__1, "SORMLQ", ch__1, &i__1, n, &i__2, &c_n1);
            }
            else
            {
                i__1 = *n - 1;
                i__2 = *n - 1;
                nb = ilaenv_(&c__1, "SORMLQ", ch__1, m, &i__1, &i__2, &c_n1);
            }
        }
        lwkopt = fla_max(1,nw) * nb;
        work[1] = (float) lwkopt;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORMBR", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    work[1] = 1.f;
    if (*m == 0 || *n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
