#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;

int cungbr_check(char *vect, integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    integer mn;
    integer iinfo;
    logical wantq;
    extern integer 
      cunglq_check( integer *, integer *, integer *, scomplex *, 
                    integer *, scomplex *, scomplex *, integer *, integer *), 
      cungqr_check(integer *, integer *, integer *, scomplex *, 
              integer *, scomplex *, scomplex *, integer *, integer *);
    integer lwkopt;
    logical lquery;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "cungbr inputs: vect %c, m %d, n %d, k %d, lda %d\n", *vect, *m, *n, *k, *lda);
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
    wantq = lsame_(vect, "Q");
    mn = min(*m,*n);
    lquery = *lwork == -1;
    if (! wantq && ! lsame_(vect, "P"))
    {
        *info = -1;
    }
    else if (*m < 0)
    {
        *info = -2;
    }
    else if (*n < 0 || wantq && (*n > *m || *n < min(*m,*k)) || ! wantq && ( *m > *n || *m < min(*n,*k)))
    {
        *info = -3;
    }
    else if (*k < 0)
    {
        *info = -4;
    }
    else if (*lda < max(1,*m))
    {
        *info = -6;
    }
    else if (*lwork < max(1,mn) && ! lquery)
    {
        *info = -9;
    }
    if (*info == 0)
    {
        work[1].real = 1.f;
        work[1].imag = 0.f; // , expr subst
        if (wantq)
        {
            if (*m >= *k)
            {
                cungqr_check(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            }
            else
            {
                if (*m > 1)
                {
                    i__1 = *m - 1;
                    i__2 = *m - 1;
                    i__3 = *m - 1;
                    cungqr_check(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, & tau[1], &work[1], &c_n1, &iinfo);
                }
            }
        }
        else
        {
            if (*k < *n)
            {
                cunglq_check(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            }
            else
            {
                if (*n > 1)
                {
                    i__1 = *n - 1;
                    i__2 = *n - 1;
                    i__3 = *n - 1;
                    cunglq_check(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, & tau[1], &work[1], &c_n1, &iinfo);
                }
            }
        }
        lwkopt = work[1].real;
        lwkopt = max(lwkopt,mn);
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CUNGBR", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        work[1].real = (float) lwkopt;
        work[1].imag = 0.f; // , expr subst
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        work[1].real = 1.f;
        work[1].imag = 0.f; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
