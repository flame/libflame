#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static integer c_n1 = -1;

int sorgbr_check(char *vect, integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    integer mn;
    integer iinfo;
    logical wantq;
    extern int
    sorglq_check( integer *, integer *, integer *, float *, integer *,
             float *, float * , integer *, integer *),
                  sorgqr_check(integer *, integer *, integer *, float *, integer *,
                          float *, float *, integer *, integer *);
    integer lwkopt;
    logical lquery;

#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "sorgbr inputs: vect %c, m %d, n %d, k %d, lda %d\n", *vect, *m, *n, *k, *lda);
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
    mn = fla_min(*m,*n);
    lquery = *lwork == -1;
    if (! wantq && ! lsame_(vect, "P"))
    {
        *info = -1;
    }
    else if (*m < 0)
    {
        *info = -2;
    }
    else if (*n < 0 || wantq && (*n > *m || *n < fla_min(*m,*k)) || ! wantq && ( *m > *n || *m < fla_min(*n,*k)))
    {
        *info = -3;
    }
    else if (*k < 0)
    {
        *info = -4;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -6;
    }
    else if (*lwork < fla_max(1,mn) && ! lquery)
    {
        *info = -9;
    }
    if (*info == 0)
    {
        work[1] = 1.f;
        if (wantq)
        {
            if (*m >= *k)
            {
                sorgqr_check(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            }
            else
            {
                if (*m > 1)
                {
                    i__1 = *m - 1;
                    i__2 = *m - 1;
                    i__3 = *m - 1;
                    sorgqr_check(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, & tau[1], &work[1], &c_n1, &iinfo);
                }
            }
        }
        else
        {
            if (*k < *n)
            {
                sorglq_check(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            }
            else
            {
                if (*n > 1)
                {
                    i__1 = *n - 1;
                    i__2 = *n - 1;
                    i__3 = *n - 1;
                    sorglq_check(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, & tau[1], &work[1], &c_n1, &iinfo);
                }
            }
        }
        lwkopt = work[1];
        lwkopt = fla_max(lwkopt,mn);
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORGBR", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        work[1] = (float) lwkopt;
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        work[1] = 1.f;
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
