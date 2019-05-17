#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static int c_n1 = -1;

int sorgbr_check(char *vect, int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    int mn;
    int iinfo;
    logical wantq;
    extern int
    sorglq_check( int *, int *, int *, float *, int *,
             float *, float * , int *, int *),
                  sorgqr_check(int *, int *, int *, float *, int *,
                          float *, float *, int *, int *);
    int lwkopt;
    logical lquery;

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
        lwkopt = max(lwkopt,mn);
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
