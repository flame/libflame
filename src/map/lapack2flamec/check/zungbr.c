#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static int c_n1 = -1;

int zungbr_check(char *vect, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    int mn;
    int iinfo;
    logical wantq;
    int lwkopt;
    logical lquery;
    extern int 
      zunglq_check(int *, int *, int *, dcomplex *, int *, 
              dcomplex *, dcomplex *, int *, int *), 
      zungqr_check(int *, int *, int *, dcomplex *, int *, 
              dcomplex *, dcomplex *, int *, int *);

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
        work[1].real = 1.;
        work[1].imag = 0.; // , expr subst
        if (wantq)
        {
            if (*m >= *k)
            {
                zungqr_check(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            }
            else
            {
                if (*m > 1)
                {
                    i__1 = *m - 1;
                    i__2 = *m - 1;
                    i__3 = *m - 1;
                    zungqr_check(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, & tau[1], &work[1], &c_n1, &iinfo);
                }
            }
        }
        else
        {
            if (*k < *n)
            {
                zunglq_check(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            }
            else
            {
                if (*n > 1)
                {
                    i__1 = *n - 1;
                    i__2 = *n - 1;
                    i__3 = *n - 1;
                    zunglq_check(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, & tau[1], &work[1], &c_n1, &iinfo);
                }
            }
        }
        lwkopt = (int) work[1].real;
        lwkopt = max(lwkopt,mn);
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZUNGBR", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        work[1].real = (double) lwkopt;
        work[1].imag = 0.; // , expr subst
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        work[1].real = 1.;
        work[1].imag = 0.; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}

