#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static int c__6 = 6;
static int c_n1 = -1;
static int c__9 = 9;
static int c__0 = 0;
static int c__1 = 1;

int dgelsd_check(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double * s, double *rcond, int *rank, double *work, int *lwork, int *iwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    /* Builtin functions */
    double log(double);
    /* Local variables */
    int mm;
    int nlvl;
    int minmn, maxmn, mnthr;
    int wlalsd;
    int liwork, minwrk, maxwrk;
    logical lquery;
    int smlsiz;
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --s;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    minmn = min(*m,*n);
    maxmn = max(*m,*n);
    mnthr = ilaenv_(&c__6, "DGELSD", " ", m, n, nrhs, &c_n1);
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*nrhs < 0)
    {
        *info = -3;
    }
    else if (*lda < max(1,*m))
    {
        *info = -5;
    }
    else if (*ldb < max(1,maxmn))
    {
        *info = -7;
    }
    smlsiz = ilaenv_(&c__9, "DGELSD", " ", &c__0, &c__0, &c__0, &c__0);
    /* Compute workspace. */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* NB refers to the optimal block size for the immediately */
    /* following subroutine, as returned by ILAENV.) */
    minwrk = 1;
    liwork = 1;
    minmn = max(1,minmn);
    /* Computing MAX */
    i__1 = (int) (log((double) minmn / (double) (smlsiz + 1)) / log(2.)) + 1;
    nlvl = max(i__1,0);
    if (*info == 0)
    {
        maxwrk = 0;
        liwork = minmn * 3 * nlvl + minmn * 11;
        mm = *m;
        if (*m >= *n && *m >= mnthr)
        {
            /* Path 1a - overdetermined, with many more rows than columns. */
            mm = *n;
            /* Computing MAX */
            i__1 = maxwrk;
            i__2 = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1); // , expr subst
            maxwrk = max(i__1,i__2);
            /* Computing MAX */
            i__1 = maxwrk;
            i__2 = *n + *nrhs * ilaenv_(&c__1, "DORMQR", "LT", m, nrhs, n, &c_n1); // , expr subst
            maxwrk = max(i__1,i__2);
        }
        if (*m >= *n)
        {
            /* Path 1 - overdetermined or exactly determined. */
            /* Computing MAX */
            i__1 = maxwrk;
            i__2 = *n * 3 + (mm + *n) * ilaenv_(&c__1, "DGEBRD" , " ", &mm, n, &c_n1, &c_n1); // , expr subst
            maxwrk = max(i__1,i__2);
            /* Computing MAX */
            i__1 = maxwrk;
            i__2 = *n * 3 + *nrhs * ilaenv_(&c__1, "DORMBR", "QLT", &mm, nrhs, n, &c_n1); // , expr subst
            maxwrk = max(i__1,i__2);
            /* Computing MAX */
            i__1 = maxwrk;
            i__2 = *n * 3 + (*n - 1) * ilaenv_(&c__1, "DORMBR", "PLN", n, nrhs, n, &c_n1); // , expr subst
            maxwrk = max(i__1,i__2);
            /* Computing 2nd power */
            i__1 = smlsiz + 1;
            wlalsd = *n * 9 + (*n << 1) * smlsiz + (*n << 3) * nlvl + *n * * nrhs + i__1 * i__1;
            /* Computing MAX */
            i__1 = maxwrk;
            i__2 = *n * 3 + wlalsd; // , expr subst
            maxwrk = max(i__1,i__2);
            /* Computing MAX */
            i__1 = *n * 3 + mm, i__2 = *n * 3 + *nrhs;
            i__1 = max(i__1,i__2);
            i__2 = *n * 3 + wlalsd; // ; expr subst
            minwrk = max(i__1,i__2);
        }
        if (*n > *m)
        {
            /* Computing 2nd power */
            i__1 = smlsiz + 1;
            wlalsd = *m * 9 + (*m << 1) * smlsiz + (*m << 3) * nlvl + *m * * nrhs + i__1 * i__1;
            if (*n >= mnthr)
            {
                /* Path 2a - underdetermined, with many more columns */
                /* than rows. */
                maxwrk = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *m * *m + (*m << 2) + (*m << 1) * ilaenv_(&c__1, "DGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                maxwrk = max(i__1,i__2);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *m * *m + (*m << 2) + *nrhs * ilaenv_(& c__1, "DORMBR", "QLT", m, nrhs, m, &c_n1); // , expr subst
                maxwrk = max(i__1,i__2);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *m * *m + (*m << 2) + (*m - 1) * ilaenv_(&c__1, "DORMBR", "PLN", m, nrhs, m, &c_n1); // , expr subst
                maxwrk = max(i__1,i__2);
                if (*nrhs > 1)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = *m * *m + *m + *m * *nrhs; // , expr subst
                    maxwrk = max(i__1,i__2);
                }
                else
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = *m * *m + (*m << 1); // , expr subst
                    maxwrk = max(i__1,i__2);
                }
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *m + *nrhs * ilaenv_(&c__1, "DORMLQ", "LT", n, nrhs, m, &c_n1); // , expr subst
                maxwrk = max(i__1,i__2);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *m * *m + (*m << 2) + wlalsd; // , expr subst
                maxwrk = max(i__1,i__2);
                /* XXX: Ensure the Path 2a case below is triggered. The workspace */
                /* calculation should use queries for all routines eventually. */
                /* Computing MAX */
                /* Computing MAX */
                i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4);
                i__3 = max(i__3,*nrhs);
                i__4 = *n - *m * 3; // ; expr subst
                i__1 = maxwrk;
                i__2 = (*m << 2) + *m * *m + max(i__3,i__4); // , expr subst
                maxwrk = max(i__1,i__2);
            }
            else
            {
                /* Path 2 - remaining underdetermined cases. */
                maxwrk = *m * 3 + (*n + *m) * ilaenv_(&c__1, "DGEBRD", " ", m, n, &c_n1, &c_n1);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *m * 3 + *nrhs * ilaenv_(&c__1, "DORMBR" , "QLT", m, nrhs, n, &c_n1); // , expr subst
                maxwrk = max(i__1,i__2);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR", "PLN", n, nrhs, m, &c_n1); // , expr subst
                maxwrk = max(i__1,i__2);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *m * 3 + wlalsd; // , expr subst
                maxwrk = max(i__1,i__2);
            }
            /* Computing MAX */
            i__1 = *m * 3 + *nrhs, i__2 = *m * 3 + *m;
            i__1 = max(i__1,i__2);
            i__2 = *m * 3 + wlalsd; // ; expr subst
            minwrk = max(i__1,i__2);
        }
        minwrk = min(minwrk,maxwrk);
        work[1] = (double) maxwrk;
        iwork[1] = liwork;
        if (*lwork < minwrk && ! lquery)
        {
            *info = -12;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGELSD", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        work[1] = (double) maxwrk;
        iwork[1] = liwork;
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible. */
    if (*m == 0 || *n == 0)
    {
        *rank = 0;
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
