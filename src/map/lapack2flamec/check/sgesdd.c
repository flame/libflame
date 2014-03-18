#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static int c__1 = 1;
static int c_n1 = -1;

int sgesdd_check(char *jobz, int *m, int *n, real *a, int *lda, real *s, real *u, int *ldu, real *vt, int *ldvt, real *work, int *lwork, int *iwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    /* Local variables */
    int minmn, wrkbl, mnthr;
    logical wntqa;
    logical wntqn, wntqo, wntqs;
    int bdspac;
    int minwrk,  maxwrk;
    logical wntqas;
    logical lquery;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    minmn = min(*m,*n);
    wntqa = lsame_(jobz, "A");
    wntqs = lsame_(jobz, "S");
    wntqas = wntqa || wntqs;
    wntqo = lsame_(jobz, "O");
    wntqn = lsame_(jobz, "N");
    lquery = *lwork == -1;
    if (! (wntqa || wntqs || wntqo || wntqn))
    {
        *info = -1;
    }
    else if (*m < 0)
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*lda < max(1,*m))
    {
        *info = -5;
    }
    else if (*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < * m)
    {
        *info = -8;
    }
    else if (*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn || wntqo && *m >= *n && *ldvt < *n)
    {
        *info = -10;
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* NB refers to the optimal block size for the immediately */
    /* following subroutine, as returned by ILAENV.) */
    if (*info == 0)
    {
        minwrk = 1;
        maxwrk = 1;
        if (*m >= *n && minmn > 0)
        {
            /* Compute space needed for SBDSDC */
            mnthr = (int) (minmn * 11.f / 6.f);
            if (wntqn)
            {
                bdspac = *n * 7;
            }
            else
            {
                bdspac = *n * 3 * *n + (*n << 2);
            }
            if (*m >= mnthr)
            {
                if (wntqn)
                {
                    /* Path 1 (M much larger than N, JOBZ='N') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, "SGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n; // , expr subst
                    maxwrk = max(i__1,i__2);
                    minwrk = bdspac + *n;
                }
                else if (wntqo)
                {
                    /* Path 2 (M much larger than N, JOBZ='O') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n + *n * ilaenv_(&c__1, "SORGQR", " ", m, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, "SGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "QLN", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = wrkbl + (*n << 1) * *n;
                    minwrk = bdspac + (*n << 1) * *n + *n * 3;
                }
                else if (wntqs)
                {
                    /* Path 3 (M much larger than N, JOBZ='S') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n + *n * ilaenv_(&c__1, "SORGQR", " ", m, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, "SGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "QLN", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = wrkbl + *n * *n;
                    minwrk = bdspac + *n * *n + *n * 3;
                }
                else if (wntqa)
                {
                    /* Path 4 (M much larger than N, JOBZ='A') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n + *m * ilaenv_(&c__1, "SORGQR", " ", m, m, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, "SGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "QLN", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = wrkbl + *n * *n;
                    minwrk = bdspac + *n * *n + (*n << 1) + *m;
                }
            }
            else
            {
                /* Path 5 (M at least N, but not much larger) */
                wrkbl = *n * 3 + (*m + *n) * ilaenv_(&c__1, "SGEBRD", " ", m, n, &c_n1, &c_n1);
                if (wntqn)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    maxwrk = max(i__1,i__2);
                    minwrk = *n * 3 + max(*m,bdspac);
                }
                else if (wntqo)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "QLN", m, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = wrkbl + *m * *n;
                    /* Computing MAX */
                    i__1 = *m;
                    i__2 = *n * *n + bdspac; // , expr subst
                    minwrk = *n * 3 + max(i__1,i__2);
                }
                else if (wntqs)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "QLN", m, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    maxwrk = max(i__1,i__2);
                    minwrk = *n * 3 + max(*m,bdspac);
                }
                else if (wntqa)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *m * ilaenv_(&c__1, "SORMBR" , "QLN", m, m, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR" , "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = bdspac + *n * 3; // , expr subst
                    maxwrk = max(i__1,i__2);
                    minwrk = *n * 3 + max(*m,bdspac);
                }
            }
        }
        else if (minmn > 0)
        {
            /* Compute space needed for SBDSDC */
            mnthr = (int) (minmn * 11.f / 6.f);
            if (wntqn)
            {
                bdspac = *m * 7;
            }
            else
            {
                bdspac = *m * 3 * *m + (*m << 2);
            }
            if (*n >= mnthr)
            {
                if (wntqn)
                {
                    /* Path 1t (N much larger than M, JOBZ='N') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, "SGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m; // , expr subst
                    maxwrk = max(i__1,i__2);
                    minwrk = bdspac + *m;
                }
                else if (wntqo)
                {
                    /* Path 2t (N much larger than M, JOBZ='O') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m + *m * ilaenv_(&c__1, "SORGLQ", " ", m, n, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, "SGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "QLN", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "PRT", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = wrkbl + (*m << 1) * *m;
                    minwrk = bdspac + (*m << 1) * *m + *m * 3;
                }
                else if (wntqs)
                {
                    /* Path 3t (N much larger than M, JOBZ='S') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m + *m * ilaenv_(&c__1, "SORGLQ", " ", m, n, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, "SGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "QLN", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "PRT", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = wrkbl + *m * *m;
                    minwrk = bdspac + *m * *m + *m * 3;
                }
                else if (wntqa)
                {
                    /* Path 4t (N much larger than M, JOBZ='A') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m + *n * ilaenv_(&c__1, "SORGLQ", " ", n, n, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, "SGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "QLN", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "PRT", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = wrkbl + *m * *m;
                    minwrk = bdspac + *m * *m + *m * 3;
                }
            }
            else
            {
                /* Path 5t (N greater than M, but not much larger) */
                wrkbl = *m * 3 + (*m + *n) * ilaenv_(&c__1, "SGEBRD", " ", m, n, &c_n1, &c_n1);
                if (wntqn)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    maxwrk = max(i__1,i__2);
                    minwrk = *m * 3 + max(*n,bdspac);
                }
                else if (wntqo)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "QLN", m, m, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "PRT", m, n, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = wrkbl + *m * *n;
                    /* Computing MAX */
                    i__1 = *n;
                    i__2 = *m * *m + bdspac; // , expr subst
                    minwrk = *m * 3 + max(i__1,i__2);
                }
                else if (wntqs)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "QLN", m, m, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "PRT", m, n, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    maxwrk = max(i__1,i__2);
                    minwrk = *m * 3 + max(*n,bdspac);
                }
                else if (wntqa)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "QLN", m, m, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR" , "PRT", n, n, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    maxwrk = max(i__1,i__2);
                    minwrk = *m * 3 + max(*n,bdspac);
                }
            }
        }
        maxwrk = max(maxwrk,minwrk);
        work[1] = (real) maxwrk;
        if (*lwork < minwrk && ! lquery)
        {
            *info = -12;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGESDD", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
