#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;

int zgesdd_check(char *jobz, int *m, int *n, dcomplex *a, int *lda, double *s, dcomplex *u, int *ldu, dcomplex *vt, int *ldvt, dcomplex *work, int *lwork, double *rwork, int *iwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    /* Local variables */
    int minmn;
    int wrkbl;
    logical wntqa;
    logical wntqn, wntqo, wntqs;
    int mnthr1, mnthr2;
    int minwrk, maxwrk;
    logical wntqas;

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
    --rwork;
    --iwork;
    /* Function Body */
    *info = 0;
    minmn = min(*m,*n);
    mnthr1 = (int) (minmn * 17. / 9.);
    mnthr2 = (int) (minmn * 5. / 3.);
    wntqa = lsame_(jobz, "A");
    wntqs = lsame_(jobz, "S");
    wntqas = wntqa || wntqs;
    wntqo = lsame_(jobz, "O");
    wntqn = lsame_(jobz, "N");
    minwrk = 1;
    maxwrk = 1;
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
    /* CWorkspace refers to complex workspace, and RWorkspace to */
    /* real workspace. NB refers to the optimal block size for the */
    /* immediately following subroutine, as returned by ILAENV.) */
    if (*info == 0 && *m > 0 && *n > 0)
    {
        if (*m >= *n)
        {
            /* There is no complex work space needed for bidiagonal SVD */
            /* The real work space needed for bidiagonal SVD is BDSPAC */
            /* for computing singular values and singular vectors;
            BDSPAN */
            /* for computing singular values only. */
            /* BDSPAC = 5*N*N + 7*N */
            /* BDSPAN = MAX(7*N+4, 3*N+2+SMLSIZ*(SMLSIZ+8)) */
            if (*m >= mnthr1)
            {
                if (wntqn)
                {
                    /* Path 1 (M much larger than N, JOBZ='N') */
                    maxwrk = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + (*n << 1) * ilaenv_(& c__1, "ZGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    minwrk = *n * 3;
                }
                else if (wntqo)
                {
                    /* Path 2 (M much larger than N, JOBZ='O') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n + *n * ilaenv_(&c__1, "ZUNGQR", " ", m, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*n << 1) + (*n << 1) * ilaenv_(& c__1, "ZGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNMBR", "QLN", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNMBR", "PRC", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = *m * *n + *n * *n + wrkbl;
                    minwrk = (*n << 1) * *n + *n * 3;
                }
                else if (wntqs)
                {
                    /* Path 3 (M much larger than N, JOBZ='S') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n + *n * ilaenv_(&c__1, "ZUNGQR", " ", m, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*n << 1) + (*n << 1) * ilaenv_(& c__1, "ZGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNMBR", "QLN", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNMBR", "PRC", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = *n * *n + wrkbl;
                    minwrk = *n * *n + *n * 3;
                }
                else if (wntqa)
                {
                    /* Path 4 (M much larger than N, JOBZ='A') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n + *m * ilaenv_(&c__1, "ZUNGQR", " ", m, m, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*n << 1) + (*n << 1) * ilaenv_(& c__1, "ZGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNMBR", "QLN", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNMBR", "PRC", n, n, n, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = *n * *n + wrkbl;
                    minwrk = *n * *n + (*n << 1) + *m;
                }
            }
            else if (*m >= mnthr2)
            {
                /* Path 5 (M much larger than N, but not as much as MNTHR1) */
                maxwrk = (*n << 1) + (*m + *n) * ilaenv_(&c__1, "ZGEBRD", " ", m, n, &c_n1, &c_n1);
                minwrk = (*n << 1) + *m;
                if (wntqo)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNGBR", "P", n, n, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNGBR", "Q", m, n, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    maxwrk += *m * *n;
                    minwrk += *n * *n;
                }
                else if (wntqs)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNGBR", "P", n, n, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNGBR", "Q", m, n, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                }
                else if (wntqa)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNGBR", "P", n, n, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *m * ilaenv_(&c__1, "ZUNGBR", "Q", m, m, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                }
            }
            else
            {
                /* Path 6 (M at least N, but not much larger) */
                maxwrk = (*n << 1) + (*m + *n) * ilaenv_(&c__1, "ZGEBRD", " ", m, n, &c_n1, &c_n1);
                minwrk = (*n << 1) + *m;
                if (wntqo)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNMBR", "PRC", n, n, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNMBR", "QLN", m, n, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    maxwrk += *m * *n;
                    minwrk += *n * *n;
                }
                else if (wntqs)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNMBR", "PRC", n, n, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNMBR", "QLN", m, n, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                }
                else if (wntqa)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *n * ilaenv_(&c__1, "ZUNGBR", "PRC", n, n, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*n << 1) + *m * ilaenv_(&c__1, "ZUNGBR", "QLN", m, m, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                }
            }
        }
        else
        {
            /* There is no complex work space needed for bidiagonal SVD */
            /* The real work space needed for bidiagonal SVD is BDSPAC */
            /* for computing singular values and singular vectors;
            BDSPAN */
            /* for computing singular values only. */
            /* BDSPAC = 5*M*M + 7*M */
            /* BDSPAN = MAX(7*M+4, 3*M+2+SMLSIZ*(SMLSIZ+8)) */
            if (*n >= mnthr1)
            {
                if (wntqn)
                {
                    /* Path 1t (N much larger than M, JOBZ='N') */
                    maxwrk = *m + *m * ilaenv_(&c__1, "ZGELQF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + (*m << 1) * ilaenv_(& c__1, "ZGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    minwrk = *m * 3;
                }
                else if (wntqo)
                {
                    /* Path 2t (N much larger than M, JOBZ='O') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "ZGELQF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m + *m * ilaenv_(&c__1, "ZUNGLQ", " ", m, n, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*m << 1) + (*m << 1) * ilaenv_(& c__1, "ZGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNMBR", "PRC", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNMBR", "QLN", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = *m * *n + *m * *m + wrkbl;
                    minwrk = (*m << 1) * *m + *m * 3;
                }
                else if (wntqs)
                {
                    /* Path 3t (N much larger than M, JOBZ='S') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "ZGELQF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m + *m * ilaenv_(&c__1, "ZUNGLQ", " ", m, n, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*m << 1) + (*m << 1) * ilaenv_(& c__1, "ZGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNMBR", "PRC", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNMBR", "QLN", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = *m * *m + wrkbl;
                    minwrk = *m * *m + *m * 3;
                }
                else if (wntqa)
                {
                    /* Path 4t (N much larger than M, JOBZ='A') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "ZGELQF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m + *n * ilaenv_(&c__1, "ZUNGLQ", " ", n, n, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*m << 1) + (*m << 1) * ilaenv_(& c__1, "ZGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNMBR", "PRC", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNMBR", "QLN", m, m, m, &c_n1); // , expr subst
                    wrkbl = max(i__1,i__2);
                    maxwrk = *m * *m + wrkbl;
                    minwrk = *m * *m + (*m << 1) + *n;
                }
            }
            else if (*n >= mnthr2)
            {
                /* Path 5t (N much larger than M, but not as much as MNTHR1) */
                maxwrk = (*m << 1) + (*m + *n) * ilaenv_(&c__1, "ZGEBRD", " ", m, n, &c_n1, &c_n1);
                minwrk = (*m << 1) + *n;
                if (wntqo)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNGBR", "P", m, n, m, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNGBR", "Q", m, m, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    maxwrk += *m * *n;
                    minwrk += *m * *m;
                }
                else if (wntqs)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNGBR", "P", m, n, m, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNGBR", "Q", m, m, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                }
                else if (wntqa)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *n * ilaenv_(&c__1, "ZUNGBR", "P", n, n, m, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNGBR", "Q", m, m, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                }
            }
            else
            {
                /* Path 6t (N greater than M, but not much larger) */
                maxwrk = (*m << 1) + (*m + *n) * ilaenv_(&c__1, "ZGEBRD", " ", m, n, &c_n1, &c_n1);
                minwrk = (*m << 1) + *n;
                if (wntqo)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNMBR", "PRC", m, n, m, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNMBR", "QLN", m, m, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    maxwrk += *m * *n;
                    minwrk += *m * *m;
                }
                else if (wntqs)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNGBR", "PRC", m, n, m, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNGBR", "QLN", m, m, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                }
                else if (wntqa)
                {
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *n * ilaenv_(&c__1, "ZUNGBR", "PRC", n, n, m, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = (*m << 1) + *m * ilaenv_(&c__1, "ZUNGBR", "QLN", m, m, n, &c_n1); // , expr subst
                    maxwrk = max(i__1,i__2);
                }
            }
        }
        maxwrk = max(maxwrk,minwrk);
    }
    if (*info == 0)
    {
        work[1].real = (double) maxwrk;
        work[1].imag = 0.; // , expr subst
        if (*lwork < minwrk && *lwork != -1)
        {
            *info = -13;
        }
    }
    /* Quick returns */
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGESDD", &i__1);
        return LAPACK_FAILURE;
    }
    if (*lwork == -1)
    {
        return LAPACK_QUERY_RETURN;
    }
    if (*m == 0 || *n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}

