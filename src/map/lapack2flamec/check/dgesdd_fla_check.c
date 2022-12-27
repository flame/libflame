#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static integer c__6 = 6;
static integer c__0 = 0;
static integer c_n1 = -1;

int dgesdd_fla_check(char *jobu, char *jobvt, integer *m, integer *n, double * a, integer *lda, double *s, double *u, integer *ldu, double *vt, integer *ldvt, double *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__2, i__3;
    char ch__1[2];
    /* Local variables */
    double dum[1];
    integer ierr, lwork_dgebrd__, lwork_dgelqf__, lwork_dgeqrf__;
    integer minmn, wrkbl, mnthr;
    logical wntua, wntva, wntun, wntuo, wntvn, wntvo, wntus, wntvs;
    integer bdspac;
    extern int 
      dgebrd_(integer *, integer *, double *, integer *, double *, double *, double *, double *, double *, integer *, integer *),
      dgelqf_(integer *, integer *, double *, integer *, double *, double *, integer *, integer *), 
      dgeqrf_(integer *, integer *, double *, integer *, double *, double *, integer *, integer *), 
      dorgbr_(char *, integer *, integer *, integer *, double *, integer *, double *, double *, integer *, integer *),
      dorglq_(integer *, integer *, integer *, double *, integer *, double *, double *, integer *, integer *), 
      dorgqr_(integer *, integer *, integer *, double *, integer *, double *, double *, integer *, integer *);
    integer minwrk, maxwrk;
    logical lquery, wntuas, wntvas;
    integer lwork_dorgbr_p__, lwork_dorgbr_q__, lwork_dorglq_m__, lwork_dorglq_n__, lwork_dorgqr_m__, lwork_dorgqr_n__;

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
    /* Function Body */
    *info = 0;
    minmn = fla_min(*m,*n);
    wntua = lsame_(jobu, "A");
    wntus = lsame_(jobu, "S");
    wntuas = wntua || wntus;
    wntuo = lsame_(jobu, "O");
    wntun = lsame_(jobu, "N");
    wntva = lsame_(jobvt, "A");
    wntvs = lsame_(jobvt, "S");
    wntvas = wntva || wntvs;
    wntvo = lsame_(jobvt, "O");
    wntvn = lsame_(jobvt, "N");
    lquery = *lwork == -1;
    if ((! (wntua || wntus || wntuo || wntun)) || (! (wntva || wntvs || wntvo || wntvn) || wntvo && wntuo))
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
    else if (*lda < fla_max(1,*m))
    {
        *info = -5;
    }
    else if (*ldu < 1 || wntuas && *ldu < *m || wntuo && *m < *n && *ldu < * m)
    {
        *info = -8;
    }
    else if (*ldvt < 1 || wntua && *ldvt < *n || wntus && *ldvt < minmn || wntuo && *m >= *n && *ldvt < *n)
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
            /* Compute space needed for DBDSQR */
            mnthr = ilaenv_(&c__6, "DGESVD", ch__1, m, n, &c__0, &c__0);
            bdspac = *n * 5;
            /* Compute space needed for DGEQRF */
            dgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dgeqrf__ = (integer) dum[0];
            /* Compute space needed for DORGQR */
            dorgqr_(m, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorgqr_n__ = (integer) dum[0];
            dorgqr_(m, m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorgqr_m__ = (integer) dum[0];
            /* Compute space needed for DGEBRD */
            dgebrd_(n, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1, &ierr);
            lwork_dgebrd__ = (integer) dum[0];
            /* Compute space needed for DORGBR P */
            dorgbr_("P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorgbr_p__ = (integer) dum[0];
            /* Compute space needed for DORGBR Q */
            dorgbr_("Q", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorgbr_q__ = (integer) dum[0];
            if (*m >= mnthr)
            {
                if (wntun)
                {
                    /* Path 1 (M much larger than N, JOBU='N') */
                    maxwrk = *n + lwork_dgeqrf__;
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = *n * 3 + lwork_dgebrd__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    if (wntvo || wntvas)
                    {
                        /* Computing MAX */
                        i__2 = maxwrk;
                        i__3 = *n * 3 + lwork_dorgbr_p__; // , expr subst
                        maxwrk = fla_max(i__2,i__3);
                    }
                    maxwrk = fla_max(maxwrk,bdspac);
                    /* Computing MAX */
                    i__2 = *n << 2;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntuo && wntvn)
                {
                    /* Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
                    wrkbl = *n + lwork_dgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_dorgqr_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    /* Computing MAX */
                    i__2 = *n * *n + wrkbl;
                    i__3 = *n * *n + *m * *n + *n; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntuo && wntvas)
                {
                    /* Path 3 (M much larger than N, JOBU='O', JOBVT='S' or */
                    /* 'A') */
                    wrkbl = *n + lwork_dgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_dorgqr_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    /* Computing MAX */
                    i__2 = *n * *n + wrkbl;
                    i__3 = *n * *n + *m * *n + *n; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntus && wntvn)
                {
                    /* Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
                    wrkbl = *n + lwork_dgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_dorgqr_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = *n * *n + wrkbl;
                    /* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntus && wntvo)
                {
                    /* Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
                    wrkbl = *n + lwork_dgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_dorgqr_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = (*n << 1) * *n + wrkbl;
                    /* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntus && wntvas)
                {
                    /* Path 6 (M much larger than N, JOBU='S', JOBVT='S' or */
                    /* 'A') */
                    wrkbl = *n + lwork_dgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_dorgqr_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = *n * *n + wrkbl;
                    /* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntua && wntvn)
                {
                    /* Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
                    wrkbl = *n + lwork_dgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_dorgqr_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = *n * *n + wrkbl;
                    /* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntua && wntvo)
                {
                    /* Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
                    wrkbl = *n + lwork_dgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_dorgqr_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = (*n << 1) * *n + wrkbl;
                    /* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntua && wntvas)
                {
                    /* Path 9 (M much larger than N, JOBU='A', JOBVT='S' or */
                    /* 'A') */
                    wrkbl = *n + lwork_dgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_dorgqr_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = *n * *n + wrkbl;
                    /* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = fla_max(i__2,bdspac);
                }
            }
            else
            {
                /* Path 10 (M at least N, but not much larger) */
                dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, & c_n1, &ierr);
                lwork_dgebrd__ = (integer) dum[0];
                maxwrk = *n * 3 + lwork_dgebrd__;
                if (wntus || wntuo)
                {
                    dorgbr_("Q", m, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
                    lwork_dorgbr_q__ = (integer) dum[0];
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = *n * 3 + lwork_dorgbr_q__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                }
                if (wntua)
                {
                    dorgbr_("Q", m, m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
                    lwork_dorgbr_q__ = (integer) dum[0];
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = *n * 3 + lwork_dorgbr_q__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                }
                if (! wntvn)
                {
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = *n * 3 + lwork_dorgbr_p__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                }
                maxwrk = fla_max(maxwrk,bdspac);
                /* Computing MAX */
                i__2 = *n * 3 + *m;
                minwrk = fla_max(i__2,bdspac);
            }
        }
        else if (minmn > 0)
        {
            /* Compute space needed for DBDSQR */
            mnthr = ilaenv_(&c__6, "DGESVD", ch__1, m, n, &c__0, &c__0);
            bdspac = *m * 5;
            /* Compute space needed for DGELQF */
            dgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dgelqf__ = (integer) dum[0];
            /* Compute space needed for DORGLQ */
            dorglq_(n, n, m, dum, n, dum, dum, &c_n1, &ierr);
            lwork_dorglq_n__ = (integer) dum[0];
            dorglq_(m, n, m, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorglq_m__ = (integer) dum[0];
            /* Compute space needed for DGEBRD */
            dgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1, &ierr);
            lwork_dgebrd__ = (integer) dum[0];
            /* Compute space needed for DORGBR P */
            dorgbr_("P", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr);
            lwork_dorgbr_p__ = (integer) dum[0];
            /* Compute space needed for DORGBR Q */
            dorgbr_("Q", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr);
            lwork_dorgbr_q__ = (integer) dum[0];
            if (*n >= mnthr)
            {
                if (wntvn)
                {
                    /* Path 1t(N much larger than M, JOBVT='N') */
                    maxwrk = *m + lwork_dgelqf__;
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = *m * 3 + lwork_dgebrd__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    if (wntuo || wntuas)
                    {
                        /* Computing MAX */
                        i__2 = maxwrk;
                        i__3 = *m * 3 + lwork_dorgbr_q__; // , expr subst
                        maxwrk = fla_max(i__2,i__3);
                    }
                    maxwrk = fla_max(maxwrk,bdspac);
                    /* Computing MAX */
                    i__2 = *m << 2;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntvo && wntun)
                {
                    /* Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
                    wrkbl = *m + lwork_dgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_dorglq_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    /* Computing MAX */
                    i__2 = *m * *m + wrkbl;
                    i__3 = *m * *m + *m * *n + *m; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntvo && wntuas)
                {
                    /* Path 3t(N much larger than M, JOBU='S' or 'A', */
                    /* JOBVT='O') */
                    wrkbl = *m + lwork_dgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_dorglq_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    /* Computing MAX */
                    i__2 = *m * *m + wrkbl;
                    i__3 = *m * *m + *m * *n + *m; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntvs && wntun)
                {
                    /* Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
                    wrkbl = *m + lwork_dgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_dorglq_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = *m * *m + wrkbl;
                    /* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntvs && wntuo)
                {
                    /* Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
                    wrkbl = *m + lwork_dgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_dorglq_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = (*m << 1) * *m + wrkbl;
                    /* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntvs && wntuas)
                {
                    /* Path 6t(N much larger than M, JOBU='S' or 'A', */
                    /* JOBVT='S') */
                    wrkbl = *m + lwork_dgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_dorglq_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = *m * *m + wrkbl;
                    /* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntva && wntun)
                {
                    /* Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
                    wrkbl = *m + lwork_dgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_dorglq_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = *m * *m + wrkbl;
                    /* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntva && wntuo)
                {
                    /* Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
                    wrkbl = *m + lwork_dgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_dorglq_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = (*m << 1) * *m + wrkbl;
                    /* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = fla_max(i__2,bdspac);
                }
                else if (wntva && wntuas)
                {
                    /* Path 9t(N much larger than M, JOBU='S' or 'A', */
                    /* JOBVT='A') */
                    wrkbl = *m + lwork_dgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_dorglq_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m * 3 + lwork_dorgbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    wrkbl = fla_max(wrkbl,bdspac);
                    maxwrk = *m * *m + wrkbl;
                    /* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = fla_max(i__2,bdspac);
                }
            }
            else
            {
                /* Path 10t(N greater than M, but not much larger) */
                dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, & c_n1, &ierr);
                lwork_dgebrd__ = (integer) dum[0];
                maxwrk = *m * 3 + lwork_dgebrd__;
                if (wntvs || wntvo)
                {
                    /* Compute space needed for DORGBR P */
                    dorgbr_("P", m, n, m, &a[a_offset], n, dum, dum, &c_n1, & ierr);
                    lwork_dorgbr_p__ = (integer) dum[0];
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = *m * 3 + lwork_dorgbr_p__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                }
                if (wntva)
                {
                    dorgbr_("P", n, n, m, &a[a_offset], n, dum, dum, &c_n1, & ierr);
                    lwork_dorgbr_p__ = (integer) dum[0];
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = *m * 3 + lwork_dorgbr_p__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                }
                if (! wntun)
                {
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = *m * 3 + lwork_dorgbr_q__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                }
                maxwrk = fla_max(maxwrk,bdspac);
                /* Computing MAX */
                i__2 = *m * 3 + *n;
                minwrk = fla_max(i__2,bdspac);
            }
        }
        maxwrk = fla_max(maxwrk,minwrk);
        work[1] = (double) maxwrk;
        if (*lwork < minwrk && ! lquery)
        {
            *info = -13;
        }
    }
    if (*info != 0)
    {
        i__2 = -(*info);
        xerbla_("DGESVD", &i__2);
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
