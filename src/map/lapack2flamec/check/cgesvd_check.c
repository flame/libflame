#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static integer c__6 = 6;
static integer c__0 = 0;
static integer c_n1 = -1;

int cgesvd_check(char *jobu, char *jobvt, integer *m, integer *n, scomplex *a, integer *lda, float *s, scomplex *u, integer *ldu, scomplex * vt, integer *ldvt, scomplex *work, integer *lwork, float *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__2, i__3;
    char ch__1[2];
    /* Local variables */
    float dum[2];
    integer ierr, lwork_cgebrd__, lwork_cgelqf__, lwork_cgeqrf__;
    integer minmn, wrkbl, mnthr;
    logical wntua, wntva, wntun, wntuo, wntvn, wntvo, wntus, wntvs;
    extern int 
      cgebrd_(integer*, integer*, scomplex*, integer*, float*, float*, scomplex*, scomplex*, scomplex*, integer*, integer*),
      cgelqf_(integer*, integer*, scomplex*, integer*, scomplex*, scomplex*, integer*, integer* ), 
      cgeqrf_(integer*, integer*, scomplex*, integer*, scomplex*, scomplex*, integer*, integer* ),
      cungbr_(char*, integer*, integer*, integer*, scomplex*, integer*, scomplex*, scomplex*, integer*, integer*),
      cunglq_(integer*, integer*, integer*, scomplex*, integer*, scomplex*, scomplex*, integer*, integer*), 
      cungqr_(integer*, integer*, integer*, scomplex*, integer*, scomplex*, scomplex*, integer*, integer*);
    integer minwrk, maxwrk;
    logical lquery, wntuas, wntvas;
    integer lwork_cungbr_p__, lwork_cungbr_q__, lwork_cunglq_m__, lwork_cunglq_n__, lwork_cungqr_m__, lwork_cungqr_n__;

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
    if (! (wntua || wntus || wntuo || wntun))
    {
        *info = -1;
    }
    else if (! (wntva || wntvs || wntvo || wntvn) || wntvo && wntuo)
    {
        *info = -2;
    }
    else if (*m < 0)
    {
        *info = -3;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -6;
    }
    else if (*ldu < 1 || wntuas && *ldu < *m)
    {
        *info = -9;
    }
    else if (*ldvt < 1 || wntva && *ldvt < *n || wntvs && *ldvt < minmn)
    {
        *info = -11;
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* CWorkspace refers to scomplex workspace, and RWorkspace to */
    /* float workspace. NB refers to the optimal block size for the */
    /* immediately following subroutine, as returned by ILAENV.) */
    if (*info == 0)
    {
        minwrk = 1;
        maxwrk = 1;
        if (*m >= *n && minmn > 0)
        {
            /* Space needed for ZBDSQR is BDSPAC = 5*N */
            mnthr = ilaenv_(&c__6, "CGESVD", ch__1, m, n, &c__0, &c__0);
            /* Compute space needed for CGEQRF */
            cgeqrf_(m, n, &a[a_offset], lda, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cgeqrf__ = (integer)dum[0];
            /* Compute space needed for CUNGQR */
            cungqr_(m, n, n, &a[a_offset], lda, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cungqr_n__ = (integer)dum[0];
            cungqr_(m, m, n, &a[a_offset], lda, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cungqr_m__ = (integer)dum[0];
            /* Compute space needed for CGEBRD */
            cgebrd_(n, n, &a[a_offset], lda, &s[1], dum, (scomplex*)dum, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cgebrd__ = (integer)dum[0];
            /* Compute space needed for CUNGBR */
            cungbr_("P", n, n, n, &a[a_offset], lda, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cungbr_p__ = (integer)dum[0];
            cungbr_("Q", n, n, n, &a[a_offset], lda, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cungbr_q__ = (integer)dum[0];
            mnthr = ilaenv_(&c__6, "CGESVD", ch__1, m, n, &c__0, &c__0);
            if (*m >= mnthr)
            {
                if (wntun)
                {
                    /* Path 1 (M much larger than N, JOBU='N') */
                    maxwrk = *n + lwork_cgeqrf__;
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = (*n << 1) + lwork_cgebrd__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    if (wntvo || wntvas)
                    {
                        /* Computing MAX */
                        i__2 = maxwrk;
                        i__3 = (*n << 1) + lwork_cungbr_p__; // , expr subst
                        maxwrk = fla_max(i__2,i__3);
                    }
                    minwrk = *n * 3;
                }
                else if (wntuo && wntvn)
                {
                    /* Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
                    wrkbl = *n + lwork_cgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_cungqr_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = *n * *n + wrkbl;
                    i__3 = *n * *n + *m * *n; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    minwrk = (*n << 1) + *m;
                }
                else if (wntuo && wntvas)
                {
                    /* Path 3 (M much larger than N, JOBU='O', JOBVT='S' or */
                    /* 'A') */
                    wrkbl = *n + lwork_cgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_cungqr_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = *n * *n + wrkbl;
                    i__3 = *n * *n + *m * *n; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    minwrk = (*n << 1) + *m;
                }
                else if (wntus && wntvn)
                {
                    /* Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
                    wrkbl = *n + lwork_cgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_cungqr_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = *n * *n + wrkbl;
                    minwrk = (*n << 1) + *m;
                }
                else if (wntus && wntvo)
                {
                    /* Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
                    wrkbl = *n + lwork_cgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_cungqr_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = (*n << 1) * *n + wrkbl;
                    minwrk = (*n << 1) + *m;
                }
                else if (wntus && wntvas)
                {
                    /* Path 6 (M much larger than N, JOBU='S', JOBVT='S' or */
                    /* 'A') */
                    wrkbl = *n + lwork_cgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_cungqr_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = *n * *n + wrkbl;
                    minwrk = (*n << 1) + *m;
                }
                else if (wntua && wntvn)
                {
                    /* Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
                    wrkbl = *n + lwork_cgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_cungqr_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = *n * *n + wrkbl;
                    minwrk = (*n << 1) + *m;
                }
                else if (wntua && wntvo)
                {
                    /* Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
                    wrkbl = *n + lwork_cgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_cungqr_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = (*n << 1) * *n + wrkbl;
                    minwrk = (*n << 1) + *m;
                }
                else if (wntua && wntvas)
                {
                    /* Path 9 (M much larger than N, JOBU='A', JOBVT='S' or */
                    /* 'A') */
                    wrkbl = *n + lwork_cgeqrf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *n + lwork_cungqr_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*n << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = *n * *n + wrkbl;
                    minwrk = (*n << 1) + *m;
                }
            }
            else
            {
                /* Path 10 (M at least N, but not much larger) */
                cgebrd_(m, n, &a[a_offset], lda, &s[1], dum, (scomplex*)dum, (scomplex*)dum, (scomplex*)dum, & c_n1, &ierr);
                lwork_cgebrd__ = (integer)dum[0];
                maxwrk = (*n << 1) + lwork_cgebrd__;
                if (wntus || wntuo)
                {
                    cungbr_("Q", m, n, n, &a[a_offset], lda, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
                    lwork_cungbr_q__ = (integer)dum[0];
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = (*n << 1) + lwork_cungbr_q__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                }
                if (wntua)
                {
                    cungbr_("Q", m, m, n, &a[a_offset], lda, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
                    lwork_cungbr_q__ = (integer)dum[0];
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = (*n << 1) + lwork_cungbr_q__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                }
                if (! wntvn)
                {
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = (*n << 1) + lwork_cungbr_p__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    minwrk = (*n << 1) + *m;
                }
            }
        }
        else if (minmn > 0)
        {
            /* Space needed for CBDSQR is BDSPAC = 5*M */
            mnthr = ilaenv_(&c__6, "CGESVD", ch__1, m, n, &c__0, &c__0);
            /* Compute space needed for CGELQF */
            cgelqf_(m, n, &a[a_offset], lda, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cgelqf__ = (integer)dum[0];
            /* Compute space needed for CUNGLQ */
            cunglq_(n, n, m, (scomplex*)dum, n, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cunglq_n__ = (integer)dum[0];
            cunglq_(m, n, m, &a[a_offset], lda, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cunglq_m__ = (integer)dum[0];
            /* Compute space needed for CGEBRD */
            cgebrd_(m, m, &a[a_offset], lda, &s[1], dum, (scomplex*)dum, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cgebrd__ = (integer)dum[0];
            /* Compute space needed for CUNGBR P */
            cungbr_("P", m, m, m, &a[a_offset], n, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cungbr_p__ = (integer)dum[0];
            /* Compute space needed for CUNGBR Q */
            cungbr_("Q", m, m, m, &a[a_offset], n, (scomplex*)dum, (scomplex*)dum, &c_n1, &ierr);
            lwork_cungbr_q__ = (integer)dum[0];
            if (*n >= mnthr)
            {
                if (wntvn)
                {
                    /* Path 1t(N much larger than M, JOBVT='N') */
                    maxwrk = *m + lwork_cgelqf__;
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = (*m << 1) + lwork_cgebrd__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    if (wntuo || wntuas)
                    {
                        /* Computing MAX */
                        i__2 = maxwrk;
                        i__3 = (*m << 1) + lwork_cungbr_q__; // , expr subst
                        maxwrk = fla_max(i__2,i__3);
                    }
                    minwrk = *m * 3;
                }
                else if (wntvo && wntun)
                {
                    /* Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
                    wrkbl = *m + lwork_cgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_cunglq_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = *m * *m + wrkbl;
                    i__3 = *m * *m + *m * *n; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    minwrk = (*m << 1) + *n;
                }
                else if (wntvo && wntuas)
                {
                    /* Path 3t(N much larger than M, JOBU='S' or 'A', */
                    /* JOBVT='O') */
                    wrkbl = *m + lwork_cgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_cunglq_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = *m * *m + wrkbl;
                    i__3 = *m * *m + *m * *n; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    minwrk = (*m << 1) + *n;
                }
                else if (wntvs && wntun)
                {
                    /* Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
                    wrkbl = *m + lwork_cgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_cunglq_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = *m * *m + wrkbl;
                    minwrk = (*m << 1) + *n;
                }
                else if (wntvs && wntuo)
                {
                    /* Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
                    wrkbl = *m + lwork_cgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_cunglq_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = (*m << 1) * *m + wrkbl;
                    minwrk = (*m << 1) + *n;
                }
                else if (wntvs && wntuas)
                {
                    /* Path 6t(N much larger than M, JOBU='S' or 'A', */
                    /* JOBVT='S') */
                    wrkbl = *m + lwork_cgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_cunglq_m__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = *m * *m + wrkbl;
                    minwrk = (*m << 1) + *n;
                }
                else if (wntva && wntun)
                {
                    /* Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
                    wrkbl = *m + lwork_cgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_cunglq_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = *m * *m + wrkbl;
                    minwrk = (*m << 1) + *n;
                }
                else if (wntva && wntuo)
                {
                    /* Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
                    wrkbl = *m + lwork_cgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_cunglq_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = (*m << 1) * *m + wrkbl;
                    minwrk = (*m << 1) + *n;
                }
                else if (wntva && wntuas)
                {
                    /* Path 9t(N much larger than M, JOBU='S' or 'A', */
                    /* JOBVT='A') */
                    wrkbl = *m + lwork_cgelqf__;
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = *m + lwork_cunglq_n__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cgebrd__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_p__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    /* Computing MAX */
                    i__2 = wrkbl;
                    i__3 = (*m << 1) + lwork_cungbr_q__; // , expr subst
                    wrkbl = fla_max(i__2,i__3);
                    maxwrk = *m * *m + wrkbl;
                    minwrk = (*m << 1) + *n;
                }
            }
            else
            {
                /* Path 10t(N greater than M, but not much larger) */
                cgebrd_(m, n, &a[a_offset], lda, &s[1], dum, (scomplex*)dum, (scomplex*)dum, (scomplex*)dum, & c_n1, &ierr);
                lwork_cgebrd__ = (integer)dum[0];
                maxwrk = (*m << 1) + lwork_cgebrd__;
                if (wntvs || wntvo)
                {
                    /* Compute space needed for CUNGBR P */
                    cungbr_("P", m, n, m, &a[a_offset], n, (scomplex*)dum, (scomplex*)dum, &c_n1, & ierr);
                    lwork_cungbr_p__ = (integer)dum[0];
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = (*m << 1) + lwork_cungbr_p__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                }
                if (wntva)
                {
                    cungbr_("P", n, n, m, &a[a_offset], n, (scomplex*)dum, (scomplex*)dum, &c_n1, & ierr);
                    lwork_cungbr_p__ = (integer)dum[0];
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = (*m << 1) + lwork_cungbr_p__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                }
                if (! wntun)
                {
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = (*m << 1) + lwork_cungbr_q__; // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    minwrk = (*m << 1) + *n;
                }
            }
        }
        maxwrk = fla_max(minwrk,maxwrk);
        work[1].real = (float) maxwrk;
        work[1].imag = 0.f; // , expr subst
        if (*lwork < minwrk && ! lquery)
        {
            *info = -13;
        }
    }
    if (*info != 0)
    {
        i__2 = -(*info);
        xerbla_("CGESVD", &i__2);
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
