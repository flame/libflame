/* ../netlib/slaln2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLALN2 solves a 1-by-1 or 2-by-2 linear system of equations of the specified form. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLALN2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaln2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaln2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaln2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, */
/* LDB, WR, WI, X, LDX, SCALE, XNORM, INFO ) */
/* .. Scalar Arguments .. */
/* LOGICAL LTRANS */
/* INTEGER INFO, LDA, LDB, LDX, NA, NW */
/* REAL CA, D1, D2, SCALE, SMIN, WI, WR, XNORM */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLALN2 solves a system of the form (ca A - w D ) X = s B */
/* > or (ca A**T - w D) X = s B with possible scaling ("s") and */
/* > perturbation of A. (A**T means A-transpose.) */
/* > */
/* > A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA */
/* > real diagonal matrix, w is a real or complex value, and X and B are */
/* > NA x 1 matrices -- real if w is real, complex if w is complex. NA */
/* > may be 1 or 2. */
/* > */
/* > If w is complex, X and B are represented as NA x 2 matrices, */
/* > the first column of each being the real part and the second */
/* > being the imaginary part. */
/* > */
/* > "s" is a scaling factor (.LE. 1), computed by SLALN2, which is */
/* > so chosen that X can be computed without overflow. X is further */
/* > scaled if necessary to assure that norm(ca A - w D)*norm(X) is less */
/* > than overflow. */
/* > */
/* > If both singular values of (ca A - w D) are less than SMIN, */
/* > SMIN*identity will be used instead of (ca A - w D). If only one */
/* > singular value is less than SMIN, one element of (ca A - w D) will be */
/* > perturbed enough to make the smallest singular value roughly SMIN. */
/* > If both singular values are at least SMIN, (ca A - w D) will not be */
/* > perturbed. In any case, the perturbation will be at most some small */
/* > multiple of fla_max( SMIN, ulp*norm(ca A - w D) ). The singular values */
/* > are computed by infinity-norm approximations, and thus will only be */
/* > correct to a factor of 2 or so. */
/* > */
/* > Note: all input quantities are assumed to be smaller than overflow */
/* > by a reasonable factor. (See BIGNUM.) */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] LTRANS */
/* > \verbatim */
/* > LTRANS is LOGICAL */
/* > =.TRUE.: A-transpose will be used. */
/* > =.FALSE.: A will be used (not transposed.) */
/* > \endverbatim */
/* > */
/* > \param[in] NA */
/* > \verbatim */
/* > NA is INTEGER */
/* > The size of the matrix A. It may (only) be 1 or 2. */
/* > \endverbatim */
/* > */
/* > \param[in] NW */
/* > \verbatim */
/* > NW is INTEGER */
/* > 1 if "w" is real, 2 if "w" is complex. It may only be 1 */
/* > or 2. */
/* > \endverbatim */
/* > */
/* > \param[in] SMIN */
/* > \verbatim */
/* > SMIN is REAL */
/* > The desired lower bound on the singular values of A. This */
/* > should be a safe distance away from underflow or overflow, */
/* > say, between (underflow/machine precision) and (machine */
/* > precision * overflow ). (See BIGNUM and ULP.) */
/* > \endverbatim */
/* > */
/* > \param[in] CA */
/* > \verbatim */
/* > CA is REAL */
/* > The coefficient c, which A is multiplied by. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,NA) */
/* > The NA x NA matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of A. It must be at least NA. */
/* > \endverbatim */
/* > */
/* > \param[in] D1 */
/* > \verbatim */
/* > D1 is REAL */
/* > The 1,1 element in the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] D2 */
/* > \verbatim */
/* > D2 is REAL */
/* > The 2,2 element in the diagonal matrix D. Not used if NW=1. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,NW) */
/* > The NA x NW matrix B (right-hand side). If NW=2 ("w" is */
/* > complex), column 1 contains the real part of B and column 2 */
/* > contains the imaginary part. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B. It must be at least NA. */
/* > \endverbatim */
/* > */
/* > \param[in] WR */
/* > \verbatim */
/* > WR is REAL */
/* > The real part of the scalar "w". */
/* > \endverbatim */
/* > */
/* > \param[in] WI */
/* > \verbatim */
/* > WI is REAL */
/* > The imaginary part of the scalar "w". Not used if NW=1. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* > X is REAL array, dimension (LDX,NW) */
/* > The NA x NW matrix X (unknowns), as computed by SLALN2. */
/* > If NW=2 ("w" is complex), on exit, column 1 will contain */
/* > the real part of X and column 2 will contain the imaginary */
/* > part. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of X. It must be at least NA. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is REAL */
/* > The scale factor that B must be multiplied by to insure */
/* > that overflow does not occur when computing X. Thus, */
/* > (ca A - w D) X will be SCALE*B, not B (ignoring */
/* > perturbations of A.) It will be at most 1. */
/* > \endverbatim */
/* > */
/* > \param[out] XNORM */
/* > \verbatim */
/* > XNORM is REAL */
/* > The infinity-norm of X, when X is regarded as an NA x NW */
/* > real matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > An error flag. It will be set to zero if no error occurs, */
/* > a negative number if an argument is in error, or a positive */
/* > number if ca A - w D had to be perturbed. */
/* > The possible values are: */
/* > = 0: No error occurred, and (ca A - w D) did not have to be */
/* > perturbed. */
/* > = 1: (ca A - w D) had to be perturbed to make its smallest */
/* > (or only) singular value greater than SMIN. */
/* > NOTE: In the interests of speed, this routine does not */
/* > check the inputs for errors. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */

/* Set FP_contract to off as floating point contraction instructions */
/* like fma is seen to lead precision issues in this function and */
/* causing few netlib tests to fail */
#pragma STDC FP_CONTRACT OFF

int slaln2_(logical *ltrans, integer *na, integer *nw, real * smin, real *ca, real *a, integer *lda, real *d1, real *d2, real *b, integer *ldb, real *wr, real *wi, real *x, integer *ldx, real *scale, real *xnorm, integer *info)
{
    /* Initialized data */
    static const logical cswap[4] =
    {
        FALSE_,FALSE_,TRUE_,TRUE_
    }
    ;
    static const logical rswap[4] =
    {
        FALSE_,TRUE_,FALSE_,TRUE_
    }
    ;
    static const integer ipivot[16] /* was [4][4] */
    =
    {
        1,2,3,4,2,1,4,3,3,4,1,2, 4,3,2,1
    }
    ;
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    real equiv_0[4], equiv_1[4];
    /* Local variables */
    integer j;
#define ci (equiv_0)
#define cr (equiv_1)
    real bi1, bi2, br1, br2, xi1, xi2, xr1, xr2, ci21, ci22, cr21, cr22, li21, csi, ui11, lr21, ui12, ui22;
#define civ (equiv_0)
    real csr, ur11, ur12, ur22;
#define crv (equiv_1)
    real bbnd, cmax, ui11r, ui12s, temp, ur11r, ur12s, u22abs;
    integer icmax;
    real bnorm, cnorm, smini;
    extern real slamch_(char *);
    real bignum;
    extern /* Subroutine */
    int sladiv_(real *, real *, real *, real *, real *, real *);
    real smlnum;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Equivalences .. */
    /* .. */
    /* .. Data statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    /* Function Body */
    /* .. */
    /* .. Executable Statements .. */
    /* Compute BIGNUM */
    smlnum = 2.f * slamch_("Safe minimum");
    bignum = 1.f / smlnum;
    smini = fla_max(*smin,smlnum);
    /* Don't check for input errors */
    *info = 0;
    /* Standard Initializations */
    *scale = 1.f;
    if (*na == 1)
    {
        /* 1 x 1 (i.e., scalar) system C X = B */
        if (*nw == 1)
        {
            /* Real 1x1 system. */
            /* C = ca A - w D */
            csr = *ca * a[a_dim1 + 1] - *wr * *d1;
            cnorm = f2c_abs(csr);
            /* If | C | < SMINI, use C = SMINI */
            if (cnorm < smini)
            {
                csr = smini;
                cnorm = smini;
                *info = 1;
            }
            /* Check scaling for X = B / C */
            bnorm = (r__1 = b[b_dim1 + 1], f2c_abs(r__1));
            if (cnorm < 1.f && bnorm > 1.f)
            {
                if (bnorm > bignum * cnorm)
                {
                    *scale = 1.f / bnorm;
                }
            }
            /* Compute X */
            x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / csr;
            *xnorm = (r__1 = x[x_dim1 + 1], f2c_abs(r__1));
        }
        else
        {
            /* Complex 1x1 system (w is complex) */
            /* C = ca A - w D */
            csr = *ca * a[a_dim1 + 1] - *wr * *d1;
            csi = -(*wi) * *d1;
            cnorm = f2c_abs(csr) + f2c_abs(csi);
            /* If | C | < SMINI, use C = SMINI */
            if (cnorm < smini)
            {
                csr = smini;
                csi = 0.f;
                cnorm = smini;
                *info = 1;
            }
            /* Check scaling for X = B / C */
            bnorm = (r__1 = b[b_dim1 + 1], f2c_abs(r__1)) + (r__2 = b[(b_dim1 << 1) + 1], f2c_abs(r__2));
            if (cnorm < 1.f && bnorm > 1.f)
            {
                if (bnorm > bignum * cnorm)
                {
                    *scale = 1.f / bnorm;
                }
            }
            /* Compute X */
            r__1 = *scale * b[b_dim1 + 1];
            r__2 = *scale * b[(b_dim1 << 1) + 1];
            sladiv_(&r__1, &r__2, &csr, &csi, &x[x_dim1 + 1], &x[(x_dim1 << 1) + 1]);
            *xnorm = (r__1 = x[x_dim1 + 1], f2c_abs(r__1)) + (r__2 = x[(x_dim1 << 1) + 1], f2c_abs(r__2));
        }
    }
    else
    {
        /* 2x2 System */
        /* Compute the real part of C = ca A - w D (or ca A**T - w D ) */
        cr[0] = *ca * a[a_dim1 + 1] - *wr * *d1;
        cr[3] = *ca * a[(a_dim1 << 1) + 2] - *wr * *d2;
        if (*ltrans)
        {
            cr[2] = *ca * a[a_dim1 + 2];
            cr[1] = *ca * a[(a_dim1 << 1) + 1];
        }
        else
        {
            cr[1] = *ca * a[a_dim1 + 2];
            cr[2] = *ca * a[(a_dim1 << 1) + 1];
        }
        if (*nw == 1)
        {
            /* Real 2x2 system (w is real) */
            /* Find the largest element in C */
            cmax = 0.f;
            icmax = 0;
            for (j = 1;
                    j <= 4;
                    ++j)
            {
                if ((r__1 = crv[j - 1], f2c_abs(r__1)) > cmax)
                {
                    cmax = (r__1 = crv[j - 1], f2c_abs(r__1));
                    icmax = j;
                }
                /* L10: */
            }
            /* If norm(C) < SMINI, use SMINI*identity. */
            if (cmax < smini)
            {
                /* Computing MAX */
                r__3 = (r__1 = b[b_dim1 + 1], f2c_abs(r__1));
                r__4 = (r__2 = b[ b_dim1 + 2], f2c_abs(r__2)); // , expr subst
                bnorm = fla_max(r__3,r__4);
                if (smini < 1.f && bnorm > 1.f)
                {
                    if (bnorm > bignum * smini)
                    {
                        *scale = 1.f / bnorm;
                    }
                }
                temp = *scale / smini;
                x[x_dim1 + 1] = temp * b[b_dim1 + 1];
                x[x_dim1 + 2] = temp * b[b_dim1 + 2];
                *xnorm = temp * bnorm;
                *info = 1;
                return 0;
            }
            /* Gaussian elimination with complete pivoting. */
            ur11 = crv[icmax - 1];
            cr21 = crv[ipivot[(icmax << 2) - 3] - 1];
            ur12 = crv[ipivot[(icmax << 2) - 2] - 1];
            cr22 = crv[ipivot[(icmax << 2) - 1] - 1];
            ur11r = 1.f / ur11;
            lr21 = ur11r * cr21;
            ur22 = cr22 - ur12 * lr21;
            /* If smaller pivot < SMINI, use SMINI */
            if (f2c_abs(ur22) < smini)
            {
                ur22 = smini;
                *info = 1;
            }
            if (rswap[icmax - 1])
            {
                br1 = b[b_dim1 + 2];
                br2 = b[b_dim1 + 1];
            }
            else
            {
                br1 = b[b_dim1 + 1];
                br2 = b[b_dim1 + 2];
            }
            br2 -= lr21 * br1;
            /* Computing MAX */
            r__2 = (r__1 = br1 * (ur22 * ur11r), f2c_abs(r__1));
            r__3 = f2c_abs(br2); // , expr subst
            bbnd = fla_max(r__2,r__3);
            if (bbnd > 1.f && f2c_abs(ur22) < 1.f)
            {
                if (bbnd >= bignum * f2c_abs(ur22))
                {
                    *scale = 1.f / bbnd;
                }
            }
            xr2 = br2 * *scale / ur22;
            xr1 = *scale * br1 * ur11r - xr2 * (ur11r * ur12);
            if (cswap[icmax - 1])
            {
                x[x_dim1 + 1] = xr2;
                x[x_dim1 + 2] = xr1;
            }
            else
            {
                x[x_dim1 + 1] = xr1;
                x[x_dim1 + 2] = xr2;
            }
            /* Computing MAX */
            r__1 = f2c_abs(xr1);
            r__2 = f2c_abs(xr2); // , expr subst
            *xnorm = fla_max(r__1,r__2);
            /* Further scaling if norm(A) norm(X) > overflow */
            if (*xnorm > 1.f && cmax > 1.f)
            {
                if (*xnorm > bignum / cmax)
                {
                    temp = cmax / bignum;
                    x[x_dim1 + 1] = temp * x[x_dim1 + 1];
                    x[x_dim1 + 2] = temp * x[x_dim1 + 2];
                    *xnorm = temp * *xnorm;
                    *scale = temp * *scale;
                }
            }
        }
        else
        {
            /* Complex 2x2 system (w is complex) */
            /* Find the largest element in C */
            ci[0] = -(*wi) * *d1;
            ci[1] = 0.f;
            ci[2] = 0.f;
            ci[3] = -(*wi) * *d2;
            cmax = 0.f;
            icmax = 0;
            for (j = 1;
                    j <= 4;
                    ++j)
            {
                if ((r__1 = crv[j - 1], f2c_abs(r__1)) + (r__2 = civ[j - 1], f2c_abs( r__2)) > cmax)
                {
                    cmax = (r__1 = crv[j - 1], f2c_abs(r__1)) + (r__2 = civ[j - 1], f2c_abs(r__2));
                    icmax = j;
                }
                /* L20: */
            }
            /* If norm(C) < SMINI, use SMINI*identity. */
            if (cmax < smini)
            {
                /* Computing MAX */
                r__5 = (r__1 = b[b_dim1 + 1], f2c_abs(r__1)) + (r__2 = b[(b_dim1 << 1) + 1], f2c_abs(r__2));
                r__6 = (r__3 = b[b_dim1 + 2], f2c_abs(r__3)) + (r__4 = b[(b_dim1 << 1) + 2], f2c_abs(r__4)); // , expr subst
                bnorm = fla_max(r__5,r__6);
                if (smini < 1.f && bnorm > 1.f)
                {
                    if (bnorm > bignum * smini)
                    {
                        *scale = 1.f / bnorm;
                    }
                }
                temp = *scale / smini;
                x[x_dim1 + 1] = temp * b[b_dim1 + 1];
                x[x_dim1 + 2] = temp * b[b_dim1 + 2];
                x[(x_dim1 << 1) + 1] = temp * b[(b_dim1 << 1) + 1];
                x[(x_dim1 << 1) + 2] = temp * b[(b_dim1 << 1) + 2];
                *xnorm = temp * bnorm;
                *info = 1;
                return 0;
            }
            /* Gaussian elimination with complete pivoting. */
            ur11 = crv[icmax - 1];
            ui11 = civ[icmax - 1];
            cr21 = crv[ipivot[(icmax << 2) - 3] - 1];
            ci21 = civ[ipivot[(icmax << 2) - 3] - 1];
            ur12 = crv[ipivot[(icmax << 2) - 2] - 1];
            ui12 = civ[ipivot[(icmax << 2) - 2] - 1];
            cr22 = crv[ipivot[(icmax << 2) - 1] - 1];
            ci22 = civ[ipivot[(icmax << 2) - 1] - 1];
            if (icmax == 1 || icmax == 4)
            {
                /* Code when off-diagonals of pivoted C are real */
                if (f2c_abs(ur11) > f2c_abs(ui11))
                {
                    temp = ui11 / ur11;
                    /* Computing 2nd power */
                    r__1 = temp;
                    ur11r = 1.f / (ur11 * (r__1 * r__1 + 1.f));
                    ui11r = -temp * ur11r;
                }
                else
                {
                    temp = ur11 / ui11;
                    /* Computing 2nd power */
                    r__1 = temp;
                    ui11r = -1.f / (ui11 * (r__1 * r__1 + 1.f));
                    ur11r = -temp * ui11r;
                }
                lr21 = cr21 * ur11r;
                li21 = cr21 * ui11r;
                ur12s = ur12 * ur11r;
                ui12s = ur12 * ui11r;
                ur22 = cr22 - ur12 * lr21;
                ui22 = ci22 - ur12 * li21;
            }
            else
            {
                /* Code when diagonals of pivoted C are real */
                ur11r = 1.f / ur11;
                ui11r = 0.f;
                lr21 = cr21 * ur11r;
                li21 = ci21 * ur11r;
                ur12s = ur12 * ur11r;
                ui12s = ui12 * ur11r;
                ur22 = cr22 - ur12 * lr21 + ui12 * li21;
                ui22 = -ur12 * li21 - ui12 * lr21;
            }
            u22abs = f2c_abs(ur22) + f2c_abs(ui22);
            /* If smaller pivot < SMINI, use SMINI */
            if (u22abs < smini)
            {
                ur22 = smini;
                ui22 = 0.f;
                *info = 1;
            }
            if (rswap[icmax - 1])
            {
                br2 = b[b_dim1 + 1];
                br1 = b[b_dim1 + 2];
                bi2 = b[(b_dim1 << 1) + 1];
                bi1 = b[(b_dim1 << 1) + 2];
            }
            else
            {
                br1 = b[b_dim1 + 1];
                br2 = b[b_dim1 + 2];
                bi1 = b[(b_dim1 << 1) + 1];
                bi2 = b[(b_dim1 << 1) + 2];
            }
            br2 = br2 - lr21 * br1 + li21 * bi1;
            bi2 = bi2 - li21 * br1 - lr21 * bi1;
            /* Computing MAX */
            r__1 = (f2c_abs(br1) + f2c_abs(bi1)) * (u22abs * (f2c_abs(ur11r) + f2c_abs(ui11r)) );
            r__2 = f2c_abs(br2) + f2c_abs(bi2); // , expr subst
            bbnd = fla_max(r__1,r__2);
            if (bbnd > 1.f && u22abs < 1.f)
            {
                if (bbnd >= bignum * u22abs)
                {
                    *scale = 1.f / bbnd;
                    br1 = *scale * br1;
                    bi1 = *scale * bi1;
                    br2 = *scale * br2;
                    bi2 = *scale * bi2;
                }
            }
            sladiv_(&br2, &bi2, &ur22, &ui22, &xr2, &xi2);
            xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
            xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
            if (cswap[icmax - 1])
            {
                x[x_dim1 + 1] = xr2;
                x[x_dim1 + 2] = xr1;
                x[(x_dim1 << 1) + 1] = xi2;
                x[(x_dim1 << 1) + 2] = xi1;
            }
            else
            {
                x[x_dim1 + 1] = xr1;
                x[x_dim1 + 2] = xr2;
                x[(x_dim1 << 1) + 1] = xi1;
                x[(x_dim1 << 1) + 2] = xi2;
            }
            /* Computing MAX */
            r__1 = f2c_abs(xr1) + f2c_abs(xi1);
            r__2 = f2c_abs(xr2) + f2c_abs(xi2); // , expr subst
            *xnorm = fla_max(r__1,r__2);
            /* Further scaling if norm(A) norm(X) > overflow */
            if (*xnorm > 1.f && cmax > 1.f)
            {
                if (*xnorm > bignum / cmax)
                {
                    temp = cmax / bignum;
                    x[x_dim1 + 1] = temp * x[x_dim1 + 1];
                    x[x_dim1 + 2] = temp * x[x_dim1 + 2];
                    x[(x_dim1 << 1) + 1] = temp * x[(x_dim1 << 1) + 1];
                    x[(x_dim1 << 1) + 2] = temp * x[(x_dim1 << 1) + 2];
                    *xnorm = temp * *xnorm;
                    *scale = temp * *scale;
                }
            }
        }
    }
    return 0;
    /* End of SLALN2 */
}
/* slaln2_ */
#undef crv
#undef civ
#undef cr
#undef ci
