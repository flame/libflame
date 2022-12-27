#ifdef FLA_ENABLE_XBLAS
/* ../netlib/cla_gerfsx_extended.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static complex c_b6 =
{
    -1.f,0.f
    }
;
static complex c_b8 =
{
    1.f,0.f
}
;
static real c_b31 = 1.f;
/* > \brief \b CLA_GERFSX_EXTENDED */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLA_GERFSX_EXTENDED + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_ger fsx_extended.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_ger fsx_extended.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_ger fsx_extended.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, NRHS, A, */
/* LDA, AF, LDAF, IPIV, COLEQU, C, B, */
/* LDB, Y, LDY, BERR_OUT, N_NORMS, */
/* ERRS_N, ERRS_C, RES, AYB, DY, */
/* Y_TAIL, RCOND, ITHRESH, RTHRESH, */
/* DZ_UB, IGNORE_CWISE, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE, */
/* $ TRANS_TYPE, N_NORMS */
/* LOGICAL COLEQU, IGNORE_CWISE */
/* INTEGER ITHRESH */
/* REAL RTHRESH, DZ_UB */
/* .. */
/* .. Array Arguments */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/* $ Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * ) */
/* REAL C( * ), AYB( * ), RCOND, BERR_OUT( * ), */
/* $ ERRS_N( NRHS, * ), ERRS_C( NRHS, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > CLA_GERFSX_EXTENDED improves the computed solution to a system of */
/* > linear equations by performing extra-precise iterative refinement */
/* > and provides error bounds and backward error estimates for the solution. */
/* > This subroutine is called by CGERFSX to perform iterative refinement. */
/* > In addition to normwise error bound, the code provides maximum */
/* > componentwise error bound if possible. See comments for ERRS_N */
/* > and ERRS_C for details of the error bounds. Note that this */
/* > subroutine is only resonsible for setting the second fields of */
/* > ERRS_N and ERRS_C. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] PREC_TYPE */
/* > \verbatim */
/* > PREC_TYPE is INTEGER */
/* > Specifies the intermediate precision to be used in refinement. */
/* > The value is defined by ILAPREC(P) where P is a CHARACTER and */
/* > P = 'S': Single */
/* > = 'D': Double */
/* > = 'I': Indigenous */
/* > = 'X', 'E': Extra */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS_TYPE */
/* > \verbatim */
/* > TRANS_TYPE is INTEGER */
/* > Specifies the transposition operation on A. */
/* > The value is defined by ILATRANS(T) where T is a CHARACTER and */
/* > T = 'N': No transpose */
/* > = 'T': Transpose */
/* > = 'C': Conjugate transpose */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right-hand-sides, i.e., the number of columns of the */
/* > matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* > AF is COMPLEX array, dimension (LDAF,N) */
/* > The factors L and U from the factorization */
/* > A = P*L*U as computed by CGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices from the factorization A = P*L*U */
/* > as computed by CGETRF;
row i of the matrix was interchanged */
/* > with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] COLEQU */
/* > \verbatim */
/* > COLEQU is LOGICAL */
/* > If .TRUE. then column equilibration was done to A before calling */
/* > this routine. This is needed to compute the solution and error */
/* > bounds correctly. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is REAL array, dimension (N) */
/* > The column scale factors for A. If COLEQU = .FALSE., C */
/* > is not accessed. If C is input, each element of C should be a power */
/* > of the radix to ensure a reliable solution and error estimates. */
/* > Scaling by powers of the radix does not cause rounding errors unless */
/* > the result underflows or overflows. Rounding errors during scaling */
/* > lead to refining with a matrix that is not equivalent to the */
/* > input matrix, producing error estimates that may not be */
/* > reliable. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > The right-hand-side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is COMPLEX array, dimension (LDY,NRHS) */
/* > On entry, the solution matrix X, as computed by CGETRS. */
/* > On exit, the improved solution matrix Y. */
/* > \endverbatim */
/* > */
/* > \param[in] LDY */
/* > \verbatim */
/* > LDY is INTEGER */
/* > The leading dimension of the array Y. LDY >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] BERR_OUT */
/* > \verbatim */
/* > BERR_OUT is REAL array, dimension (NRHS) */
/* > On exit, BERR_OUT(j) contains the componentwise relative backward */
/* > error for right-hand-side j from the formula */
/* > fla_max(i) ( f2c_abs(RES(i)) / ( f2c_abs(op(A_s))*f2c_abs(Y) + f2c_abs(B_s) )(i) ) */
/* > where f2c_abs(Z) is the componentwise absolute value of the matrix */
/* > or vector Z. This is computed by CLA_LIN_BERR. */
/* > \endverbatim */
/* > */
/* > \param[in] N_NORMS */
/* > \verbatim */
/* > N_NORMS is INTEGER */
/* > Determines which error bounds to return (see ERRS_N */
/* > and ERRS_C). */
/* > If N_NORMS >= 1 return normwise error bounds. */
/* > If N_NORMS >= 2 return componentwise error bounds. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ERRS_N */
/* > \verbatim */
/* > ERRS_N is REAL array, dimension (NRHS, N_ERR_BNDS) */
/* > For each right-hand side, this array contains information about */
/* > various error bounds and condition numbers corresponding to the */
/* > normwise relative error, which is defined as follows: */
/* > */
/* > Normwise relative error in the ith solution vector: */
/* > max_j (f2c_abs(XTRUE(j,i) - X(j,i))) */
/* > ------------------------------ */
/* > max_j f2c_abs(X(j,i)) */
/* > */
/* > The array is indexed by the type of error information as described */
/* > below. There currently are up to three pieces of information */
/* > returned. */
/* > */
/* > The first index in ERRS_N(i,:) corresponds to the ith */
/* > right-hand side. */
/* > */
/* > The second index in ERRS_N(:,err) contains the following */
/* > three fields: */
/* > err = 1 "Trust/don't trust" boolean. Trust the answer if the */
/* > reciprocal condition number is less than the threshold */
/* > sqrt(n) * slamch('Epsilon'). */
/* > */
/* > err = 2 "Guaranteed" error bound: The estimated forward error, */
/* > almost certainly within a factor of 10 of the true error */
/* > so long as the next entry is greater than the threshold */
/* > sqrt(n) * slamch('Epsilon'). This error bound should only */
/* > be trusted if the previous boolean is true. */
/* > */
/* > err = 3 Reciprocal condition number: Estimated normwise */
/* > reciprocal condition number. Compared with the threshold */
/* > sqrt(n) * slamch('Epsilon') to determine if the error */
/* > estimate is "guaranteed". These reciprocal condition */
/* > numbers are 1 / (norm(Z^{
-1}
,inf) * norm(Z,inf)) for some */
/* > appropriately scaled matrix Z. */
/* > Let Z = S*A, where S scales each row by a power of the */
/* > radix so all absolute row sums of Z are approximately 1. */
/* > */
/* > This subroutine is only responsible for setting the second field */
/* > above. */
/* > See Lapack Working Note 165 for further details and extra */
/* > cautions. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ERRS_C */
/* > \verbatim */
/* > ERRS_C is REAL array, dimension (NRHS, N_ERR_BNDS) */
/* > For each right-hand side, this array contains information about */
/* > various error bounds and condition numbers corresponding to the */
/* > componentwise relative error, which is defined as follows: */
/* > */
/* > Componentwise relative error in the ith solution vector: */
/* > f2c_abs(XTRUE(j,i) - X(j,i)) */
/* > max_j ---------------------- */
/* > f2c_abs(X(j,i)) */
/* > */
/* > The array is indexed by the right-hand side i (on which the */
/* > componentwise relative error depends), and the type of error */
/* > information as described below. There currently are up to three */
/* > pieces of information returned for each right-hand side. If */
/* > componentwise accuracy is not requested (PARAMS(3) = 0.0), then */
/* > ERRS_C is not accessed. If N_ERR_BNDS .LT. 3, then at most */
/* > the first (:,N_ERR_BNDS) entries are returned. */
/* > */
/* > The first index in ERRS_C(i,:) corresponds to the ith */
/* > right-hand side. */
/* > */
/* > The second index in ERRS_C(:,err) contains the following */
/* > three fields: */
/* > err = 1 "Trust/don't trust" boolean. Trust the answer if the */
/* > reciprocal condition number is less than the threshold */
/* > sqrt(n) * slamch('Epsilon'). */
/* > */
/* > err = 2 "Guaranteed" error bound: The estimated forward error, */
/* > almost certainly within a factor of 10 of the true error */
/* > so long as the next entry is greater than the threshold */
/* > sqrt(n) * slamch('Epsilon'). This error bound should only */
/* > be trusted if the previous boolean is true. */
/* > */
/* > err = 3 Reciprocal condition number: Estimated componentwise */
/* > reciprocal condition number. Compared with the threshold */
/* > sqrt(n) * slamch('Epsilon') to determine if the error */
/* > estimate is "guaranteed". These reciprocal condition */
/* > numbers are 1 / (norm(Z^{
-1}
,inf) * norm(Z,inf)) for some */
/* > appropriately scaled matrix Z. */
/* > Let Z = S*(A*diag(x)), where x is the solution for the */
/* > current right-hand side and S scales each row of */
/* > A*diag(x) by a power of the radix so all absolute row */
/* > sums of Z are approximately 1. */
/* > */
/* > This subroutine is only responsible for setting the second field */
/* > above. */
/* > See Lapack Working Note 165 for further details and extra */
/* > cautions. */
/* > \endverbatim */
/* > */
/* > \param[in] RES */
/* > \verbatim */
/* > RES is COMPLEX array, dimension (N) */
/* > Workspace to hold the intermediate residual. */
/* > \endverbatim */
/* > */
/* > \param[in] AYB */
/* > \verbatim */
/* > AYB is REAL array, dimension (N) */
/* > Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] DY */
/* > \verbatim */
/* > DY is COMPLEX array, dimension (N) */
/* > Workspace to hold the intermediate solution. */
/* > \endverbatim */
/* > */
/* > \param[in] Y_TAIL */
/* > \verbatim */
/* > Y_TAIL is COMPLEX array, dimension (N) */
/* > Workspace to hold the trailing bits of the intermediate solution. */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > Reciprocal scaled condition number. This is an estimate of the */
/* > reciprocal Skeel condition number of the matrix A after */
/* > equilibration (if done). If this is less than the machine */
/* > precision (in particular, if it is zero), the matrix is singular */
/* > to working precision. Note that the error may still be small even */
/* > if this number is very small and the matrix appears ill- */
/* > conditioned. */
/* > \endverbatim */
/* > */
/* > \param[in] ITHRESH */
/* > \verbatim */
/* > ITHRESH is INTEGER */
/* > The maximum number of residual computations allowed for */
/* > refinement. The default is 10. For 'aggressive' set to 100 to */
/* > permit convergence using approximate factorizations or */
/* > factorizations other than LU. If the factorization uses a */
/* > technique other than Gaussian elimination, the guarantees in */
/* > ERRS_N and ERRS_C may no longer be trustworthy. */
/* > \endverbatim */
/* > */
/* > \param[in] RTHRESH */
/* > \verbatim */
/* > RTHRESH is REAL */
/* > Determines when to stop refinement if the error estimate stops */
/* > decreasing. Refinement will stop when the next solution no longer */
/* > satisfies norm(dx_{
i+1}
) < RTHRESH * norm(dx_i) where norm(Z) is */
/* > the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The */
/* > default value is 0.5. For 'aggressive' set to 0.9 to permit */
/* > convergence on extremely ill-conditioned matrices. See LAWN 165 */
/* > for more details. */
/* > \endverbatim */
/* > */
/* > \param[in] DZ_UB */
/* > \verbatim */
/* > DZ_UB is REAL */
/* > Determines when to start considering componentwise convergence. */
/* > Componentwise convergence is only considered after each component */
/* > of the solution Y is stable, which we definte as the relative */
/* > change in each component being less than DZ_UB. The default value */
/* > is 0.25, requiring the first bit to be stable. See LAWN 165 for */
/* > more details. */
/* > \endverbatim */
/* > */
/* > \param[in] IGNORE_CWISE */
/* > \verbatim */
/* > IGNORE_CWISE is LOGICAL */
/* > If .TRUE. then ignore componentwise convergence. Default value */
/* > is .FALSE.. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: Successful exit. */
/* > < 0: if INFO = -i, the ith argument to CGETRS had an illegal */
/* > value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexGEcomputational */
/* ===================================================================== */
/* Subroutine */
int cla_gerfsx_extended_(integer *prec_type__, integer * trans_type__, integer *n, integer *nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *ipiv, logical *colequ, real *c__, complex *b, integer *ldb, complex *y, integer *ldy, real *berr_out__, integer *n_norms__, real *errs_n__, real *errs_c__, complex *res, real *ayb, complex *dy, complex *y_tail__, real *rcond, integer * ithresh, real *rthresh, real *dz_ub__, logical *ignore_cwise__, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"cla_gerfsx_extended inputs: prec_type__ %lld, trans_type__ %lld, n %lld, nrhs %lld, lda %lld, ldaf %lld, ldb %lld, ldy %lld, n_norms__ %lld, ithresh %lld",*prec_type__, *trans_type__, *n, *nrhs, *lda, *ldaf, *ldb, *ldy, *n_norms__, *ithresh);
#else
    snprintf(buffer, 256,"cla_gerfsx_extended inputs: prec_type__ %d, trans_type__ %d, n %d, nrhs %d, lda %d, ldaf %d, ldb %d, ldy %d, n_norms__ %d, ithresh %d",*prec_type__, *trans_type__, *n, *nrhs, *lda, *ldaf, *ldb, *ldy, *n_norms__, *ithresh);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, y_dim1, y_offset, errs_n_dim1, errs_n_offset, errs_c_dim1, errs_c_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    char ch__1[1];
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    real dxratmax, dzratmax;
    integer i__, j;
    extern /* Subroutine */
    int cla_geamv_(integer *, integer *, integer *, real *, complex *, integer *, complex *, integer *, real *, real *, integer *);
    logical incr_prec__;
    real prev_dz_z__, yk, final_dx_x__;
    extern /* Subroutine */
    int cla_wwaddw_(integer *, complex *, complex *, complex *);
    real final_dz_z__, prevnormdx;
    integer cnt;
    real dyk, eps, incr_thresh__, dx_x__, dz_z__;
    extern /* Subroutine */
    int cla_lin_berr_(integer *, integer *, integer *, complex *, real *, real *);
    real ymin;
    extern /* Subroutine */
    int blas_cgemv_x_(integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *, integer *);
    integer y_prec_state__;
    extern /* Subroutine */
    int blas_cgemv2_x_(integer *, integer *, integer *, complex *, complex *, integer *, complex *, complex *, integer *, complex *, complex *, integer *, integer *), cgemv_(char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *), ccopy_( integer *, complex *, integer *, complex *, integer *);
    real dxrat, dzrat;
    extern /* Subroutine */
    int caxpy_(integer *, complex *, complex *, integer *, complex *, integer *);
    char trans[1];
    real normx, normy;
    extern real slamch_(char *);
    extern /* Subroutine */
    int cgetrs_(char *, integer *, integer *, complex *, integer *, integer *, complex *, integer *, integer *);
    real normdx;
    extern /* Character */
    VOID chla_transtype_(char *, integer *);
    real hugeval;
    integer x_state__, z_state__;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Parameters .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function Definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    errs_c_dim1 = *nrhs;
    errs_c_offset = 1 + errs_c_dim1;
    errs_c__ -= errs_c_offset;
    errs_n_dim1 = *nrhs;
    errs_n_offset = 1 + errs_n_dim1;
    errs_n__ -= errs_n_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    --ipiv;
    --c__;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --berr_out__;
    --res;
    --ayb;
    --dy;
    --y_tail__;
    /* Function Body */
    if (*info != 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    chla_transtype_(ch__1, trans_type__);
    *(unsigned char *)trans = *(unsigned char *)&ch__1[0];
    eps = slamch_("Epsilon");
    hugeval = slamch_("Overflow");
    /* Force HUGEVAL to Inf */
    hugeval *= hugeval;
    /* Using HUGEVAL may lead to spurious underflows. */
    incr_thresh__ = (real) (*n) * eps;
    i__1 = *nrhs;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        y_prec_state__ = 1;
        if (y_prec_state__ == 2)
        {
            i__2 = *n;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__;
                y_tail__[i__3].r = 0.f;
                y_tail__[i__3].i = 0.f; // , expr subst
            }
        }
        dxrat = 0.f;
        dxratmax = 0.f;
        dzrat = 0.f;
        dzratmax = 0.f;
        final_dx_x__ = hugeval;
        final_dz_z__ = hugeval;
        prevnormdx = hugeval;
        prev_dz_z__ = hugeval;
        dz_z__ = hugeval;
        dx_x__ = hugeval;
        x_state__ = 1;
        z_state__ = 0;
        incr_prec__ = FALSE_;
        i__2 = *ithresh;
        for (cnt = 1;
                cnt <= i__2;
                ++cnt)
        {
            /* Compute residual RES = B_s - op(A_s) * Y, */
            /* op(A) = A, A**T, or A**H depending on TRANS (and type). */
            ccopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
            if (y_prec_state__ == 0)
            {
                cgemv_(trans, n, n, &c_b6, &a[a_offset], lda, &y[j * y_dim1 + 1], &c__1, &c_b8, &res[1], &c__1);
            }
            else if (y_prec_state__ == 1)
            {
                blas_cgemv_x_(trans_type__, n, n, &c_b6, &a[a_offset], lda, & y[j * y_dim1 + 1], &c__1, &c_b8, &res[1], &c__1, prec_type__);
            }
            else
            {
                blas_cgemv2_x_(trans_type__, n, n, &c_b6, &a[a_offset], lda, &y[j * y_dim1 + 1], &y_tail__[1], &c__1, &c_b8, &res[ 1], &c__1, prec_type__);
            }
            /* XXX: RES is no longer needed. */
            ccopy_(n, &res[1], &c__1, &dy[1], &c__1);
            cgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &dy[1], n, info);
            /* Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT. */
            normx = 0.f;
            normy = 0.f;
            normdx = 0.f;
            dz_z__ = 0.f;
            ymin = hugeval;
            i__3 = *n;
            for (i__ = 1;
                    i__ <= i__3;
                    ++i__)
            {
                i__4 = i__ + j * y_dim1;
                yk = (r__1 = y[i__4].r, f2c_abs(r__1)) + (r__2 = r_imag(&y[i__ + j * y_dim1]), f2c_abs(r__2));
                i__4 = i__;
                dyk = (r__1 = dy[i__4].r, f2c_abs(r__1)) + (r__2 = r_imag(&dy[i__] ), f2c_abs(r__2));
                if (yk != 0.f)
                {
                    /* Computing MAX */
                    r__1 = dz_z__;
                    r__2 = dyk / yk; // , expr subst
                    dz_z__ = fla_max(r__1,r__2);
                }
                else if (dyk != 0.f)
                {
                    dz_z__ = hugeval;
                }
                ymin = fla_min(ymin,yk);
                normy = fla_max(normy,yk);
                if (*colequ)
                {
                    /* Computing MAX */
                    r__1 = normx;
                    r__2 = yk * c__[i__]; // , expr subst
                    normx = fla_max(r__1,r__2);
                    /* Computing MAX */
                    r__1 = normdx;
                    r__2 = dyk * c__[i__]; // , expr subst
                    normdx = fla_max(r__1,r__2);
                }
                else
                {
                    normx = normy;
                    normdx = fla_max(normdx,dyk);
                }
            }
            if (normx != 0.f)
            {
                dx_x__ = normdx / normx;
            }
            else if (normdx == 0.f)
            {
                dx_x__ = 0.f;
            }
            else
            {
                dx_x__ = hugeval;
            }
            dxrat = normdx / prevnormdx;
            dzrat = dz_z__ / prev_dz_z__;
            /* Check termination criteria */
            if (! (*ignore_cwise__) && ymin * *rcond < incr_thresh__ * normy && y_prec_state__ < 2)
            {
                incr_prec__ = TRUE_;
            }
            if (x_state__ == 3 && dxrat <= *rthresh)
            {
                x_state__ = 1;
            }
            if (x_state__ == 1)
            {
                if (dx_x__ <= eps)
                {
                    x_state__ = 2;
                }
                else if (dxrat > *rthresh)
                {
                    if (y_prec_state__ != 2)
                    {
                        incr_prec__ = TRUE_;
                    }
                    else
                    {
                        x_state__ = 3;
                    }
                }
                else
                {
                    if (dxrat > dxratmax)
                    {
                        dxratmax = dxrat;
                    }
                }
                if (x_state__ > 1)
                {
                    final_dx_x__ = dx_x__;
                }
            }
            if (z_state__ == 0 && dz_z__ <= *dz_ub__)
            {
                z_state__ = 1;
            }
            if (z_state__ == 3 && dzrat <= *rthresh)
            {
                z_state__ = 1;
            }
            if (z_state__ == 1)
            {
                if (dz_z__ <= eps)
                {
                    z_state__ = 2;
                }
                else if (dz_z__ > *dz_ub__)
                {
                    z_state__ = 0;
                    dzratmax = 0.f;
                    final_dz_z__ = hugeval;
                }
                else if (dzrat > *rthresh)
                {
                    if (y_prec_state__ != 2)
                    {
                        incr_prec__ = TRUE_;
                    }
                    else
                    {
                        z_state__ = 3;
                    }
                }
                else
                {
                    if (dzrat > dzratmax)
                    {
                        dzratmax = dzrat;
                    }
                }
                if (z_state__ > 1)
                {
                    final_dz_z__ = dz_z__;
                }
            }
            /* Exit if both normwise and componentwise stopped working, */
            /* but if componentwise is unstable, let it go at least two */
            /* iterations. */
            if (x_state__ != 1)
            {
                if (*ignore_cwise__)
                {
                    goto L666;
                }
                if (z_state__ == 3 || z_state__ == 2)
                {
                    goto L666;
                }
                if (z_state__ == 0 && cnt > 1)
                {
                    goto L666;
                }
            }
            if (incr_prec__)
            {
                incr_prec__ = FALSE_;
                ++y_prec_state__;
                i__3 = *n;
                for (i__ = 1;
                        i__ <= i__3;
                        ++i__)
                {
                    i__4 = i__;
                    y_tail__[i__4].r = 0.f;
                    y_tail__[i__4].i = 0.f; // , expr subst
                }
            }
            prevnormdx = normdx;
            prev_dz_z__ = dz_z__;
            /* Update soluton. */
            if (y_prec_state__ < 2)
            {
                caxpy_(n, &c_b8, &dy[1], &c__1, &y[j * y_dim1 + 1], &c__1);
            }
            else
            {
                cla_wwaddw_(n, &y[j * y_dim1 + 1], &y_tail__[1], &dy[1]);
            }
        }
        /* Target of "IF (Z_STOP .AND. X_STOP)". Sun's f77 won't CALL F90_EXIT. */
L666: /* Set final_* when cnt hits ithresh */
        if (x_state__ == 1)
        {
            final_dx_x__ = dx_x__;
        }
        if (z_state__ == 1)
        {
            final_dz_z__ = dz_z__;
        }
        /* Compute error bounds */
        if (*n_norms__ >= 1)
        {
            errs_n__[j + (errs_n_dim1 << 1)] = final_dx_x__ / (1 - dxratmax);
        }
        if (*n_norms__ >= 2)
        {
            errs_c__[j + (errs_c_dim1 << 1)] = final_dz_z__ / (1 - dzratmax);
        }
        /* Compute componentwise relative backward error from formula */
        /* fla_max(i) ( f2c_abs(R(i)) / ( f2c_abs(op(A_s))*f2c_abs(Y) + f2c_abs(B_s) )(i) ) */
        /* where f2c_abs(Z) is the componentwise absolute value of the matrix */
        /* or vector Z. */
        /* Compute residual RES = B_s - op(A_s) * Y, */
        /* op(A) = A, A**T, or A**H depending on TRANS (and type). */
        ccopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
        cgemv_(trans, n, n, &c_b6, &a[a_offset], lda, &y[j * y_dim1 + 1], & c__1, &c_b8, &res[1], &c__1);
        i__2 = *n;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            i__3 = i__ + j * b_dim1;
            ayb[i__] = (r__1 = b[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&b[i__ + j * b_dim1]), f2c_abs(r__2));
        }
        /* Compute f2c_abs(op(A_s))*f2c_abs(Y) + f2c_abs(B_s). */
        cla_geamv_(trans_type__, n, n, &c_b31, &a[a_offset], lda, &y[j * y_dim1 + 1], &c__1, &c_b31, &ayb[1], &c__1);
        cla_lin_berr_(n, n, &c__1, &res[1], &ayb[1], &berr_out__[j]);
        /* End of loop for each RHS. */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
}
/* cla_gerfsx_extended__ */
#endif
