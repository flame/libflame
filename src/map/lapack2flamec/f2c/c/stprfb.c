/* ../netlib/stprfb.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b12 = 1.f;
static real c_b20 = 0.f;
static real c_b27 = -1.f;
/* > \brief \b STPRFB applies a real or complex "triangular-pentagonal" blocked reflector to a real or complex matrix, which is composed of two blocks. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STPRFB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stprfb. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stprfb. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stprfb. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STPRFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, */
/* V, LDV, T, LDT, A, LDA, B, LDB, WORK, LDWORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIRECT, SIDE, STOREV, TRANS */
/* INTEGER K, L, LDA, LDB, LDT, LDV, LDWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ), T( LDT, * ), */
/* $ V( LDV, * ), WORK( LDWORK, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPRFB applies a real "triangular-pentagonal" block reflector H or its */
/* > conjugate transpose H^H to a real matrix C, which is composed of two */
/* > blocks A and B, either from the left or right. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': apply H or H^H from the Left */
/* > = 'R': apply H or H^H from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': apply H (No transpose) */
/* > = 'C': apply H^H (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIRECT */
/* > \verbatim */
/* > DIRECT is CHARACTER*1 */
/* > Indicates how H is formed from a product of elementary */
/* > reflectors */
/* > = 'F': H = H(1) H(2) . . . H(k) (Forward) */
/* > = 'B': H = H(k) . . . H(2) H(1) (Backward) */
/* > \endverbatim */
/* > */
/* > \param[in] STOREV */
/* > \verbatim */
/* > STOREV is CHARACTER*1 */
/* > Indicates how the vectors which define the elementary */
/* > reflectors are stored: */
/* > = 'C': Columns */
/* > = 'R': Rows */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix B. */
/* > M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix B. */
/* > N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The order of the matrix T, i.e. the number of elementary */
/* > reflectors whose product defines the block reflector. */
/* > K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is INTEGER */
/* > The order of the trapezoidal part of V. */
/* > K >= L >= 0. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is REAL array, dimension */
/* > (LDV,K) if STOREV = 'C' */
/* > (LDV,M) if STOREV = 'R' and SIDE = 'L' */
/* > (LDV,N) if STOREV = 'R' and SIDE = 'R' */
/* > The pentagonal matrix V, which contains the elementary reflectors */
/* > H(1), H(2), ..., H(K). See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V. */
/* > If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*/
/* > if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*/
/* > if STOREV = 'R', LDV >= K. */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* > T is REAL array, dimension (LDT,K) */
/* > The triangular K-by-K matrix T in the representation of the */
/* > block reflector. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. */
/* > LDT >= K. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension */
/* > (LDA,N) if SIDE = 'L' or (LDA,K) if SIDE = 'R' */
/* > On entry, the K-by-N or M-by-K matrix A. */
/* > On exit, A is overwritten by the corresponding block of */
/* > H*C or H^H*C or C*H or C*H^H. See Futher Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. */
/* > If SIDE = 'L', LDC >= max(1,K);
*/
/* > If SIDE = 'R', LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,N) */
/* > On entry, the M-by-N matrix B. */
/* > On exit, B is overwritten by the corresponding block of */
/* > H*C or H^H*C or C*H or C*H^H. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. */
/* > LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension */
/* > (LDWORK,N) if SIDE = 'L', */
/* > (LDWORK,K) if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDWORK */
/* > \verbatim */
/* > LDWORK is INTEGER */
/* > The leading dimension of the array WORK. */
/* > If SIDE = 'L', LDWORK >= K;
*/
/* > if SIDE = 'R', LDWORK >= M. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix C is a composite matrix formed from blocks A and B. */
/* > The block B is of size M-by-N;
if SIDE = 'R', A is of size M-by-K, */
/* > and if SIDE = 'L', A is of size K-by-N. */
/* > */
/* > If SIDE = 'R' and DIRECT = 'F', C = [A B]. */
/* > */
/* > If SIDE = 'L' and DIRECT = 'F', C = [A] */
/* > [B]. */
/* > */
/* > If SIDE = 'R' and DIRECT = 'B', C = [B A]. */
/* > */
/* > If SIDE = 'L' and DIRECT = 'B', C = [B] */
/* > [A]. */
/* > */
/* > The pentagonal matrix V is composed of a rectangular block V1 and a */
/* > trapezoidal block V2. The size of the trapezoidal block is determined by */
/* > the parameter L, where 0<=L<=K. If L=K, the V2 block of V is triangular;
*/
/* > if L=0, there is no trapezoidal block, thus V = V1 is rectangular. */
/* > */
/* > If DIRECT = 'F' and STOREV = 'C': V = [V1] */
/* > [V2] */
/* > - V2 is upper trapezoidal (first L rows of K-by-K upper triangular) */
/* > */
/* > If DIRECT = 'F' and STOREV = 'R': V = [V1 V2] */
/* > */
/* > - V2 is lower trapezoidal (first L columns of K-by-K lower triangular) */
/* > */
/* > If DIRECT = 'B' and STOREV = 'C': V = [V2] */
/* > [V1] */
/* > - V2 is lower trapezoidal (last L rows of K-by-K lower triangular) */
/* > */
/* > If DIRECT = 'B' and STOREV = 'R': V = [V2 V1] */
/* > */
/* > - V2 is upper trapezoidal (last L columns of K-by-K upper triangular) */
/* > */
/* > If STOREV = 'C' and SIDE = 'L', V is M-by-K with V2 L-by-K. */
/* > */
/* > If STOREV = 'C' and SIDE = 'R', V is N-by-K with V2 L-by-K. */
/* > */
/* > If STOREV = 'R' and SIDE = 'L', V is K-by-M with V2 K-by-L. */
/* > */
/* > If STOREV = 'R' and SIDE = 'R', V is K-by-N with V2 K-by-L. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int stprfb_(char *side, char *trans, char *direct, char * storev, integer *m, integer *n, integer *k, integer *l, real *v, integer *ldv, real *t, integer *ldt, real *a, integer *lda, real *b, integer *ldb, real *work, integer *ldwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, work_offset, i__1, i__2;
    /* Local variables */
    logical backward;
    integer i__, j, kp, mp, np;
    logical row, left;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    logical right;
    extern /* Subroutine */
    int strmm_(char *, char *, char *, char *, integer *, integer *, real *, real *, integer *, real *, integer * );
    logical column, forward;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ========================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1;
    work -= work_offset;
    /* Function Body */
    if (*m <= 0 || *n <= 0 || *k <= 0 || *l < 0)
    {
        return 0;
    }
    if (lsame_(storev, "C"))
    {
        column = TRUE_;
        row = FALSE_;
    }
    else if (lsame_(storev, "R"))
    {
        column = FALSE_;
        row = TRUE_;
    }
    else
    {
        column = FALSE_;
        row = FALSE_;
    }
    if (lsame_(side, "L"))
    {
        left = TRUE_;
        right = FALSE_;
    }
    else if (lsame_(side, "R"))
    {
        left = FALSE_;
        right = TRUE_;
    }
    else
    {
        left = FALSE_;
        right = FALSE_;
    }
    if (lsame_(direct, "F"))
    {
        forward = TRUE_;
        backward = FALSE_;
    }
    else if (lsame_(direct, "B"))
    {
        forward = FALSE_;
        backward = TRUE_;
    }
    else
    {
        forward = FALSE_;
        backward = FALSE_;
    }
    /* --------------------------------------------------------------------------- */
    if (column && forward && left)
    {
        /* --------------------------------------------------------------------------- */
        /* Let W = [ I ] (K-by-K) */
        /* [ V ] (M-by-K) */
        /* Form H C or H^H C where C = [ A ] (K-by-N) */
        /* [ B ] (M-by-N) */
        /* H = I - W T W^H or H^H = I - W T^H W^H */
        /* A = A - T (A + V^H B) or A = A - T^H (A + V^H B) */
        /* B = B - V T (A + V^H B) or B = B - V T^H (A + V^H B) */
        /* --------------------------------------------------------------------------- */
        /* Computing MIN */
        i__1 = *m - *l + 1;
        mp = min(i__1,*m);
        /* Computing MIN */
        i__1 = *l + 1;
        kp = min(i__1,*k);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *l;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] = b[*m - *l + i__ + j * b_dim1];
            }
        }
        strmm_("L", "U", "T", "N", l, n, &c_b12, &v[mp + v_dim1], ldv, &work[ work_offset], ldwork);
        i__1 = *m - *l;
        sgemm_("T", "N", l, n, &i__1, &c_b12, &v[v_offset], ldv, &b[b_offset], ldb, &c_b12, &work[work_offset], ldwork);
        i__1 = *k - *l;
        sgemm_("T", "N", &i__1, n, m, &c_b12, &v[kp * v_dim1 + 1], ldv, &b[ b_offset], ldb, &c_b20, &work[kp + work_dim1], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *k;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
            }
        }
        strmm_("L", "U", trans, "N", k, n, &c_b12, &t[t_offset], ldt, &work[ work_offset], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *k;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
            }
        }
        i__1 = *m - *l;
        sgemm_("N", "N", &i__1, n, k, &c_b27, &v[v_offset], ldv, &work[ work_offset], ldwork, &c_b12, &b[b_offset], ldb);
        i__1 = *k - *l;
        sgemm_("N", "N", l, n, &i__1, &c_b27, &v[mp + kp * v_dim1], ldv, & work[kp + work_dim1], ldwork, &c_b12, &b[mp + b_dim1], ldb);
        strmm_("L", "U", "N", "N", l, n, &c_b12, &v[mp + v_dim1], ldv, &work[ work_offset], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *l;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                b[*m - *l + i__ + j * b_dim1] -= work[i__ + j * work_dim1];
            }
        }
        /* --------------------------------------------------------------------------- */
    }
    else if (column && forward && right)
    {
        /* --------------------------------------------------------------------------- */
        /* Let W = [ I ] (K-by-K) */
        /* [ V ] (N-by-K) */
        /* Form C H or C H^H where C = [ A B ] (A is M-by-K, B is M-by-N) */
        /* H = I - W T W^H or H^H = I - W T^H W^H */
        /* A = A - (A + B V) T or A = A - (A + B V) T^H */
        /* B = B - (A + B V) T V^H or B = B - (A + B V) T^H V^H */
        /* --------------------------------------------------------------------------- */
        /* Computing MIN */
        i__1 = *n - *l + 1;
        np = min(i__1,*n);
        /* Computing MIN */
        i__1 = *l + 1;
        kp = min(i__1,*k);
        i__1 = *l;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] = b[i__ + (*n - *l + j) * b_dim1];
            }
        }
        strmm_("R", "U", "N", "N", m, l, &c_b12, &v[np + v_dim1], ldv, &work[ work_offset], ldwork);
        i__1 = *n - *l;
        sgemm_("N", "N", m, l, &i__1, &c_b12, &b[b_offset], ldb, &v[v_offset], ldv, &c_b12, &work[work_offset], ldwork);
        i__1 = *k - *l;
        sgemm_("N", "N", m, &i__1, n, &c_b12, &b[b_offset], ldb, &v[kp * v_dim1 + 1], ldv, &c_b20, &work[kp * work_dim1 + 1], ldwork);
        i__1 = *k;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
            }
        }
        strmm_("R", "U", trans, "N", m, k, &c_b12, &t[t_offset], ldt, &work[ work_offset], ldwork);
        i__1 = *k;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
            }
        }
        i__1 = *n - *l;
        sgemm_("N", "T", m, &i__1, k, &c_b27, &work[work_offset], ldwork, &v[ v_offset], ldv, &c_b12, &b[b_offset], ldb);
        i__1 = *k - *l;
        sgemm_("N", "T", m, l, &i__1, &c_b27, &work[kp * work_dim1 + 1], ldwork, &v[np + kp * v_dim1], ldv, &c_b12, &b[np * b_dim1 + 1] , ldb);
        strmm_("R", "U", "T", "N", m, l, &c_b12, &v[np + v_dim1], ldv, &work[ work_offset], ldwork);
        i__1 = *l;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                b[i__ + (*n - *l + j) * b_dim1] -= work[i__ + j * work_dim1];
            }
        }
        /* --------------------------------------------------------------------------- */
    }
    else if (column && backward && left)
    {
        /* --------------------------------------------------------------------------- */
        /* Let W = [ V ] (M-by-K) */
        /* [ I ] (K-by-K) */
        /* Form H C or H^H C where C = [ B ] (M-by-N) */
        /* [ A ] (K-by-N) */
        /* H = I - W T W^H or H^H = I - W T^H W^H */
        /* A = A - T (A + V^H B) or A = A - T^H (A + V^H B) */
        /* B = B - V T (A + V^H B) or B = B - V T^H (A + V^H B) */
        /* --------------------------------------------------------------------------- */
        /* Computing MIN */
        i__1 = *l + 1;
        mp = min(i__1,*m);
        /* Computing MIN */
        i__1 = *k - *l + 1;
        kp = min(i__1,*k);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *l;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[*k - *l + i__ + j * work_dim1] = b[i__ + j * b_dim1];
            }
        }
        strmm_("L", "L", "T", "N", l, n, &c_b12, &v[kp * v_dim1 + 1], ldv, & work[kp + work_dim1], ldwork);
        i__1 = *m - *l;
        sgemm_("T", "N", l, n, &i__1, &c_b12, &v[mp + kp * v_dim1], ldv, &b[ mp + b_dim1], ldb, &c_b12, &work[kp + work_dim1], ldwork);
        i__1 = *k - *l;
        sgemm_("T", "N", &i__1, n, m, &c_b12, &v[v_offset], ldv, &b[b_offset], ldb, &c_b20, &work[work_offset], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *k;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
            }
        }
        strmm_("L", "L", trans, "N", k, n, &c_b12, &t[t_offset], ldt, &work[ work_offset], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *k;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
            }
        }
        i__1 = *m - *l;
        sgemm_("N", "N", &i__1, n, k, &c_b27, &v[mp + v_dim1], ldv, &work[ work_offset], ldwork, &c_b12, &b[mp + b_dim1], ldb);
        i__1 = *k - *l;
        sgemm_("N", "N", l, n, &i__1, &c_b27, &v[v_offset], ldv, &work[ work_offset], ldwork, &c_b12, &b[b_offset], ldb);
        strmm_("L", "L", "N", "N", l, n, &c_b12, &v[kp * v_dim1 + 1], ldv, & work[kp + work_dim1], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *l;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                b[i__ + j * b_dim1] -= work[*k - *l + i__ + j * work_dim1];
            }
        }
        /* --------------------------------------------------------------------------- */
    }
    else if (column && backward && right)
    {
        /* --------------------------------------------------------------------------- */
        /* Let W = [ V ] (N-by-K) */
        /* [ I ] (K-by-K) */
        /* Form C H or C H^H where C = [ B A ] (B is M-by-N, A is M-by-K) */
        /* H = I - W T W^H or H^H = I - W T^H W^H */
        /* A = A - (A + B V) T or A = A - (A + B V) T^H */
        /* B = B - (A + B V) T V^H or B = B - (A + B V) T^H V^H */
        /* --------------------------------------------------------------------------- */
        /* Computing MIN */
        i__1 = *l + 1;
        np = min(i__1,*n);
        /* Computing MIN */
        i__1 = *k - *l + 1;
        kp = min(i__1,*k);
        i__1 = *l;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + (*k - *l + j) * work_dim1] = b[i__ + j * b_dim1];
            }
        }
        strmm_("R", "L", "N", "N", m, l, &c_b12, &v[kp * v_dim1 + 1], ldv, & work[kp * work_dim1 + 1], ldwork);
        i__1 = *n - *l;
        sgemm_("N", "N", m, l, &i__1, &c_b12, &b[np * b_dim1 + 1], ldb, &v[np + kp * v_dim1], ldv, &c_b12, &work[kp * work_dim1 + 1], ldwork);
        i__1 = *k - *l;
        sgemm_("N", "N", m, &i__1, n, &c_b12, &b[b_offset], ldb, &v[v_offset], ldv, &c_b20, &work[work_offset], ldwork);
        i__1 = *k;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
            }
        }
        strmm_("R", "L", trans, "N", m, k, &c_b12, &t[t_offset], ldt, &work[ work_offset], ldwork);
        i__1 = *k;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
            }
        }
        i__1 = *n - *l;
        sgemm_("N", "T", m, &i__1, k, &c_b27, &work[work_offset], ldwork, &v[ np + v_dim1], ldv, &c_b12, &b[np * b_dim1 + 1], ldb);
        i__1 = *k - *l;
        sgemm_("N", "T", m, l, &i__1, &c_b27, &work[work_offset], ldwork, &v[ v_offset], ldv, &c_b12, &b[b_offset], ldb);
        strmm_("R", "L", "T", "N", m, l, &c_b12, &v[kp * v_dim1 + 1], ldv, & work[kp * work_dim1 + 1], ldwork);
        i__1 = *l;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                b[i__ + j * b_dim1] -= work[i__ + (*k - *l + j) * work_dim1];
            }
        }
        /* --------------------------------------------------------------------------- */
    }
    else if (row && forward && left)
    {
        /* --------------------------------------------------------------------------- */
        /* Let W = [ I V ] ( I is K-by-K, V is K-by-M ) */
        /* Form H C or H^H C where C = [ A ] (K-by-N) */
        /* [ B ] (M-by-N) */
        /* H = I - W^H T W or H^H = I - W^H T^H W */
        /* A = A - T (A + V B) or A = A - T^H (A + V B) */
        /* B = B - V^H T (A + V B) or B = B - V^H T^H (A + V B) */
        /* --------------------------------------------------------------------------- */
        /* Computing MIN */
        i__1 = *m - *l + 1;
        mp = min(i__1,*m);
        /* Computing MIN */
        i__1 = *l + 1;
        kp = min(i__1,*k);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *l;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] = b[*m - *l + i__ + j * b_dim1];
            }
        }
        strmm_("L", "L", "N", "N", l, n, &c_b12, &v[mp * v_dim1 + 1], ldv, & work[work_offset], ldb);
        i__1 = *m - *l;
        sgemm_("N", "N", l, n, &i__1, &c_b12, &v[v_offset], ldv, &b[b_offset], ldb, &c_b12, &work[work_offset], ldwork);
        i__1 = *k - *l;
        sgemm_("N", "N", &i__1, n, m, &c_b12, &v[kp + v_dim1], ldv, &b[ b_offset], ldb, &c_b20, &work[kp + work_dim1], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *k;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
            }
        }
        strmm_("L", "U", trans, "N", k, n, &c_b12, &t[t_offset], ldt, &work[ work_offset], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *k;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
            }
        }
        i__1 = *m - *l;
        sgemm_("T", "N", &i__1, n, k, &c_b27, &v[v_offset], ldv, &work[ work_offset], ldwork, &c_b12, &b[b_offset], ldb);
        i__1 = *k - *l;
        sgemm_("T", "N", l, n, &i__1, &c_b27, &v[kp + mp * v_dim1], ldv, & work[kp + work_dim1], ldwork, &c_b12, &b[mp + b_dim1], ldb);
        strmm_("L", "L", "T", "N", l, n, &c_b12, &v[mp * v_dim1 + 1], ldv, & work[work_offset], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *l;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                b[*m - *l + i__ + j * b_dim1] -= work[i__ + j * work_dim1];
            }
        }
        /* --------------------------------------------------------------------------- */
    }
    else if (row && forward && right)
    {
        /* --------------------------------------------------------------------------- */
        /* Let W = [ I V ] ( I is K-by-K, V is K-by-N ) */
        /* Form C H or C H^H where C = [ A B ] (A is M-by-K, B is M-by-N) */
        /* H = I - W^H T W or H^H = I - W^H T^H W */
        /* A = A - (A + B V^H) T or A = A - (A + B V^H) T^H */
        /* B = B - (A + B V^H) T V or B = B - (A + B V^H) T^H V */
        /* --------------------------------------------------------------------------- */
        /* Computing MIN */
        i__1 = *n - *l + 1;
        np = min(i__1,*n);
        /* Computing MIN */
        i__1 = *l + 1;
        kp = min(i__1,*k);
        i__1 = *l;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] = b[i__ + (*n - *l + j) * b_dim1];
            }
        }
        strmm_("R", "L", "T", "N", m, l, &c_b12, &v[np * v_dim1 + 1], ldv, & work[work_offset], ldwork);
        i__1 = *n - *l;
        sgemm_("N", "T", m, l, &i__1, &c_b12, &b[b_offset], ldb, &v[v_offset], ldv, &c_b12, &work[work_offset], ldwork);
        i__1 = *k - *l;
        sgemm_("N", "T", m, &i__1, n, &c_b12, &b[b_offset], ldb, &v[kp + v_dim1], ldv, &c_b20, &work[kp * work_dim1 + 1], ldwork);
        i__1 = *k;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
            }
        }
        strmm_("R", "U", trans, "N", m, k, &c_b12, &t[t_offset], ldt, &work[ work_offset], ldwork);
        i__1 = *k;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
            }
        }
        i__1 = *n - *l;
        sgemm_("N", "N", m, &i__1, k, &c_b27, &work[work_offset], ldwork, &v[ v_offset], ldv, &c_b12, &b[b_offset], ldb);
        i__1 = *k - *l;
        sgemm_("N", "N", m, l, &i__1, &c_b27, &work[kp * work_dim1 + 1], ldwork, &v[kp + np * v_dim1], ldv, &c_b12, &b[np * b_dim1 + 1] , ldb);
        strmm_("R", "L", "N", "N", m, l, &c_b12, &v[np * v_dim1 + 1], ldv, & work[work_offset], ldwork);
        i__1 = *l;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                b[i__ + (*n - *l + j) * b_dim1] -= work[i__ + j * work_dim1];
            }
        }
        /* --------------------------------------------------------------------------- */
    }
    else if (row && backward && left)
    {
        /* --------------------------------------------------------------------------- */
        /* Let W = [ V I ] ( I is K-by-K, V is K-by-M ) */
        /* Form H C or H^H C where C = [ B ] (M-by-N) */
        /* [ A ] (K-by-N) */
        /* H = I - W^H T W or H^H = I - W^H T^H W */
        /* A = A - T (A + V B) or A = A - T^H (A + V B) */
        /* B = B - V^H T (A + V B) or B = B - V^H T^H (A + V B) */
        /* --------------------------------------------------------------------------- */
        /* Computing MIN */
        i__1 = *l + 1;
        mp = min(i__1,*m);
        /* Computing MIN */
        i__1 = *k - *l + 1;
        kp = min(i__1,*k);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *l;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[*k - *l + i__ + j * work_dim1] = b[i__ + j * b_dim1];
            }
        }
        strmm_("L", "U", "N", "N", l, n, &c_b12, &v[kp + v_dim1], ldv, &work[ kp + work_dim1], ldwork);
        i__1 = *m - *l;
        sgemm_("N", "N", l, n, &i__1, &c_b12, &v[kp + mp * v_dim1], ldv, &b[ mp + b_dim1], ldb, &c_b12, &work[kp + work_dim1], ldwork);
        i__1 = *k - *l;
        sgemm_("N", "N", &i__1, n, m, &c_b12, &v[v_offset], ldv, &b[b_offset], ldb, &c_b20, &work[work_offset], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *k;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
            }
        }
        strmm_("L", "L ", trans, "N", k, n, &c_b12, &t[t_offset], ldt, &work[ work_offset], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *k;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
            }
        }
        i__1 = *m - *l;
        sgemm_("T", "N", &i__1, n, k, &c_b27, &v[mp * v_dim1 + 1], ldv, &work[ work_offset], ldwork, &c_b12, &b[mp + b_dim1], ldb);
        i__1 = *k - *l;
        sgemm_("T", "N", l, n, &i__1, &c_b27, &v[v_offset], ldv, &work[ work_offset], ldwork, &c_b12, &b[b_offset], ldb);
        strmm_("L", "U", "T", "N", l, n, &c_b12, &v[kp + v_dim1], ldv, &work[ kp + work_dim1], ldwork);
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *l;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                b[i__ + j * b_dim1] -= work[*k - *l + i__ + j * work_dim1];
            }
        }
        /* --------------------------------------------------------------------------- */
    }
    else if (row && backward && right)
    {
        /* --------------------------------------------------------------------------- */
        /* Let W = [ V I ] ( I is K-by-K, V is K-by-N ) */
        /* Form C H or C H^H where C = [ B A ] (A is M-by-K, B is M-by-N) */
        /* H = I - W^H T W or H^H = I - W^H T^H W */
        /* A = A - (A + B V^H) T or A = A - (A + B V^H) T^H */
        /* B = B - (A + B V^H) T V or B = B - (A + B V^H) T^H V */
        /* --------------------------------------------------------------------------- */
        /* Computing MIN */
        i__1 = *l + 1;
        np = min(i__1,*n);
        /* Computing MIN */
        i__1 = *k - *l + 1;
        kp = min(i__1,*k);
        i__1 = *l;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + (*k - *l + j) * work_dim1] = b[i__ + j * b_dim1];
            }
        }
        strmm_("R", "U", "T", "N", m, l, &c_b12, &v[kp + v_dim1], ldv, &work[ kp * work_dim1 + 1], ldwork);
        i__1 = *n - *l;
        sgemm_("N", "T", m, l, &i__1, &c_b12, &b[np * b_dim1 + 1], ldb, &v[kp + np * v_dim1], ldv, &c_b12, &work[kp * work_dim1 + 1], ldwork);
        i__1 = *k - *l;
        sgemm_("N", "T", m, &i__1, n, &c_b12, &b[b_offset], ldb, &v[v_offset], ldv, &c_b20, &work[work_offset], ldwork);
        i__1 = *k;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
            }
        }
        strmm_("R", "L", trans, "N", m, k, &c_b12, &t[t_offset], ldt, &work[ work_offset], ldwork);
        i__1 = *k;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
            }
        }
        i__1 = *n - *l;
        sgemm_("N", "N", m, &i__1, k, &c_b27, &work[work_offset], ldwork, &v[ np * v_dim1 + 1], ldv, &c_b12, &b[np * b_dim1 + 1], ldb);
        i__1 = *k - *l;
        sgemm_("N", "N", m, l, &i__1, &c_b27, &work[work_offset], ldwork, &v[ v_offset], ldv, &c_b12, &b[b_offset], ldb);
        strmm_("R", "U", "N", "N", m, l, &c_b12, &v[kp + v_dim1], ldv, &work[ kp * work_dim1 + 1], ldwork);
        i__1 = *l;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                b[i__ + j * b_dim1] -= work[i__ + (*k - *l + j) * work_dim1];
            }
        }
    }
    return 0;
    /* End of STPRFB */
}
/* stprfb_ */
