/*
 *     Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.
 */

/* dbdsqr.f -- translated by f2c (version 20160102). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/*
 *  Modifications Copyright (c) 2024 Advanced Micro Devices, Inc.  All rights reserved.
 */
#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
#include "FLA_f2c.h" /* Table of constant values */
#include "fla_bdsqr_small_defs.h"

#if FLA_ENABLE_AMD_OPT
static doublereal c_b15 = -.125;
static aocl_int64_t c__1 = 1;
static doublereal c_b49 = 1.;
static doublereal c_b72 = -1.;
/* > \brief \b DBDSQR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DBDSQR + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbdsqr.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dbdsqr.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dbdsqr.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, */
/* LDU, C, LDC, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION C( LDC, * ), D( * ), E( * ), U( LDU, * ), */
/* $ VT( LDVT, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DBDSQR computes the singular values and, optionally, the right and/or */
/* > left singular vectors from the singular value decomposition (SVD) of */
/* > a real N-by-N (upper or lower) bidiagonal matrix B using the implicit */
/* > zero-shift QR algorithm. The SVD of B has the form */
/* > */
/* > B = Q * S * P**T */
/* > */
/* > where S is the diagonal matrix of singular values, Q is an orthogonal */
/* > matrix of left singular vectors, and P is an orthogonal matrix of */
/* > right singular vectors. If left singular vectors are requested, this */
/* > subroutine actually returns U*Q instead of Q, and, if right singular */
/* > vectors are requested, this subroutine returns P**T*VT instead of */
/* > P**T, for given real input matrices U and VT. When U and VT are the */
/* > orthogonal matrices that reduce a general matrix A to bidiagonal */
/* > form: A = U*B*VT, as computed by DGEBRD, then */
/* > */
/* > A = (U*Q) * S * (P**T*VT) */
/* > */
/* > is the SVD of A. Optionally, the subroutine may also compute Q**T*C */
/* > for a given real input matrix C. */
/* > */
/* > See "Computing Small Singular Values of Bidiagonal Matrices With */
/* > Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan, */
/* > LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11, */
/* > no. 5, pp. 873-912, Sept 1990) and */
/* > "Accurate singular values and differential qd algorithms," by */
/* > B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics */
/* > Department, University of California at Berkeley, July 1992 */
/* > for a detailed description of the algorithm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': B is upper bidiagonal;
 */
/* > = 'L': B is lower bidiagonal. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCVT */
/* > \verbatim */
/* > NCVT is INTEGER */
/* > The number of columns of the matrix VT. NCVT >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRU */
/* > \verbatim */
/* > NRU is INTEGER */
/* > The number of rows of the matrix U. NRU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCC */
/* > \verbatim */
/* > NCC is INTEGER */
/* > The number of columns of the matrix C. NCC >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > On entry, the n diagonal elements of the bidiagonal matrix B. */
/* > On exit, if INFO=0, the singular values of B in decreasing */
/* > order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (N-1) */
/* > On entry, the N-1 offdiagonal elements of the bidiagonal */
/* > matrix B. */
/* > On exit, if INFO = 0, E is destroyed;
if INFO > 0, D and E */
/* > will contain the diagonal and superdiagonal elements of a */
/* > bidiagonal matrix orthogonally equivalent to the one given */
/* > as input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VT */
/* > \verbatim */
/* > VT is DOUBLE PRECISION array, dimension (LDVT, NCVT) */
/* > On entry, an N-by-NCVT matrix VT. */
/* > On exit, VT is overwritten by P**T * VT. */
/* > Not referenced if NCVT = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* > LDVT is INTEGER */
/* > The leading dimension of the array VT. */
/* > LDVT >= fla_max(1,N) if NCVT > 0;
LDVT >= 1 if NCVT = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U */
/* > \verbatim */
/* > U is DOUBLE PRECISION array, dimension (LDU, N) */
/* > On entry, an NRU-by-N matrix U. */
/* > On exit, U is overwritten by U * Q. */
/* > Not referenced if NRU = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is INTEGER */
/* > The leading dimension of the array U. LDU >= fla_max(1,NRU). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (LDC, NCC) */
/* > On entry, an N-by-NCC matrix C. */
/* > On exit, C is overwritten by Q**T * C. */
/* > Not referenced if NCC = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. */
/* > LDC >= fla_max(1,N) if NCC > 0;
LDC >=1 if NCC = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (4*(N-1)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: If INFO = -i, the i-th argument had an illegal value */
/* > > 0: */
/* > if NCVT = NRU = NCC = 0, */
/* > = 1, a split was marked by a positive value in E */
/* > = 2, current block of Z not diagonalized after 30*N */
/* > iterations (in inner while loop) */
/* > = 3, termination criterion of outer while loop not met */
/* > (program created more than N unreduced blocks) */
/* > else NCVT = NRU = NCC = 0, */
/* > the algorithm did not converge;
D and E contain the */
/* > elements of a bidiagonal matrix which is orthogonally */
/* > similar to the input matrix B;
if INFO = i, i */
/* > elements of E have not converged to zero. */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > TOLMUL DOUBLE PRECISION, default = fla_max(10,fla_min(100,EPS**(-1/8))) */
/* > TOLMUL controls the convergence criterion of the QR loop. */
/* > If it is positive, TOLMUL*EPS is the desired relative */
/* > precision in the computed singular values. */
/* > If it is negative, f2c_abs(TOLMUL*EPS*sigma_max) is the */
/* > desired absolute accuracy in the computed singular */
/* > values (corresponds to relative accuracy */
/* > f2c_abs(TOLMUL*EPS) in the largest singular value. */
/* > f2c_abs(TOLMUL) should be between 1 and 1/EPS, and preferably */
/* > between 10 (for fast convergence) and .1/EPS */
/* > (for there to be some accuracy in the results). */
/* > Default is to lose at either one eighth or 2 of the */
/* > available decimal digits in each computed singular value */
/* > (whichever is smaller). */
/* > */
/* > MAXITR INTEGER, default = 6 */
/* > MAXITR controls the maximum number of passes of the */
/* > algorithm through its inner loop. The algorithms stops */
/* > (and so fails to converge) if the number of passes */
/* > through the inner loop exceeds MAXITR*N**2. */
/* > */
/* > \endverbatim */
/* > \par Note: */
/* =========== */
/* > */
/* > \verbatim */
/* > Bug report from Cezary Dendek. */
/* > On March 23rd 2017, the INTEGER variable MAXIT = MAXITR*N**2 is */
/* > removed since it can overflow pretty easily (for N larger or equal */
/* > than 18,919). We instead use MAXITDIVN = MAXITR*N. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup auxOTHERcomputational */
/* ===================================================================== */
#define BDSQR_PRE d
#include "fla_bdsqr_small_kernel.h"

#endif
