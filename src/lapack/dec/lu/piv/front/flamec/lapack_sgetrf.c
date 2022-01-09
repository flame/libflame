/*
    Copyright (c) 2021 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static real c_b11 = -1.f;
static real c_b12 = 1.f;

/* Subroutine */ integer lapack_sgetrf(integer *m, integer *n, real *a, integer *lda,
	integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    integer i__, j, jb, nb, iinfo;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
	    integer *, integer *);
    extern /* Subroutine */ int slaswp_(integer *, real *, integer *, integer
	    *, integer *, integer *, integer *);


/*  ======= */

/*  SGETRF computes an LU factorization of a general M-by-N matrix A */
/*  using partial pivoting with row interchanges. */

/*  The factorization has the form */
/*     A = P * L * U */
/*  where P is a permutation matrix, L is lower triangular with unit */
/*  diagonal elements (lower trapezoidal if m > n), and U is upper */
/*  triangular (upper trapezoidal if m < n). */

/*  This is the Crout Level 3 BLAS version of the algorithm. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  A       (input/output) REAL array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix to be factored. */
/*          On exit, the factors L and U from the factorization */
/*          A = P*L*U; the unit diagonal elements of L are not stored. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  IPIV    (output) INTEGER array, dimension (min(M,N)) */
/*          The pivot indices; for 1 <= i <= min(M,N), row i of the */
/*          matrix was interchanged with row IPIV(i). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization */
/*                has been completed, but the factor U is exactly */
/*                singular, and division by zero will occur if it is used */
/*                to solve a system of equations. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("LAPACK_SGETRF", &i__1);
	return *info;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return *info;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, "SGETRF", " ", m, n, &c_n1, &c_n1);
    if (nb <= 1 || nb >= min(*m,*n)) {

/*        Use unblocked code. */

	lapack_sgetf2(m, n, &a[a_offset], lda, &ipiv[1], info);
    } else {

/*        Use blocked code. */

	i__1 = min(*m,*n);
	i__2 = nb;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = min(*m,*n) - j + 1;
	    jb = min(i__3,nb);

/*           Update current block. */

	    i__3 = *m - j + 1;
	    i__4 = j - 1;
	    sgemm_("No transpose", "No transpose", &i__3, &jb, &i__4, &c_b11,
		    &a[j + a_dim1], lda, &a[j * a_dim1 + 1], lda, &c_b12, &a[
		    j + j * a_dim1], lda);

/*           Factor diagonal and subdiagonal blocks and test for exact */
/*           singularity. */

	    i__3 = *m - j + 1;
	    sgetrf2_(&i__3, &jb, &a[j + j * a_dim1], lda, &ipiv[j], &iinfo);

/*           Adjust INFO and the pivot indices. */

	    if (*info == 0 && iinfo > 0) {
		*info = iinfo + j - 1;
	    }
/* Computing MIN */
	    i__4 = *m, i__5 = j + jb - 1;
	    i__3 = min(i__4,i__5);
	    for (i__ = j; i__ <= i__3; ++i__) {
		ipiv[i__] = j - 1 + ipiv[i__];
/* L10: */
	    }

/*           Apply interchanges to column 1:J-1 */

	    i__3 = j - 1;
	    i__4 = j + jb - 1;
	    slaswp_(&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);

	    if (j + jb <= *n) {

/*              Apply interchanges to column J+JB:N */

		i__3 = *n - j - jb + 1;
		i__4 = j + jb - 1;
		slaswp_(&i__3, &a[(j + jb) * a_dim1 + 1], lda, &j, &i__4, &
			ipiv[1], &c__1);

		i__3 = *n - j - jb + 1;
		i__4 = j - 1;
		sgemm_("No transpose", "No transpose", &jb, &i__3, &i__4, &
			c_b11, &a[j + a_dim1], lda, &a[(j + jb) * a_dim1 + 1],
			 lda, &c_b12, &a[j + (j + jb) * a_dim1], lda);

/*              Compute block row of U. */

		i__3 = *n - j - jb + 1;
		strsm_("Left", "Lower", "No transpose", "Unit", &jb, &i__3, &
			c_b12, &a[j + j * a_dim1], lda, &a[j + (j + jb) *
			a_dim1], lda);
	    }
/* L20: */
	}
    }
    return *info;

/*     End of SGETRF */

} /* sgetrf_ */
