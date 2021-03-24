/*
    Copyright (c) 2021 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"

/* Subroutine */ int lapack_sgetf2(integer *m, integer *n, real *a, integer *lda,
	integer *ipiv, integer *info)
{


/*    Purpose
    =======

    SGETF2 computes an LU factorization of a general m-by-n matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 2 BLAS version of the algorithm.

    Arguments
    =========

    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) REAL array, dimension (LDA,N)
            On entry, the m by n matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    IPIV    (output) INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0: successful exit
            < 0: if INFO = -k, the k-th argument had an illegal value
            > 0: if INFO = k, U(k,k) is exactly zero. The factorization
                 has been completed, but the factor U is exactly
                 singular, and division by zero will occur if it is used
                 to solve a system of equations.

    =====================================================================


       Test the input parameters.

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static real c_b6 = -1.f;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;
    /* Local variables */
	static integer j;
    static integer jp;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer isamax_(integer *, real *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
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
	xerbla_("LAPACK_SGETF2", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = min(*m,*n);
    for (j = 1; j <= i__1; ++j) {

/*        Find pivot and test for singularity. */

	i__2 = *m - j + 1;
	jp = j - 1 + isamax_(&i__2, &a_ref(j, j), &c__1);
	ipiv[j] = jp;
	if (a_ref(jp, j) != 0.f) {

/*           Apply the interchange to columns 1:N. */

	    if (jp != j) {
		sswap_(n, &a_ref(j, 1), lda, &a_ref(jp, 1), lda);
	    }

/*           Compute elements J+1:M of J-th column. */

	    if (j < *m) {
		i__2 = *m - j;
		r__1 = 1.f / a_ref(j, j);
		sscal_(&i__2, &r__1, &a_ref(j + 1, j), &c__1);
	    }

	} else if (*info == 0) {

	    *info = j;
	}

	if (j < min(*m,*n)) {

/*           Update trailing submatrix. */

	    i__2 = *m - j;
	    i__3 = *n - j;
	    sger_(&i__2, &i__3, &c_b6, &a_ref(j + 1, j), &c__1, &a_ref(j, j +
		    1), lda, &a_ref(j + 1, j + 1), lda);
	}
/* L10: */
    }
    return 0;

/*     End of SGETF2 */

} /* sgetf2_ */

#undef a_ref

