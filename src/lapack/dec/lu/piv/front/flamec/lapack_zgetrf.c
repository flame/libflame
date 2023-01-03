/*
    Copyright (c) 2021-2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"

/* Subroutine */ integer lapack_zgetrf(integer *m, integer *n, dcomplex *a,
	integer *lda, integer *ipiv, integer *info)
{



/*    Purpose
    =======

    ZGETRF computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.

    Arguments
    =========

    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= fla_max(1,M).

    IPIV    (output) INTEGER array, dimension (fla_min(M,N))
            The pivot indices; for 1 <= i <= fla_min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    =====================================================================


       Test the input parameters.

       Parameter adjustments */
    /* Table of constant values */
    static TLS_CLASS_SPEC dcomplex c_b1 = {1.,0.};
    static TLS_CLASS_SPEC integer c__1 = 1;
    static TLS_CLASS_SPEC integer c_n1 = -1;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    dcomplex z__1;
    /* Local variables */
    static TLS_CLASS_SPEC integer i__, j, iinfo;
    static TLS_CLASS_SPEC integer jb, nb;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,integer *, integer *);
#define a_subscr(a_1,a_2) (a_2)*a_dim1 + a_1
#define a_ref(a_1,a_2) a[a_subscr(a_1,a_2)]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --ipiv;
    #if AOCL_FLA_PROGRESS_H
        AOCL_FLA_PROGRESS_VAR;
    #endif

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < fla_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("LAPACK_ZGETRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, "ZGETRF", " ", m, n, &c_n1, &c_n1);
    if (nb <= 1 || nb >= fla_min(*m,*n)) {

/*        Use unblocked code. */
                #if AOCL_FLA_PROGRESS_H

                    #ifndef FLA_ENABLE_WINDOWS_BUILD
	    		if(!aocl_fla_progress_ptr)
                                aocl_fla_progress_ptr=aocl_fla_progress;
		    #endif

                        if(aocl_fla_progress_ptr){
                                step_count= fla_min(*m,*n);
                                AOCL_FLA_PROGRESS_FUNC_PTR("ZGETRF",6,&step_count,&thread_id,&total_threads);
                        }
                #endif

	lapack_zgetf2(m, n, &a[a_offset], lda, &ipiv[1], info);
    } else {

/*        Use blocked code. */

	#if AOCL_FLA_PROGRESS_H
                step_count =0;
        #endif
    
	i__1 = fla_min(*m,*n);
	i__2 = nb;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = fla_min(*m,*n) - j + 1;
	    jb = fla_min(i__3,nb);


            #if AOCL_FLA_PROGRESS_H
		#ifndef FLA_ENABLE_WINDOWS_BUILD
	    	    if(!aocl_fla_progress_ptr)
                        aocl_fla_progress_ptr=aocl_fla_progress;
		#endif

                    if(aocl_fla_progress_ptr){
                        step_count+=jb;
                        AOCL_FLA_PROGRESS_FUNC_PTR("ZGETRF",6,&step_count,&thread_id,&total_threads);
                    }
            #endif


/*           Factor diagonal and subdiagonal blocks and test for exact
             singularity. */

	    i__3 = *m - j + 1;
	    zgetrf2_(&i__3, &jb, &a_ref(j, j), lda, &ipiv[j], &iinfo);

/*           Adjust INFO and the pivot indices. */

	    if (*info == 0 && iinfo > 0) {
		*info = iinfo + j - 1;
	    }
/* Computing MIN */
	    i__4 = *m, i__5 = j + jb - 1;
	    i__3 = fla_min(i__4,i__5);
	    for (i__ = j; i__ <= i__3; ++i__) {
		ipiv[i__] = j - 1 + ipiv[i__];
/* L10: */
	    }

/*           Apply interchanges to columns 1:J-1. */

	    i__3 = j - 1;
	    i__4 = j + jb - 1;
	    zlaswp_(&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);

	    if (j + jb <= *n) {

/*              Apply interchanges to columns J+JB:N. */

		i__3 = *n - j - jb + 1;
		i__4 = j + jb - 1;
		zlaswp_(&i__3, &a_ref(1, j + jb), lda, &j, &i__4, &ipiv[1], &
			c__1);

/*              Compute block row of U. */

		i__3 = *n - j - jb + 1;
		ztrsm_("Left", "Lower", "No transpose", "Unit", &jb, &i__3, &
			c_b1, &a_ref(j, j), lda, &a_ref(j, j + jb), lda);
		if (j + jb <= *m) {

/*                 Update trailing submatrix. */

		    i__3 = *m - j - jb + 1;
		    i__4 = *n - j - jb + 1;
		    z__1.real = -1., z__1.imag = 0.;
		    zgemm_("No transpose", "No transpose", &i__3, &i__4, &jb,
			    &z__1, &a_ref(j + jb, j), lda, &a_ref(j, j + jb),
			    lda, &c_b1, &a_ref(j + jb, j + jb), lda);
		}
	    }
/* L20: */
	}
    }
    return *info;

/*     End of ZGETRF */

} /* zgetrf_ */

#undef a_ref
#undef a_subscr

