/*
    Copyright (c) 2021-2022 Advanced Micro Devices, Inc. All rights reserved.
*/

#include "FLAME.h"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static real c_b13 = -1.f;
static real c_b14 = 1.f;

/* Subroutine */ int lapack_spotrf(char *uplo, integer *n, real *a, integer *lda,
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer j, jb, nb;
    logical upper;

	int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
	logical lsame_(char *ca, char *cb);
	int lapack_spotf2(char *uplo, integer *n, real *a, integer *lda, integer *info);

/*  SPOTRF computes the Cholesky factorization of a real symmetric */
/*  positive definite matrix A. */

/*  The factorization has the form */
/*     A = U**T * U,  if UPLO = 'U', or */
/*     A = L  * L**T,  if UPLO = 'L', */
/*  where U is an upper triangular matrix and L is lower triangular. */

/*  This is the block version of the algorithm, calling Level 3 BLAS. */

/*  Arguments */
/*  ========= */

/*  UPLO    (input) CHARACTER*1 */
/*          = 'U':  Upper triangle of A is stored; */
/*          = 'L':  Lower triangle of A is stored. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  A       (input/output) REAL array, dimension (LDA,N) */
/*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/*          N-by-N upper triangular part of A contains the upper */
/*          triangular part of the matrix A, and the strictly lower */
/*          triangular part of A is not referenced.  If UPLO = 'L', the */
/*          leading N-by-N lower triangular part of A contains the lower */
/*          triangular part of the matrix A, and the strictly upper */
/*          triangular part of A is not referenced. */

/*          On exit, if INFO = 0, the factor U or L from the Cholesky */
/*          factorization A = U**T*U or A = L*L**T. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= fla_max(1,N). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, the leading minor of order i is not */
/*                positive definite, and the factorization could not be */
/*                completed. */

/*  ===================================================================== */
 /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    #if AOCL_FLA_PROGRESS_H
        AOCL_FLA_PROGRESS_VAR;
	progress_step_count=0;
      #ifndef FLA_ENABLE_WINDOWS_BUILD
	if(!aocl_fla_progress_ptr)
            aocl_fla_progress_ptr=aocl_fla_progress;
      #endif 

    #endif
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < fla_max(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("LAPACK_SPOTRF", &i__1, (ftnlen)13);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, "SPOTRF", uplo, n, &c_n1, &c_n1, &c_n1);
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code. */

	lapack_spotf2(uplo, n, &a[a_offset], lda, info);
    } else {

/*        Use blocked code. */
	if (upper) {

/*           Compute the Cholesky factorization A = U'*U. */

	    i__1 = *n;
	    i__2 = nb;
	    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

/* Computing MIN */
		i__3 = nb, i__4 = *n - j + 1;
		jb = fla_min(i__3,i__4);
		i__3 = j - 1;
		#if AOCL_FLA_PROGRESS_H
		    if(aocl_fla_progress_ptr){
                	progress_step_count+=jb;
                	AOCL_FLA_PROGRESS_FUNC_PTR("SPOTRF",6,&progress_step_count,&progress_thread_id,&progress_total_threads);
            	    }
        	#endif 
		ssyrk_("Upper", "Transpose", &jb, &i__3, &c_b13, &a[j *
			a_dim1 + 1], lda, &c_b14, &a[j + j * a_dim1], lda);
		lapack_spotf2("Upper", &jb, &a[j + j * a_dim1], lda, info);
		if (*info != 0) {
		    goto L30;
		}
		if (j + jb <= *n) {

/*                 Compute the current block row. */

		    i__3 = *n - j - jb + 1;
		    i__4 = j - 1;
		    sgemm_("Transpose", "No transpose", &jb, &i__3, &i__4, &
			    c_b13, &a[j * a_dim1 + 1], lda, &a[(j + jb) *
			    a_dim1 + 1], lda, &c_b14, &a[j + (j + jb) *
			    a_dim1], lda);
		    i__3 = *n - j - jb + 1;
		    strsm_("Left", "Upper", "Transpose", "Non-unit", &jb, &
			    i__3, &c_b14, &a[j + j * a_dim1], lda, &a[j + (j
			    + jb) * a_dim1], lda);
		}
/* L10: */
	    }

	} else {

/*           Compute the Cholesky factorization A = L*L'. */

	    i__2 = *n;
	    i__1 = nb;
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

/* Computing MIN */
		i__3 = nb, i__4 = *n - j + 1;
		jb = fla_min(i__3,i__4);
		i__3 = j - 1;
		#if AOCL_FLA_PROGRESS_H
		    if(aocl_fla_progress_ptr){
                	progress_step_count+=jb;
                	AOCL_FLA_PROGRESS_FUNC_PTR("SPOTRF",6,&progress_step_count,&progress_thread_id,&progress_total_threads);
                    }
                #endif
		ssyrk_("Lower", "No transpose", &jb, &i__3, &c_b13, &a[j +
			a_dim1], lda, &c_b14, &a[j + j * a_dim1], lda);
		lapack_spotf2("Lower", &jb, &a[j + j * a_dim1], lda, info);
		if (*info != 0) {
		    goto L30;
		}
		if (j + jb <= *n) {

/*                 Compute the current block column. */

		    i__3 = *n - j - jb + 1;
		    i__4 = j - 1;
		    sgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &
			    c_b13, &a[j + jb + a_dim1], lda, &a[j + a_dim1],
			    lda, &c_b14, &a[j + jb + j * a_dim1], lda);
		    i__3 = *n - j - jb + 1;
		    strsm_("Right", "Lower", "Transpose", "Non-unit", &i__3, &
			    jb, &c_b14, &a[j + j * a_dim1], lda, &a[j + jb +
			    j * a_dim1], lda);
		}
/* L20: */
	    }
	}
    }
    goto L40;

L30:
    *info = *info + j - 1;

L40:
    return 0;

/*     End of SPOTRF */

} /* spotrf_ */
