/* ../netlib/dlacpy.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/******************************************************************************
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#include "FLA_f2c.h" /* > \brief \b DLACPY copies all or part of one two-dimensional array to another. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLACPY + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacpy.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacpy.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacpy.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER LDA, LDB, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLACPY copies all or part of a two-dimensional matrix A to another */
/* > matrix B. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies the part of the matrix A to be copied to B. */
/* > = 'U': Upper triangular part */
/* > = 'L': Lower triangular part */
/* > Otherwise: All of the matrix A */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > The m by n matrix A. If UPLO = 'U', only the upper triangle */
/* > or trapezoid is accessed;
if UPLO = 'L', only the lower */
/* > triangle or trapezoid is accessed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,N) */
/* > On exit, B = A in the locations specified by UPLO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,M). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void dlacpy_(char *uplo, aocl_int_t *m, aocl_int_t *n, doublereal *a, aocl_int_t *lda,
             doublereal *b, aocl_int_t *ldb)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dlacpy(uplo, m, n, a, lda, b, ldb);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;

    aocl_lapack_dlacpy(uplo, &m_64, &n_64, a, &lda_64, b, &ldb_64);
#endif
}

void aocl_lapack_dlacpy(char *uplo, aocl_int64_t *m, aocl_int64_t *n, doublereal *a,
                        aocl_int64_t *lda, doublereal *b, aocl_int64_t *ldb)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlacpy inputs: uplo %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *uplo, *m, *n, *lda, *ldb);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
#ifdef FLA_ENABLE_AMD_OPT
    /* Handling invalid cases */
    if(*m <= 0 || *n <= 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }

    if(lsame_(uplo, "U", 1, 1))
    {
        /* Upper triangular part */
        for(aocl_int64_t i = 1; i <= *n; i++)
        {
            memcpy(b + (i * b_dim1) + 1, a + (i * a_dim1) + 1, fla_min(i, *m) * sizeof(double));
        }
    }
    else if(lsame_(uplo, "L", 1, 1))
    {
        /* Lower triangular part */
        for(aocl_int64_t i = 1; i <= *n && i <= *m; i++)
        {
            memcpy(b + (i * b_dim1 + i), a + (i * a_dim1 + i), (*m - i + 1) * sizeof(double));
        }
    }
    else
    {
        /* All of the matrix A */
        /* If the leading dimension of A and B are the same, and the number of rows of A and B are
         * the same, then we can use memcpy to copy the entire matrix. */
        if(a_dim1 == b_dim1 && a_dim1 == *m)
        {
            memcpy(b + b_offset, a + a_offset, a_dim1 * *n * sizeof(double));
        }
        else
        {
            /* Otherwise, we need to copy the matrix column by column. */
            for(aocl_int64_t i = 1; i <= *n; i++)
            {
                memcpy(b + (i * b_dim1) + 1, a + (i * a_dim1) + 1, *m * sizeof(double));
            }
        }
    }
#else
    aocl_int64_t i__1, i__2, i__, j;
    if(lsame_(uplo, "U", 1, 1))
    {
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = fla_min(j, *m);
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
                /* L10: */
            }
            /* L20: */
        }
    }
    else if(lsame_(uplo, "L", 1, 1))
    {
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *m;
            for(i__ = j; i__ <= i__2; ++i__)
            {
                b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
                /* L30: */
            }
            /* L40: */
        }
    }
    else
    {
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *m;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
                /* L50: */
            }
            /* L60: */
        }
    }
#endif
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLACPY */
}
/* dlacpy_ */
