/*
    Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
*/
/* ../netlib/dgtsv.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief <b> DGTSV computes the solution to system of linear equations A * X = B for GT matrices <b> */
#if FLA_ENABLE_AMD_OPT
/* Forward-declare only the kernels used below. The AVX2/AVX512 kernel headers
   are intentionally not included here: they unconditionally pull in
   x86-specific headers (e.g. immintrin.h), which would break non-AMD-opt
   (and potentially non-x86) builds. */
int fla_dgtsv_kernel_avx2(aocl_int64_t *n, aocl_int64_t *nrhs, doublereal *dl, doublereal *d__,
                          doublereal *du, doublereal *b, aocl_int64_t *ldb, aocl_int64_t *info);
int fla_dgtsv_kernel_avx512(aocl_int64_t *n, aocl_int64_t *nrhs, doublereal *dl, doublereal *d__,
                            doublereal *du, doublereal *b, aocl_int64_t *ldb, aocl_int64_t *info);
#endif
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGTSV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtsv.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtsv.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtsv.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION B( LDB, * ), D( * ), DL( * ), DU( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGTSV solves the equation */
/* > */
/* > A*X = B, */
/* > */
/* > where A is an n by n tridiagonal matrix, by Gaussian elimination with */
/* > partial pivoting. */
/* > */
/* > Note that the equation A**T*X = B may be solved by interchanging the */
/* > order of the arguments DU and DL. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DL */
/* > \verbatim */
/* > DL is DOUBLE PRECISION array, dimension (N-1) */
/* > On entry, DL must contain the (n-1) sub-diagonal elements of */
/* > A. */
/* > */
/* > On exit, DL is overwritten by the (n-2) elements of the */
/* > second super-diagonal of the upper triangular matrix U from */
/* > the LU factorization of A, in DL(1), ..., DL(n-2). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > On entry, D must contain the diagonal elements of A. */
/* > */
/* > On exit, D is overwritten by the n diagonal elements of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* > DU is DOUBLE PRECISION array, dimension (N-1) */
/* > On entry, DU must contain the (n-1) super-diagonal elements */
/* > of A. */
/* > */
/* > On exit, DU is overwritten by the (n-1) elements of the first */
/* > super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* > On entry, the N by NRHS matrix of right hand side matrix B. */
/* > On exit, if INFO = 0, the N by NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, U(i,i) is exactly zero, and the solution */
/* > has not been computed. The factorization has not been */
/* > completed unless i = N. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleGTsolve */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void dgtsv_(aocl_int_t *n, aocl_int_t *nrhs, doublereal *dl, doublereal *d__, doublereal *du,
            doublereal *b, aocl_int_t *ldb, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgtsv(n, nrhs, dl, d__, du, b, ldb, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgtsv(&n_64, &nrhs_64, dl, d__, du, b, &ldb_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_dgtsv(aocl_int64_t *n, aocl_int64_t *nrhs, doublereal *dl, doublereal *d__,
                       doublereal *du, doublereal *b, aocl_int64_t *ldb, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgtsv inputs: n %" FLA_IS ", nrhs %" FLA_IS ", ldb %" FLA_IS "", *n, *nrhs,
                      *ldb);
    /* System generated locals */
    aocl_int64_t b_dim1, b_offset, i__1, i__2;
    doublereal d__1, d__2;
    /* Local variables */
    aocl_int64_t i__, j;
    doublereal fact, temp;
    /* -- LAPACK driver routine (version 3.4.2) -- */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Function Body */
    *info = 0;
    if(*n < 0)
    {
        *info = -1;
    }
    else if(*nrhs < 0)
    {
        *info = -2;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -7;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DGTSV ", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*nrhs == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#if FLA_ENABLE_AMD_OPT
    /* Initialize the global context (architecture detection) and dispatch to the
       AVX512/AVX2 optimized kernels based on problem size. The kernels vectorize
       across the NRHS columns, so the selection is driven by *nrhs:
         - AVX512 (8 columns/iteration) once *nrhs is large enough to fill the
           512-bit registers,
         - AVX2 (4 columns/iteration) for the mid range,
         - the scalar reference path below for very small *nrhs, where SIMD gather
           / scatter overheads outweigh any benefit. */
    aocl_fla_init();
    if(*nrhs >= FLA_DGTSV_AVX512_MIN_NRHS && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_dgtsv_kernel_avx512(n, nrhs, dl, d__, du, b, ldb, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*nrhs >= FLA_DGTSV_AVX2_MIN_NRHS && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgtsv_kernel_avx2(n, nrhs, dl, d__, du, b, ldb, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    if(*nrhs == 1)
    {
        i__1 = *n - 2;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            if((d__1 = d__[i__], f2c_dabs(d__1)) >= (d__2 = dl[i__], f2c_dabs(d__2)))
            {
                /* No row interchange required */
                if(d__[i__] != 0.)
                {
                    fact = dl[i__] / d__[i__];
                    d__[i__ + 1] -= fact * du[i__];
                    b[i__ + 1 + b_dim1] -= fact * b[i__ + b_dim1];
                }
                else
                {
                    *info = i__;
                    AOCL_DTL_TRACE_LOG_EXIT
                    return;
                }
                dl[i__] = 0.;
            }
            else
            {
                /* Interchange rows I and I+1 */
                fact = d__[i__] / dl[i__];
                d__[i__] = dl[i__];
                temp = d__[i__ + 1];
                d__[i__ + 1] = du[i__] - fact * temp;
                dl[i__] = du[i__ + 1];
                du[i__ + 1] = -fact * dl[i__];
                du[i__] = temp;
                temp = b[i__ + b_dim1];
                b[i__ + b_dim1] = b[i__ + 1 + b_dim1];
                b[i__ + 1 + b_dim1] = temp - fact * b[i__ + 1 + b_dim1];
            }
            /* L10: */
        }
        if(*n > 1)
        {
            i__ = *n - 1;
            if((d__1 = d__[i__], f2c_dabs(d__1)) >= (d__2 = dl[i__], f2c_dabs(d__2)))
            {
                if(d__[i__] != 0.)
                {
                    fact = dl[i__] / d__[i__];
                    d__[i__ + 1] -= fact * du[i__];
                    b[i__ + 1 + b_dim1] -= fact * b[i__ + b_dim1];
                }
                else
                {
                    *info = i__;
                    AOCL_DTL_TRACE_LOG_EXIT
                    return;
                }
            }
            else
            {
                fact = d__[i__] / dl[i__];
                d__[i__] = dl[i__];
                temp = d__[i__ + 1];
                d__[i__ + 1] = du[i__] - fact * temp;
                du[i__] = temp;
                temp = b[i__ + b_dim1];
                b[i__ + b_dim1] = b[i__ + 1 + b_dim1];
                b[i__ + 1 + b_dim1] = temp - fact * b[i__ + 1 + b_dim1];
            }
        }
        if(d__[*n] == 0.)
        {
            *info = *n;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    else
    {
        i__1 = *n - 2;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            if((d__1 = d__[i__], f2c_dabs(d__1)) >= (d__2 = dl[i__], f2c_dabs(d__2)))
            {
                /* No row interchange required */
                if(d__[i__] != 0.)
                {
                    fact = dl[i__] / d__[i__];
                    d__[i__ + 1] -= fact * du[i__];
                    i__2 = *nrhs;
                    for(j = 1; j <= i__2; ++j)
                    {
                        b[i__ + 1 + j * b_dim1] -= fact * b[i__ + j * b_dim1];
                        /* L20: */
                    }
                }
                else
                {
                    *info = i__;
                    AOCL_DTL_TRACE_LOG_EXIT
                    return;
                }
                dl[i__] = 0.;
            }
            else
            {
                /* Interchange rows I and I+1 */
                fact = d__[i__] / dl[i__];
                d__[i__] = dl[i__];
                temp = d__[i__ + 1];
                d__[i__ + 1] = du[i__] - fact * temp;
                dl[i__] = du[i__ + 1];
                du[i__ + 1] = -fact * dl[i__];
                du[i__] = temp;
                i__2 = *nrhs;
                for(j = 1; j <= i__2; ++j)
                {
                    temp = b[i__ + j * b_dim1];
                    b[i__ + j * b_dim1] = b[i__ + 1 + j * b_dim1];
                    b[i__ + 1 + j * b_dim1] = temp - fact * b[i__ + 1 + j * b_dim1];
                    /* L30: */
                }
            }
            /* L40: */
        }
        if(*n > 1)
        {
            i__ = *n - 1;
            if((d__1 = d__[i__], f2c_dabs(d__1)) >= (d__2 = dl[i__], f2c_dabs(d__2)))
            {
                if(d__[i__] != 0.)
                {
                    fact = dl[i__] / d__[i__];
                    d__[i__ + 1] -= fact * du[i__];
                    i__1 = *nrhs;
                    for(j = 1; j <= i__1; ++j)
                    {
                        b[i__ + 1 + j * b_dim1] -= fact * b[i__ + j * b_dim1];
                        /* L50: */
                    }
                }
                else
                {
                    *info = i__;
                    AOCL_DTL_TRACE_LOG_EXIT
                    return;
                }
            }
            else
            {
                fact = d__[i__] / dl[i__];
                d__[i__] = dl[i__];
                temp = d__[i__ + 1];
                d__[i__ + 1] = du[i__] - fact * temp;
                du[i__] = temp;
                i__1 = *nrhs;
                for(j = 1; j <= i__1; ++j)
                {
                    temp = b[i__ + j * b_dim1];
                    b[i__ + j * b_dim1] = b[i__ + 1 + j * b_dim1];
                    b[i__ + 1 + j * b_dim1] = temp - fact * b[i__ + 1 + j * b_dim1];
                    /* L60: */
                }
            }
        }
        if(d__[*n] == 0.)
        {
            *info = *n;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    /* Back solve with the matrix U from the factorization. */
    if(*nrhs <= 2)
    {
        j = 1;
    L70:
        b[*n + j * b_dim1] /= d__[*n];
        if(*n > 1)
        {
            b[*n - 1 + j * b_dim1]
                = (b[*n - 1 + j * b_dim1] - du[*n - 1] * b[*n + j * b_dim1]) / d__[*n - 1];
        }
        for(i__ = *n - 2; i__ >= 1; --i__)
        {
            b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__ + 1 + j * b_dim1]
                                   - dl[i__] * b[i__ + 2 + j * b_dim1])
                                  / d__[i__];
            /* L80: */
        }
        if(j < *nrhs)
        {
            ++j;
            goto L70;
        }
    }
    else
    {
        i__1 = *nrhs;
        for(j = 1; j <= i__1; ++j)
        {
            b[*n + j * b_dim1] /= d__[*n];
            if(*n > 1)
            {
                b[*n - 1 + j * b_dim1]
                    = (b[*n - 1 + j * b_dim1] - du[*n - 1] * b[*n + j * b_dim1]) / d__[*n - 1];
            }
            for(i__ = *n - 2; i__ >= 1; --i__)
            {
                b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__ + 1 + j * b_dim1]
                                       - dl[i__] * b[i__ + 2 + j * b_dim1])
                                      / d__[i__];
                /* L90: */
            }
            /* L100: */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DGTSV */
}
/* dgtsv_ */
