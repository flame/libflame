/*
    Copyright (c) 2021-2022 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLA_f2c.h"

extern integer dspr_( char *, integer *, doublereal *, doublereal *, integer *, doublereal * );
extern int dgemm_(char *transa, char *transb, integer *m, integer * n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, integer *ldc);

#ifdef FLA_ENABLE_BLAS_EXT_GEMMT
extern integer dgemmt_( char *, char *, char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer * );
#endif

static void dspffrt2_fla_def( doublereal *ap, integer *n, integer *ncolm, doublereal *work );
static void dspffrt2_fla_unp_var1( doublereal *ap, integer *n, integer *ncolm, doublereal *work );
static void dspffrt2_fla_unp_var2( doublereal *ap, integer *n, integer *ncolm, doublereal *work );

extern void DTL_Trace(
    uint8 ui8LogLevel,
    uint8 ui8LogType,
    const int8 *pi8FileName,
    const int8 *pi8FunctionName,
    uint32 ui32LineNumber,
    const int8 *pi8Message);

/*! @brief Partial LDL' factorization without pivoting
    *
    * @details
    * \b Purpose:
    * \verbatim
        DSPFFRT2 computes the partial factorization of a real matrix A
        stored in packed format.
        The factorization has the form
            A = L*D*L**T
        where L is a lower triangular matrix, and D is a diagonal matrix.
        This is an unblocked algorithm.
        The algorthm does not do pivoting and does not handle zero diagonal elements.
        Hence, it may give unexpected outputs for certain inputs.
    \endverbatim

    * @param[in,out] ap
    ap is DOUBLE PRECISION array, dimension (N*(N+1)/2)
    On entry, the lower triangle of the symmetric matrix A, packed columnwise in a 
    linear array. The j-th column of A is stored in the array AP as follows:
            AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
    On exit, the block diagonal matrix D and the multipliers used
    to obtain the factor L, stored as a packed triangular matrix overwriting A
    (see below for further details).

    * @param[in] n
    n is integer*. \n
    The order of the matrix A. *n >= 0

    * @param[in] ncolm
    ncolm is integer*. \n
    The number of columns / rows to be factorized. 0 <= *ncolm <= *n

    * @param[in] work
    work is DOUBLE PRECISION array. \n
    Currently an unused buffer

    * @param[in] work2
    work2 is DOUBLE PRECISION array. \n
    Currently an unused buffer

\par Further Details:
   ===================

    * \verbatim
    If input matrix A is of the form,

        ( a  b**T )
    A = (         )
        ( b    C  ), where

    a is the first diagonal element  A, b is a column vector of size n - 1 containing the
    elements from the first column of A excluding the diagonal element,
    C is the lower-right square submatrix of A, and I is the identity matrix,
    then DSPFFRT2 performs ncolm successive factorizations of the form:

        ( a  b**T )   ( a  0 )   ( 1/a       0       )   ( a  b**T )
    A = (         ) = (      ) * (                   ) * (         )
        ( b    C  )   ( b  I )   (  0   C-b*1/a*b**T )   ( 0    I  )

    \endverbatim
    *  */
void dspffrt2_fla( doublereal *ap, integer *n, integer * ncolm, doublereal *work, doublereal *work2 )
{
    /* ncolm as fraction of n */
    integer ncolm_pc = (integer) ( ( *ncolm * 100 ) / *n );

    if( ( *n < FLA_SPFFRT2__NTHRESH1 ) || ( *n < FLA_SPFFRT2__NTHRESH2 && ncolm_pc < FLA_SPFFRT2__NCOLFRAC_THRESH1 ) || ( *ncolm <= FLA_SPFFRT2__NCOLTHRESH ) )
    {
        /* dspr based implementation for small problem sizes */
        dspffrt2_fla_def( ap, n, ncolm, work );
    }
    else if( ncolm_pc < FLA_SPFFRT2__NCOLFRAC_THRESH2 )
    {
        /* Unpacking/packing based variant for large ncolm values */
        dspffrt2_fla_unp_var1( ap, n, ncolm, work );
    }
    else
    {
        /* Unpacking/packing based variant for smaller ncolm values */
        dspffrt2_fla_unp_var2( ap, n, ncolm, work );
    }
    return;
}

/*
 *  Unpacking the input packed symmetric matrix into lower
 *  triangular part of unpacked full matrix.
 *  The strictly upper triangular part is left untouched.
 */
void dunpack_fla( doublereal *ap, doublereal *a, integer m, integer n, integer lda )
{
   integer i, j;
   doublereal *aptr = ap;

   for( i = 0; i < n; i++ )
   {
      for( j = i; j < m; j++ )
      {
         a[i * lda + j] = *aptr++;
      }
   }

   return;
}

/*
 *  Packing the lower triangular part of input symmetric matrix into
 *  packed matrix.
 *  The strictly upper triangular parts of the input and output are
 *  left unused and untouched respectiely.
 */
void dpack_fla( doublereal *a, doublereal *ap, integer m, integer n, integer lda )
{
   integer i, j;
   doublereal *aptr = ap;

   for( i = 0; i < n; i++ )
   {
      for( j = i; j < m; j++ )
      {
         *aptr++ = a[i * lda + j];
      }
   }

   return;
}

/*
 * LDLT factorization of skinny symmetric matrices (m > n)
 * in unpacked format.
 *
 * Only the lower trapezoidal part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 */
void dsffrk2_fla( doublereal *au, integer *m, integer *n, integer *lda, doublereal *bt, integer *ldbt )
{
    doublereal d__1;
    integer i__1, i__2, i__3;
    integer k, kc, kcn;
    integer c__1 = 1;
    doublereal r1;
    extern int dger_(integer *m, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *a, integer *lda);
    extern int dcopy_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
    extern int dscal_(integer *n, doublereal *da, doublereal *dx, integer *incx);

    --au;
    --bt;
    kc = 1;
    i__3 = *m - *n;
    for( k = 1; k <= *n; k ++ )
    {
        /* D(k) = -1/A(k,k) */
        r1 = 1. / au[kc];
        d__1 = -r1;

        i__1 = *n - k;
        i__2 = *m - k;
        kcn = kc + *lda + 1;

        /* Update trailing matrix with rank-1 operation */
        dger_( &i__2, &i__1, &d__1, &au[kc + 1], &c__1, &au[kc + 1], &c__1, &au[kcn], lda );

        /* Compute b**T/a */
        dcopy_( &i__3, &au[kc + *n - k + 1], &c__1, &bt[k], ldbt );
        dscal_( &i__3, &d__1, &bt[k], ldbt );

        au[kc] = r1;
        kc = kcn;
    }

    return;
}

/*
 * LDLT factorization of symmetric matrices in unpacked format.
 * Blocked algorithm employing GEMMT is used for better
 * performance.
 *
 * Only the lower triangular part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 *
 * Variant 1 does both factorization of (N x ncolm) and
 * trailing matrix update inside the main loop
 */

void dspffrt2_fla_unp_var1( doublereal *ap, integer *n, integer *ncolm, doublereal *work )
{
    doublereal d__1 = 1.0;
    integer k, kc;
    integer m, nb;

    doublereal *au, *bt;
    doublereal *mau, *mbt;

    /* Choose block size for the blocked variant */
    if( *n < FLA_SPFFRT2__BSIZE_NL1 )
        nb = FLA_SPFFRT2__BSIZE1;
    else if( *n < FLA_SPFFRT2__BSIZE_NL2 )
        nb = FLA_SPFFRT2__BSIZE2;
    else
        nb = FLA_SPFFRT2__BSIZE3;
    nb = ( nb > *ncolm ) ? *ncolm : nb;

    /* Allocate unpacked matrix and do the unpacking */
    mau = ( doublereal * ) malloc( *n * *n * sizeof( doublereal ) );
    mbt = ( doublereal * ) malloc( nb * (*n - nb) * sizeof( doublereal ) );

    dunpack_fla( ap, mau, *n, *n, *n );

    --ap;
    au = mau - 1;
    bt = mbt - 1;

    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* k is the main loop index, increasing from 1 to ncolm in steps of nb */
    m  = *n;
    for( k = 1; k <= ( *ncolm - nb ); k += nb )
    {
        kc = k * *n - *n + k;

        /* Panel factorization using unblocked variant */
        dsffrk2_fla( &au[kc], &m, &nb, n, &bt[1], &nb );
        m -= nb;

        /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        dgemm_( "N", "N", &m, &m, &nb, &d__1, &au[kc + nb], n, &bt[1], &nb, &d__1, &au[kc + nb * *n + nb], n );
#else
        dgemmt_( "L", "N", "N", &m, &nb, &d__1, &au[kc + nb], n, &bt[1], &nb, &d__1, &au[kc + nb * *n + nb], n );
#endif
    }

    /* Process the remaining columns */
    if( k <= *ncolm )
    {
        kc = k * *n - *n + k;
        nb = *ncolm - k + 1;

        /* Panel factorization using unblocked variant */
        dsffrk2_fla( &au[kc], &m, &nb, n, &bt[1], &nb );
        m -= nb;

        /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        dgemm_( "N", "N", &m, &m, &nb, &d__1, &au[kc + nb], n, &bt[1], &nb, &d__1, &au[kc + nb * *n + nb], n );
#else
        dgemmt_( "L", "N", "N", &m, &nb, &d__1, &au[kc + nb], n, &bt[1], &nb, &d__1, &au[kc + nb * *n + nb], n );
#endif
    }

    /* Pack the updated matrix */
    dpack_fla( mau, &ap[1], *n, *n, *n );

    free( mau );
    free( mbt );
    return;
}

/*
 * LDLT factorization of symmetric matrices in unpacked format.
 * Blocked algorithm employing GEMMT is used for better
 * performance.
 *
 * Only the lower triangular part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 *
 * Variant 2 does factorization of (N x ncolm) in the main loop
 * and the trailing matrix is updated outside the main loop
 */

void dspffrt2_fla_unp_var2( doublereal *ap, integer *n, integer *ncolm, doublereal *work )
{
    doublereal d__1 = 1.0;
    integer kc, mg, nb;
    integer k, ni, mp;

    doublereal *au;
    doublereal *mau;

    /* Choose block size for the blocked variant */
    if( *n < FLA_SPFFRT2__BSIZE_NL1 )
        nb = FLA_SPFFRT2__BSIZE1;
    else if( *n < FLA_SPFFRT2__BSIZE_NL2 )
        nb = FLA_SPFFRT2__BSIZE2;
    else
        nb = FLA_SPFFRT2__BSIZE3;
    nb = ( nb > *ncolm ) ? *ncolm : nb;

    /* Allocate unpacked matrix and do the unpacking */
    mau = ( doublereal * ) malloc( *n * *n * sizeof( doublereal ) );

    dunpack_fla( ap, mau, *n, *n, *n );

    --ap;
    au = mau - 1;

    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* k is the main loop index, increasing from 1 to ncolm in steps of nb */
    mp = *n + nb;
    mg = *n - *ncolm;
    ni = *ncolm;
    for( k = 1; k <= ( *ncolm - nb ); k += nb )
    {
        mp -= nb;
        kc = k * *n - *n + k;
        ni = ni - nb;

        /* Panel factorization using unblocked variant */
        dsffrk2_fla( &au[kc], &mp, &nb, n, &au[kc + nb * *n], n );

        /* Update trailing matrix within the panel */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        mg = *n - *ncolm + ni;
        dgemm_( "N", "N", &mg, &ni, &nb, &d__1, &au[kc + nb], n, &au[kc + nb * *n], n, &d__1, &au[kc + nb * *n + nb], n );
#else
        dgemmt_( "L", "N", "N", &ni, &nb, &d__1, &au[kc + nb], n, &au[kc + nb * *n], n, &d__1, &au[kc + nb * *n + nb], n );
        dgemm_( "N", "N", &mg, &ni, &nb, &d__1, &au[kc + ni + nb], n, &au[kc + nb * *n], n, &d__1, &au[kc + nb * *n + nb + ni], n );
#endif
    }

    /* Process the remaining columns */
    if( k <= *ncolm )
    {
        mp -= nb;
        kc = k * *n - *n + k;
        nb = *ncolm - k + 1;

        /* Panel factorization using unblocked variant */
        dsffrk2_fla( &au[kc], &mp, &nb, n, &au[kc + nb * *n], n );
    }

    /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
    mg = *n - *ncolm;
    dgemm_( "N", "N", &mg, &mg, ncolm, &d__1, &au[*ncolm + 1], n, &au[*ncolm * *n + 1], n, &d__1, &au[*ncolm + *ncolm * *n + 1], n );
#else
    dgemmt_( "L", "N", "N", &mg, ncolm, &d__1, &au[*ncolm + 1], n, &au[*ncolm * *n + 1], n, &d__1, &au[*ncolm + *ncolm * *n + 1], n );
#endif

    /* Pack the updated matrix */
    dpack_fla( mau, &ap[1], *n, *n, *n );

    free( mau );
    return;
}

void dspffrt2_fla_def( doublereal *ap, integer *n, integer * ncolm, doublereal *work )
{
    doublereal d__1;
    integer i__1, k, kc;
    integer c__1 = 1;
    doublereal r1;

    --ap;

    /* k is the main loop index, increasing from 1 to ncolm in steps of 1 */
    kc = 1;
    for( k = 1; k <= *ncolm; k++ )
    {
       /* D(k) = -1/A(k,k) */
       r1 = 1. / ap[kc];
       d__1 = -r1;

       /* Update the trailing submatrix */
       /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
       /* A := A - L(k)*D(k)*L(k)**T */
       i__1 = *n - k;
       dspr_( "Lower", &i__1, &d__1, &ap[kc + 1], &c__1, &ap[kc + *n - k + 1] );

       ap[kc] = r1;
       kc = kc + *n - k + 1;
    }
    return;
}

/* dspffrt2_fla */
