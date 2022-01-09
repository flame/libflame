/*
    Copyright (c) 2021 Advanced Micro Devices, Inc.  All rights reserved.
    May 13, 2021
*/

#include "FLA_f2c.h"

extern void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
extern int zspr_(char *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *);
#ifdef FLA_ENABLE_BLAS_EXT_GEMMT
extern integer zgemmt_( char *, char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer * );
#endif

void zspffrt2_fla_def(doublecomplex *ap, integer *n, integer *ncolm, doublecomplex *work );
static void zspffrt2_fla_unp_var1( doublecomplex *ap, integer *n, integer *ncolm, doublecomplex *work );
static void zspffrt2_fla_unp_var2( doublecomplex *ap, integer *n, integer *ncolm, doublecomplex *work );

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
        ZSPFFRT2 computes the partial factorization of a complex symmetric matrix A
        stored in packed format.
        The factorization has the form
            A = L*D*L**T
        where L is a lower triangular matrix, and D is a diagonal matrix.
        This is an unblocked algorithm.
        The algorthm does not do pivoting and does not handle zero diagonal elements.
        Hence, it may give unexpected outputs for certain inputs.
    \endverbatim

    * @param[in,out] ap
    ap is COMPLEX*16 array, dimension (N*(N+1)/2)
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
    work is COMPLEX*16 array. \n
    Currently an unused buffer

    * @param[in] work2
    work2 is COMPLEX*16 array. \n
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
    then ZSPFFRT2 performs ncolm successive factorizations of the form:

        ( a  b**T )   ( a  0 )   ( 1/a       0       )   ( a  b**T )
    A = (         ) = (      ) * (                   ) * (         )
        ( b    C  )   ( b  I )   (  0   C-b*1/a*b**T )   ( 0    I  )

    \endverbatim
    *  */
void  zspffrt2_fla( doublecomplex *ap, integer *n, integer *ncolm, doublecomplex *work, doublecomplex *work2 )
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256, "zspffrt2 inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* ncolm as fraction of n */
    integer ncolm_pc = (integer) ( ( *ncolm * 100 ) / *n );
    if( ( *n > ( FLA_SPFFRT2__NTHRESH1 - 1 ) ) && ( ncolm_pc >= FLA_SPFFRT2__NCOLFRAC_THRESH3 ) )
    {
        /* Unpacking/packing based variant for small n &  ncolm values */
        zspffrt2_fla_unp_var2( ap, n, ncolm, work );
    }
    else if( *n > FLA_SPFFRT2__NTHRESH3 )
    {
        /* Unpacking/packing based variant for large n values */
        zspffrt2_fla_unp_var2( ap, n, ncolm, work );
    }
    else
    {
        /* zspr based implementation for small problem sizes */
        zspffrt2_fla_def( ap, n, ncolm, work );
    }

    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}

/*
 *  Unpacking the input packed symmetric matrix into lower
 *  triangular part of unpacked full matrix.
 *  The strictly upper triangular part is left untouched.
 */
void zunpack_fla( doublecomplex *ap, doublecomplex *a, integer m, integer n, integer lda )
{
   integer i, j;
   doublecomplex *aptr = ap;

   for( i = 0; i < n; i++ )
   {
      for( j = i; j < m; j++ )
      {
         a[ i * lda + j ] = *aptr++;
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
void zpack_fla( doublecomplex *a, doublecomplex *ap, integer m, integer n, integer lda )
{
   integer i, j;
   doublecomplex *aptr = ap;

   for( i = 0; i < n; i++ )
   {
      for( j = i; j < m; j++ )
      {
         *aptr++ = a[ i * lda + j ];
      }
   }

   return;
}

/*
 * LDLT factorization of skinny symmetric matrices (M > N)
 * in unpacked format.
 *
 * Only the lower trapezoidal part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 */
void zsffrk2_fla( doublecomplex *au, integer *m, integer *n, integer *lda, doublecomplex *bt, integer *ldbt )
{
    doublecomplex z__1;
    integer i__1, i__2, i__3;
    integer k, kc, kcn;
    integer c__1 = 1;
    doublecomplex r1;
    doublecomplex c_b1 = { 1., 0. };

    --au;
    --bt;
    kc = 1;
    i__3 = *m - *n;
    for( k = 1; k <= *n; k ++ )
    {
       /* D(k) = -1/A(k,k) */
       z_div(&z__1, &c_b1, &au[kc]);

       r1.r = -z__1.r;
       r1.i = -z__1.i;

       i__1 = *n - k;
       i__2 = *m - k;
       kcn = kc + *lda + 1;

       /* Update trailing matrix with rank-1 operation */
       zgeru_( &i__2, &i__1, &r1, &au[kc + 1], &c__1, &au[kc + 1], &c__1, &au[kcn], lda );

       /* Compute b**T/a */
       zcopy_( &i__3, &au[kc + *n - k + 1], &c__1, &bt[k], ldbt );
       zscal_( &i__3, &r1, &bt[k], ldbt );

       au[kc].r = z__1.r;
       au[kc].i = z__1.i;
       kc = kcn;
    }

    return;
}

/*
 * LDLT factorization of symmetric matrices (M > N)
 * in unpacked format.
 * Blocked algorithm employing GEMMT is used for better
 * performance.
 *
 * Only the lower triangular part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 *
 * Variant 1 does both factorization of (N x ncolm) and
 * trailing matrix update inside the main loop
 */
void zspffrt2_fla_unp_var1( doublecomplex *ap, integer *n, integer *ncolm, doublecomplex *work )
{
    doublecomplex c_b1 = { 1., 0. };
    integer k, kc;
    integer m, nb;

    doublecomplex *au, *bt;
    doublecomplex *mau, *mbt;

    /* Choose block size for the blocked variant */
    if( *n < FLA_SPFFRT2__BSIZE_NL1 )
        nb = FLA_SPFFRT2__BSIZE1;
    else if( *n < FLA_SPFFRT2__BSIZE_NL2 )
        nb = FLA_SPFFRT2__BSIZE2;
    else
        nb = FLA_SPFFRT2__BSIZE3;
    nb = (nb > *ncolm) ? *ncolm : nb;

    /* Allocate unpacked matrix and do the unpacking */
    mau = ( doublecomplex * ) malloc( *n * *n * sizeof( doublecomplex ) );
    mbt = ( doublecomplex * ) malloc( nb * (*n - nb) * sizeof( doublecomplex ) );

    zunpack_fla(ap, mau, *n, *n, *n );

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
        zsffrk2_fla( &au[kc], &m, &nb, n, &bt[1], &nb );
        m -= nb;

        /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        zgemm_( "N", "N", &m, &m, &nb, &c_b1, &au[kc + nb], n, &bt[1], &nb, &c_b1, &au[kc + nb * *n + nb], n );
#else
        zgemmt_( "L", "N", "N", &m, &nb, &c_b1, &au[kc + nb], n, &bt[1], &nb, &c_b1, &au[kc + nb * *n + nb], n );
#endif
    }

    /* Process the remaining columns */
    if( k <= *ncolm )
    {
        kc = k * *n - *n + k;
        nb = *ncolm - k + 1;

        /* Panel factorization using unblocked variant */
        zsffrk2_fla( &au[kc], &m, &nb, n, &bt[1], &nb );
        m -= nb;

        /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        zgemm_( "N", "N", &m, &m, &nb, &c_b1, &au[kc + nb], n, &bt[1], &nb, &c_b1, &au[kc + nb * *n + nb], n );
#else
        zgemmt_( "L", "N", "N", &m, &nb, &c_b1, &au[kc + nb], n, &bt[1], &nb, &c_b1, &au[kc + nb * *n + nb], n );
#endif
    }

    /* Pack the updated matrix */
    zpack_fla( mau, &ap[1], *n, *n, *n );

    free( mau );
    free( mbt );
    return;
}

/*
 * LDLT factorization of symmetric matrices (M > N)
 * in unpacked format.
 * Blocked algorithm employing GEMMT is used for better
 * performance.
 *
 * Only the lower triangular part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 *
 * Variant 2 does factorization of (N x ncolm) in the main loop
 * and the trailing matrix is updated outside the main loop
 */
void zspffrt2_fla_unp_var2( doublecomplex *ap, integer *n, integer *ncolm, doublecomplex *work )
{
    doublecomplex d__1 = { 1., 0. };
    integer kc, mg, nb;
    integer k, ni, mp;

    doublecomplex *au;
    doublecomplex *mau;

    /* Choose block size for the blocked variant */
    if( *n < FLA_SPFFRT2__BSIZE_NL1 )
        nb = FLA_SPFFRT2__BSIZE1;
    else if( *n < FLA_SPFFRT2__BSIZE_NL2 )
        nb = FLA_SPFFRT2__BSIZE2;
    else
        nb = FLA_SPFFRT2__BSIZE3;
    nb = (nb > *ncolm) ? *ncolm : nb;

    /* Allocate unpacked matrix and do the unpacking */
    mau = ( doublecomplex * ) malloc( *n * *n * sizeof( doublecomplex ) );

    zunpack_fla( ap, mau, *n, *n, *n );

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
        zsffrk2_fla( &au[kc], &mp, &nb, n, &au[kc + nb * *n], n );

        /* Update trailing matrix within the panel */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        mg = *n - *ncolm + ni;
        zgemm_( "N", "N", &mg, &ni, &nb, &d__1, &au[kc + nb], n, &au[kc + nb * *n], n, &d__1, &au[kc + nb * *n + nb], n );
#else
        zgemmt_( "L", "N", "N", &ni, &nb, &d__1, &au[kc + nb], n, &au[kc + nb * *n], n, &d__1, &au[kc + nb * *n + nb], n );
        zgemm_( "N", "N", &mg, &ni, &nb, &d__1, &au[kc + ni + nb], n, &au[kc + nb * *n], n, &d__1, &au[kc + nb * *n + nb + ni], n );
#endif
    }

    /* Process the remaining columns */
    if( k <= *ncolm )
    {
        mp -= nb;
        kc = k * *n - *n + k;
        nb = *ncolm - k + 1;

        /* Panel factorization using unblocked variant */
        zsffrk2_fla( &au[kc], &mp, &nb, n, &au[kc + nb * *n], n );
    }

    /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
    mg = *n - *ncolm;
    zgemm_( "N", "N", &mg, &mg, ncolm, &d__1, &au[*ncolm + 1], n, &au[*ncolm * *n + 1], n, &d__1, &au[*ncolm + *ncolm * *n + 1], n );
#else
    zgemmt_( "L", "N", "N", &mg, ncolm, &d__1, &au[*ncolm + 1], n, &au[*ncolm * *n + 1], n, &d__1, &au[*ncolm + *ncolm * *n + 1], n );
#endif

    /* Pack the updated matrix */
    zpack_fla( mau, &ap[1], *n, *n, *n );

    free( mau );
    return;
}

/*
 * AXPY based optimized version of rank-1 update for packedm
 * matrices.
 */

int lzspr_( char *uplo, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *ap, doublecomplex *work )
{
    integer incw = 1;
    integer k, kn, nz;

    ap--;
    work--;
    x--;

    for( k = 1; k <= *n; k++ )
    {
        work[k].r = alpha->r * x[k].r - alpha->i * x[k].i;
        work[k].i = alpha->r * x[k].i + alpha->i * x[k].r;
    }

    kn = 1;
    for( k = 1; k <= *n; k++ )
    {
        nz = *n - k + 1;
        zaxpy_( &nz, &work[k], &x[k], &incw, &ap[kn], &incw );
        kn = kn + *n - k + 1;
    }

    return 0;

}

void zspffrt2_fla_def(doublecomplex *ap, integer *n, integer *ncolm, doublecomplex *work )
{
    doublecomplex z__1;
    integer i__1, k, kc;
    integer c__1 = 1;
    doublecomplex r1;
    doublecomplex c_b1 = { 1., 0. };

    --ap;
    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* K is the main loop index, increasing from 1 to ncolm in steps of 1 */
    kc = 1;
    for( k = 1; k <= *ncolm; k++ )
    {
       /* Update the trailing submatrix */
       /* W(k) = L(k)*D(k) */
       /* where L(k) is the k-th column of L */

       z_div(&z__1, &c_b1, &ap[kc]);

       r1.r = -z__1.r;
       r1.i = -z__1.i;

       /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
       /* A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */
       i__1 = *n - k;
       lzspr_("Lower", &i__1, &r1, &ap[kc + 1], &c__1, &ap[kc + *n - k + 1], work);

       ap[kc].r = z__1.r;
       ap[kc].i = z__1.i;

       kc = kc + *n - k + 1;
    }
    return;
}
/* zspffrt2_fla */
