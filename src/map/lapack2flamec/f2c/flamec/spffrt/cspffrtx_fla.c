/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    Oct 09, 2020
*/

#include "FLA_f2c.h"

extern int cspr_(char *, integer *, complex *, complex *, integer *, complex *);
extern void cscal_(integer *, complex *, complex *, integer *);
extern void c_div(complex *, complex *, complex *);

/*! @brief Partial LDL' factorization without pivoting
    *
    * @details
    * \b Purpose:
    * \verbatim
        CSPFFRTX computes the partial factorization of a complex symmetric matrix A
        stored in packed format.
        The factorization has the form
            A = L*D*L**T
        where L is a lower triangular matrix, and D is a diagonal matrix.
        This is an unblocked algorithm.
        The algorthm does not do pivoting and does not handle zero diagonal elements.
        Hence, it may give unexpected outputs for certain inputs.
    \endverbatim

    * @param[in,out] ap
    ap is COMPLEX array, dimension (N*(N+1)/2)
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
    work is COMPLEX array. \n
    Currently an unused buffer

    * @param[in] work2
    work2 is COMPLEX array. \n
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
    then CSPFFRTX performs ncolm successive factorizations of the form:

        ( a  b**T )   (  1   0 )   (  a        0       )   ( 1  (b/a)**T )
    A = (         ) = (        ) * (                   ) * (             )
        ( b    C  )   ( b/a  I )   (  0   C-b*1/a*b**T )   ( 0      I    )

    \endverbatim
    *  */

void  cspffrtx_fla(complex *ap, integer *n, integer *ncolm, complex *work, complex *work2 )
{
    complex z__1;
    integer i__1, k, kc;
    complex r1;
    complex c_b1 =
    {
        1.,0.
    };
    integer c__1 = 1;
    
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256, "cspffrtx inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
	  AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    --ap;
    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* K is the main loop index, increasing from 1 to ncolm in steps of */
    kc = 1;
    for( k = 1; k <= *ncolm; k++ )
    {
        /* Update the trailing submatrix */
        /* W(k) = L(k)*D(k) */
        /* where L(k) is the k-th column of L */

        c_div(&z__1, &c_b1, &ap[kc]);
        r1.r = -z__1.r;
        r1.i = -z__1.i;

        /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
        /* A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */
        i__1 = *n - k;
        cspr_("Lower", &i__1, &r1, &ap[kc + 1], &c__1, &ap[kc + *n - k + 1]);
        cscal_(&i__1, &z__1, &ap[kc + 1], &c__1);

        kc = kc + *n - k + 1;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}
/* cspffrtx_fla */
