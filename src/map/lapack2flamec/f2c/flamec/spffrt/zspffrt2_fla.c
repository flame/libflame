/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    Oct 09, 2020
*/

#include "FLA_f2c.h"
 
extern void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
extern integer zspr_(char *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *);

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
extern void DTL_Trace(
    uint8 ui8LogLevel,
    uint8 ui8LogType,
    const int8 *pi8FileName,
    const int8 *pi8FunctionName,
    uint32 ui32LineNumber,
    const int8 *pi8Message);

void  zspffrt2_fla(doublecomplex *ap, integer *n, integer *ncolm, doublecomplex *work, doublecomplex *work2 )
{
    doublecomplex z__1;
    integer i__1, k, kc;
    doublecomplex r1;
    doublecomplex c_b1 =
    {
        1.,0.
    };
    integer c__1 = 1;
    
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
    #if AOCL_DTL_LOG_ENABLE
    	char buffer[256];
	sprintf(buffer, "zspffrt2 inputs: n %d, ncolm %d\n", *n, *ncolm);
	AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
    #endif
    
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
       zspr_("Lower", &i__1, &r1, &ap[kc + 1], &c__1, &ap[kc + *n - k + 1]);

       ap[kc].r = z__1.r;
       ap[kc].i = z__1.i;

       kc = kc + *n - k + 1;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}
/* zspffrt2_fla */
