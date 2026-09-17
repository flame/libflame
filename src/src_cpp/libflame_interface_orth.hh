/******************************************************************************
 * Copyright (C) 2021-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file libflame_interface_orth.hh
 *  libflame_interface.hh defines all the orthogonal/unitary factors routines for
 *  Libflame CPP templated public interfaces.
 *  */
#ifndef LIBFLAME_INTERFACE_ORTH_HH
#define LIBFLAME_INTERFACE_ORTH_HH

namespace libflame
{
    /** @defgroup Orthogonal Orthogonal/Unitary factors(QR, CS, etc)
     * @ingroup LAPACK
     * @{
     */
    /** @defgroup QR QR
     * @ingroup Orthogonal
     * @{
     */

    /** @defgroup geqr geqr
     * @ingroup QR
     * @{
     */
    /*! @brief GEQR computes a QR factorization of M-by-N matrix
* @details
* \b Purpose:
   \verbatim
    GEQR computes a QR factorization of M-by-N matrix A:

    A = Q * ( R),
            ( 0)

    where:
    Q is a M-by-M orthogonal matrix;
    R is an upper-triangular N-by-N matrix;
    0 is a (M-N)-by-N zero matrix, if M > N.
   \endverbatim

* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
         N is INTEGER \n
         The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
         A is array, dimension (LDA,N) \n
         On entry, the M-by-N matrix A. \n
         On exit, the elements on and above the diagonal of the array
         contain the min(M,N)-by-N upper trapezoidal matrix R
         (R is upper triangular if M >= N);
         the elements below the diagonal are used to store part of the
         data structure to represent Q. \n
* @param[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A.  LDA >= fla_max(1,M). \n
* @param[out] T
         T is array, dimension (MAX(5,TSIZE)) \n
         On exit, if INFO = 0, T(1) returns optimal (or either minimal
         or optimal, if query is assumed) TSIZE. See TSIZE for details.
         Remaining T contains part of the data structure used to represent Q.
         If one wants to apply or construct Q, then one needs to keep T
         (in addition to A) and pass it to further subroutines. \n
* @param[in] TSIZE
         TSIZE is INTEGER \n
         If TSIZE >= 5, the dimension of the array T. \n
         If TSIZE = -1 or -2, then a workspace query is assumed. The routine
         only calculates the sizes of the T and WORK arrays, returns these
         values as the first entries of the T and WORK arrays, and no error
         message related to T or WORK is issued by XERBLA. \n
         If TSIZE = -1, the routine calculates optimal size of T for the
         optimum performance and returns this value in T(1). \n
         If TSIZE = -2, the routine calculates minimal size of T and
         returns this value in T(1). \n
* @param[out]	WORK
         (workspace) REAL array, dimension (MAX(1,LWORK)) \n
         On exit, if INFO = 0, WORK(1) contains optimal (or either minimal
         or optimal, if query was assumed) LWORK.
         See LWORK for details. \n
* @param[in]	LWORK
         LWORK is INTEGER \n
         The dimension of the array WORK. \n
         If LWORK = -1 or -2, then a workspace query is assumed. The routine
         only calculates the sizes of the T and WORK arrays, returns these
         values as the first entries of the T and WORK arrays, and no error
         message related to T or WORK is issued by XERBLA. \n
         If LWORK = -1, the routine calculates optimal size of WORK for the
         optimal performance and returns this value in WORK(1). \n
         If LWORK = -2, the routine calculates minimal size of WORK and
         returns this value in WORK(1). \n
* @param[out]	INFO
         INFO is INTEGER \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void geqr(integer *m, integer *n, T *a, integer *lda, T *t, integer *tsize, T *work,
              integer *lwork, integer *info)
    {
        geqr(m, n, a, lda, t, tsize, work, lwork, info);
    }
    /** @} */ // end of  geqr

    /** @defgroup gemqr gemqr
     * @ingroup QR
     * @{
     */
    /*! @brief Multiples a matrix C by a real orthogonal or complex unitary matrix Q, as computed by
 ?geqr,
  \n with best performance for tall and skinny matrices.
 * @details
 * \b Purpose:
    \verbatim
     GEMQR overwrites the general real M-by-N matrix C with

                    SIDE = 'L'     SIDE = 'R'
    TRANS = 'N':      Q * C          C * Q
    TRANS = 'T':      Q**T * C       C * Q**T
    where Q is a real orthogonal matrix defined as the product
    of blocked elementary reflectors computed by short wide LQ
    factorization (GEQR)
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q**T from the Left; \n
          = 'R': apply Q or Q**T from the Right. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N':  No transpose, apply Q; \n
          = 'T':  Transpose, apply Q**T. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >=0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines
          the matrix Q. \n
          If SIDE = 'L', M >= K >= 0; \n
          if SIDE = 'R', N >= K >= 0. \n
 * @param[in] A
          A is REAL array, dimension \n
                               (LDA,M) if SIDE = 'L', \n
                               (LDA,N) if SIDE = 'R' \n
          Part of the data structure to represent Q as returned by DGELQ. \n
 * @para@[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,K). \n
 * @param[in] T
          T is REAL array, dimension (MAX(5,TSIZE)). \n
          Part of the data structure to represent Q as returned by SGELQ. \n
 * @param[in] TSIZE
          TSIZE is INTEGER \n
          The dimension of the array T. TSIZE >= 5. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the M-by-N matrix C. \n
          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out]	WORK
         (workspace) COMPLEX array, dimension (MAX(1,LWORK)) \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If LWORK = -1, then a workspace query is assumed. The routine
          only calculates the size of the WORK array, returns this
          value as WORK(1), and no error message related to WORK
          is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void gemqr(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *t, integer *tsize, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        gemqr(side, trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info);
    }
    /** @} */ // end of gemqr

    /** @defgroup  geqrf geqrf
     * @ingroup QR
     * @{
     */
    /*! @brief QR factorization of a real m-by-n matrix a
    *
    * @details
    * \b Purpose:
    * \verbatim
        QR factorization of a real m-by-n matrix a
        The factorization has the form
            A = Q * R
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows of the matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the m-by-n matrix to be factored. \n
              On exit, the elements on and above the diagonal of the array
              contain the min(m,n)-by-n upper trapezoidal matrix R (R is
              upper triangular if m >= n); the elements below the diagonal,
              with the array tau, represent the orthogonal matrix Q as a
              product of min(m,n) elementary reflectors (see Further
              Details). \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,m) \n
    * @param[out] tau
              tau is float/double/COMPLEX/COMPLEX*16 array, dimension (min(m,n)) \n
              The scalar factors of the elementary reflectors (see Further
              Details). \n
              *
              * \n
              * **Further Details**
              * \verbatim
                      The matrix Q is represented as a product of elementary reflectors
                      Q = H(1) H(2) . . . H(k), where k = min(m,n).

                      Each H(i) has the form
                      H(i) = I - tau * V * V**T

                      where, tau is a real scalar, and V is a real vector with V(1:i-1) = 0 and V(i)
    = 1; V(i+1:M) is stored on exit in A(i+1:M,i), and tau in tau(i). \endverbatim
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK.  LWORK >= fla_max(1,N).
              For optimum performance LWORK >= N*NB, where NB is
              the optimal blocksize. \n

              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void geqrf(integer *m, integer *n, T *a, integer *lda, T *tau, T *work, integer *lwork,
               integer *info)
    {
        geqrf(m, n, a, lda, tau, work, lwork, info);
    }
    /** @} */ // end of geqrf

    /** @defgroup geqr2 geqr2
     * @ingroup QR
     * @{
     */
    /*! @brief QR factorization of a real m-by-n matrix a
    *
    * @details
    * \b Purpose:
    * \verbatim
        QR factorization of a real m-by-n matrix a
        The factorization has the form
            A = Q * R
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows of the matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the m-by-n matrix to be factored. \n
              On exit, the elements on and above the diagonal of the array
              contain the min(m,n)-by-n upper trapezoidal matrix R (R is
              upper triangular if m >= n); the elements below the diagonal,
              with the array tau, represent the orthogonal matrix Q as a
              product of elementary reflectors (see Further Details). \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,m) \n
    * @param[out] tau
              tau is float/double/COMPLEX/COMPLEX*16 array, dimension (min(m,n)) \n
              The scalar factors of the elementary reflectors (see Further
              Details). \n
              *
              * \n
              * **Further Details**
              * \verbatim
                      The matrix Q is represented as a product of elementary reflectors
                      Q = H(1) H(2) . . . H(k), where k = min(m,n).

                      Each H(i) has the form
                      H(i) = I - tau * V * V**T

                      where, tau is a real scalar, and V is a real vector with V(1:i-1) = 0 and V(i)
    = 1; V(i+1:M) is stored on exit in A(i+1:M,i), and tau in tau(i). \endverbatim
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (N) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void geqr2(integer *m, integer *n, T *a, integer *lda, T *tau, T *work, integer *info)
    {
        geqr2(m, n, a, lda, tau, work, info);
    }
    /** @} */ // end of geqr2

    /** @defgroup orgqr orgqr
     * @ingroup QR
     * @{
     */
    /*! @brief Form Q from QR factorization
    *
    * @details
    * \b Purpose:
    * \verbatim
        Generation of an m-by-n real matrix Q with orthonormal columns, which is defined as the
        first N columns of a product of K elementary reflectors of order M

        Q  =  H(1) H(2) . . . H(k)

        as returned by SGEQRF.
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows of the matrix Q. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix Q. m >= n >= 0. \n
    * @param[in] k
              k is integer* \n
              The number of elementary reflectors whose product defines the
              matrix Q. n >= k >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the i-th column must contain the vector which
              defines the elementary reflector H(i), for i = 1,2,...,k, as
              returned by SGEQRF in the first k columns of its array
              argument a. \n
              On exit, the m-by-n matrix Q. \n
    * @param[in] lda
              lda is integer* \n
              The first dimension of the array a. lda >= fla_max(1,m). \n
    * @param[in] tau
              tau is float/double array, dimension (k) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i), as returned by SGEQRF. \n
    * @param[out]	WORK
              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. LWORK >= fla_max(1,N).
              For optimum performance LWORK >= N*NB, where NB is the
              optimal blocksize. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument has an illegal value \n

    *     *  */
    template <typename T>
    void orgqr(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        orgqr(m, n, k, a, lda, tau, work, lwork, info);
    }
    /** @} */ // end of orgqr

    /** @defgroup  ungqr ungqr
     * @ingroup QR
     * @{
     */
    /*! @brief Form Q from QR factorization
    *
    * @details
    * \b Purpose:
    * \verbatim
        Generation of an m-by-n complex matrix Q with orthonormal columns, which is defined as the
        first N columns of a product of K elementary reflectors of order M

        Q  =  H(1) H(2) . . . H(k)

        as returned by SGEQRF.
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows of the matrix Q. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix Q. m >= n >= 0. \n
    * @param[in] k
              k is integer* \n
              The number of elementary reflectors whose product defines the
              matrix Q. n >= k >= 0. \n
    * @param[in,out] a
              a is COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the i-th column must contain the vector which
              defines the elementary reflector H(i), for i = 1,2,...,k, as
              returned by CGEQRF in the first k columns of its array
              argument a. \n
              On exit, the m-by-n matrix Q. \n
    * @param[in] lda
              lda is integer* \n
              The first dimension of the array a. lda >= fla_max(1,m). \n
    * @param[in] tau
              tau is COMPLEX/COMPLEX*16 array, dimension (k) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i), as returned by SGEQRF. \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. LWORK >= fla_max(1,N).
              For optimum performance LWORK >= N*NB, where NB is the
              optimal blocksize. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument has an illegal value \n

    *     *  */
    template <typename T>
    void ungqr(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        ungqr(m, n, k, a, lda, tau, work, lwork, info);
    }
    /** @} */ // end of ungqr

    /** @defgroup ung2r {un,or}g2r
     * @ingroup QR
     * @{
     */
    /*! @brief ORG2R generates all or part of the orthogonal matrix Q from a QR factorization
 determined by sgeqrf.

 * @details
 * \b Purpose:
    \verbatim
    ORG2R generates an m by n real matrix Q with orthonormal columns,
    which is defined as the first n columns of a product of k elementary
    reflectors of order m

          Q  =  H(1) H(2) . . . H(k)

    as   returned by GEQRF.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix Q. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix Q. M >= N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines the
          matrix Q. N >= K >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the i-th column must contain the vector which
          defines the elementary reflector H(i), for i = 1,2,...,k, as
          returned by SGEQRF in the first k columns of its array
          argument A. \n
          On exit, the m-by-n matrix Q. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The first dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[in] TAU
          TAU is REAL array, dimension (K) \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as   returned by SGEQRF. \n
 * @param[out] WORK
          WORK is REAL array, dimension (N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument has an illegal value \n

 *  * */
    template <typename T>
    void org2r(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *info)
    {
        org2r(m, n, k, a, lda, tau, work, info);
    }
    template <typename T>
    void ung2r(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *info)
    {
        ung2r(m, n, k, a, lda, tau, work, info);
    }
    /** @} */ // end of ung2r

    /** @defgroup unmqr unmqr
     * @ingroup QR
     * @{
     */
    /*! @brief Apply Q or Q' from QR factorization
    *
    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from QR factorization
        Overwrite the general complex m-by-n matrix c with

        side = 'L'  side = 'R'
        trans = 'N':   Q * C C * Q
        trans = 'H':   Q**H * C C * Q**H

        where Q is a complex unitary matrix defined as the product of k elementary reflectors

        Q = H(1) H(2) . . . H(k)

        as returned by CGEQRF. Q is of order M if side = 'L' and of order N if side = 'R'.
    \endverbatim

    * @param[in] side
              side is char* \n
              = 'L': apply Q or Q**H from the Left; \n
              = 'R': apply Q or Q**H from the Right. \n
    * @param[in] trans
              trans is char* \n
              = 'N':  No transpose, apply Q; \n
              = 'H':  Transpose, apply Q**H. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix c. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix c. n >= 0. \n
    * @param[in] k
              k is integer* \n
              The number of elementary reflectors whose product defines
              the matrix Q. \n
              If side = 'L', m >= k >= 0; \n
              if side = 'R', n >= k >= 0. \n
    * @param[in] a
              a is COMPLEX/COMPLEX*16 array, dimension (lda,k) \n
              The i-th column must contain the vector which defines the
              elementary reflector H(i), for i = 1,2,...,k, as returned by
              CGEQRF in the first k columns of its array argument a. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. \n
              If side = 'L', lda >= fla_max(1,m); \n
              if side = 'R', lda >= fla_max(1,n). \n
    * @param[in] tau
              tau is COMPLEX/COMPLEX*16 array, dimension (k) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i), as returned by SGEQRF. \n
    * @param[in,out] c
              c is COMPLEX/COMPLEX*16 array, dimension (ldc,n) \n
              On entry, the m-by-n matrix c. \n
              On exit, c is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. \n
    * @param[in] ldc
              ldc is integer* \n
              The leading dimension of the array c. ldc >= fla_max(1,m). \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. \n
              If SIDE = 'L', LWORK >= fla_max(1,N); \n
              if SIDE = 'R', LWORK >= fla_max(1,M). \n
              For good performance, LWORK should generally be larger. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void unmqr(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        unmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
    /** @} */ // end of unmqr

    /** @defgroup ormqr ormqr
     * @ingroup QR
     * @{
     */
    /*! @brief Apply Q or Q' from QR factorization
    *
    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from QR factorization
        Overwrite the general real m-by-n matrix c with

        side = 'L'  side = 'R'
        trans = 'N':   Q * C C * Q
        trans = 'T':   Q**T * C C * Q**T

        where Q is a real orthogonal matrix defined as the product of k elementary reflectors

        Q = H(1) H(2) . . . H(k)

        as returned by SGEQRF. Q is of order M if side = 'L' and of order N if side = 'R'.
    \endverbatim

    * @param[in] side
              side is char* \n
              = 'L': apply Q or Q**T from the Left; \n
              = 'R': apply Q or Q**T from the Right. \n
    * @param[in] trans
              trans is char* \n
              = 'N':  No transpose, apply Q; \n
              = 'T':  Transpose, apply Q**T. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix c. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix c. n >= 0. \n
    * @param[in] k
              k is integer* \n
              The number of elementary reflectors whose product defines
              the matrix Q. \n
              If side = 'L', m >= k >= 0; \n
              if side = 'R', n >= k >= 0. \n
    * @param[in] a
              a is float/double array, dimension (lda,k) \n
              The i-th column must contain the vector which defines the
              elementary reflector H(i), for i = 1,2,...,k, as returned by
              SGEQRF in the first k columns of its array argument a. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. \n
              If side = 'L', lda >= fla_max(1,m); \n
              if side = 'R', lda >= fla_max(1,n). \n
    * @param[in] tau
              tau is float/double array, dimension (k) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i), as returned by SGEQRF. \n
    * @param[in,out] c
              c is float/double array, dimension (ldc,n) \n
              On entry, the m-by-n matrix c. \n
              On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
    * @param[in] ldc
              ldc is integer* \n
              The leading dimension of the array c. ldc >= fla_max(1,m). \n
    * @param[out]	WORK
              WORK is REAL array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. \n
              If SIDE = 'L', LWORK >= fla_max(1,N); \n
              if SIDE = 'R', LWORK >= fla_max(1,M). \n
              For good performance, LWORK should generally be larger. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void ormqr(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        ormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
    /** @} */ // end of ormqr

    /** @defgroup orm2r orm2r
     * @ingroup QR
     * @{
     */
    /*! @brief Multiply a general matrix by the orthogonal matrix from a QR factorization
    *   determined by sgeqrf (unblocked algorithm).
    *
    * @details
    * \b Purpose:
    * \verbatim
        Multiply a general matrix by the orthogonal matrix from a QR factorization determined by
        sgeqrf (unblocked algorithm).
        Overwrite the general real m by n matrix c with

        Q * C  if side = 'L' and trans = 'N', or

        Q**T* C  if side = 'L' and trans = 'T', or

        C * Q  if side = 'R' and trans = 'N', or

        C * Q**T if side = 'R' and trans = 'T',

        where Q is a real orthogonal matrix defined as the product of k elementary reflectors

        Q = H(1) H(2) . . . H(k)

        as returned by SGEQRF. Q is of order m if side = 'L' and of order n if side = 'R'.
    \endverbatim

    * @param[in] side
              side is char* \n
              = 'L': apply Q or Q**T from the Left; \n
              = 'R': apply Q or Q**T from the Right. \n
    * @param[in] trans
              trans is char* \n
              = 'N':  No transpose, apply Q; \n
              = 'T':  Transpose, apply Q**T. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix c. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix c. n >= 0. \n
    * @param[in] k
              k is integer* \n
              The number of elementary reflectors whose product defines
              the matrix Q. \n
              If side = 'L', m >= k >= 0; \n
              if side = 'R', n >= k >= 0. \n
    * @param[in] a
              a is float/double array, dimension (lda,k) \n
              The i-th column must contain the vector which defines the
              elementary reflector H(i), for i = 1,2,...,k, as returned by
              SGEQRF in the first k columns of its array argument a. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. \n
              If side = 'L', lda >= fla_max(1,m); \n
              if side = 'R', lda >= fla_max(1,n). \n
    * @param[in] tau
              tau is float/double array, dimension (k) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i), as returned by SGEQRF. \n
    * @param[in,out] c
              c is float/double array, dimension (ldc,n) \n
              On entry, the m-by-n matrix c. \n
              On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
    * @param[in] ldc
              ldc is integer* \n
              The leading dimension of the array c. ldc >= fla_max(1,m). \n
    * @param[out]	WORK
              WORK is REAL array, dimension \n
                                       (N) if SIDE = 'L', \n
                                       (M) if SIDE = 'R' \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void orm2r(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *info)
    {
        orm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    }
    /** @} */ // end of  orm2r

    /** @defgroup unm2r unm2r
     * @ingroup QR
     * @{
     */
    /*! @brief Multiply a general matrix by the orthogonal matrix from a QR factorization
    *   determined by sgeqrf (unblocked algorithm).
    *
    * @details
    * \b Purpose:
    * \verbatim
        Multiply a general matrix by the orthogonal matrix from a QR factorization determined by
        sgeqrf (unblocked algorithm).
        Overwrite the general complex m by n matrix c with

        Q * C  if side = 'L' and trans = 'N', or

        Q**H* C  if side = 'L' and trans = 'C', or

        C * Q  if side = 'R' and trans = 'N', or

        C * Q**H if side = 'R' and trans = 'C',

        where Q is a complex unitary matrix defined as the product of k elementary reflectors

        Q = H(1) H(2) . . . H(k)

        as returned by SGEQRF. Q is of order m if side = 'L' and of order n if side = 'R'.
    \endverbatim

    * @param[in] side
              side is char* \n
              = 'L': apply Q or Q**H from the Left; \n
              = 'R': apply Q or Q**H from the Right. \n
    * @param[in] trans
              trans is char* \n
              = 'N':  No transpose, apply Q; \n
              = 'C':  Transpose, apply Q**H. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix c. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix c. n >= 0. \n
    * @param[in] k
              k is integer* \n
              The number of elementary reflectors whose product defines
              the matrix Q. \n
              If side = 'L', m >= k >= 0; \n
              if side = 'R', n >= k >= 0. \n
    * @param[in] a
              a is COMPLEX/COMPLEX*16 array, dimension (lda,k) \n
              The i-th column must contain the vector which defines the
              elementary reflector H(i), for i = 1,2,...,k, as returned by
              CGEQRF in the first k columns of its array argument a. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. \n
              If side = 'L', lda >= fla_max(1,m); \n
              if side = 'R', lda >= fla_max(1,n). \n
    * @param[in] tau
              tau is COMPLEX/COMPLEX*16 array, dimension (k) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i), as returned by CGEQRF. \n
    * @param[in,out] c
              c is COMPLEX/COMPLEX*16 array, dimension (ldc,n) \n
              On entry, the m-by-n matrix c. \n
              On exit, c is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. \n
    * @param[in] ldc
              ldc is integer* \n
              The leading dimension of the array c. ldc >= fla_max(1,m). \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension \n
                                       (N) if SIDE = 'L', \n
                                       (M) if SIDE = 'R' \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void unm2r(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *info)
    {
        unm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    }
    /** @} */ // end of unm2r

    /** @defgroup geqrt geqrt
     * @ingroup QR
     * @{
     */
    /*! @brief GEQRT computes a blocked QR factorization of a M-by-N matrix A \n
     using the compact WY representation
 *  @details
 *  \b Purpose:
    \verbatim
     GEQRT computes a blocked QR factorization of a M-by-N matrix A
     using the compact WY representation of Q.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] NB
          NB is INTEGER \n
          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit, the elements on and above the diagonal of the array
          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
          upper triangular if M >= N); the elements below the diagonal
          are the columns of V. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] T
          T is REAL array, dimension (LDT,MIN(M,N)) \n
          The upper triangular block reflectors stored in compact form
          as a sequence of upper triangular blocks.  See below
          for further details. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= NB. \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (NB*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void geqrt(integer *m, integer *n, integer *nb, T *a, integer *lda, T *t, integer *ldt, T *work,
               integer *info)
    {
        geqrt(m, n, nb, a, lda, t, ldt, work, info);
    }
    /** @} */ // end of geqrt

    /** @defgroup geqrt2 geqrt2
     * @ingroup QR
     * @{
     */
    /*! @brief GEQRT2 computes a blocked QR factorization of a M-by-N matrix A \n
     using the compact WY representation
 * @details
 * \b Purpose:
    \verbatim
     GEQRT2 computes a QR factorization of a real M-by-N matrix A,
     using the compact WY representation of Q.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= N. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the real M-by-N matrix A.  On exit, the elements on and
          above the diagonal contain the N-by-N upper triangular matrix R; the
          elements below the diagonal are the columns of V.  See below for
          further details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] T
          T is REAL array, dimension (LDT,N) \n
          The N-by-N upper triangular factor of the block reflector.
          The elements on and above the diagonal contain the block
          reflector T; the elements below the diagonal are not used.
          See below for further details. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void geqrt2(integer *m, integer *n, T *a, integer *lda, T *t, integer *ldt, integer *info)
    {
        geqrt2(m, n, a, lda, t, ldt, info);
    }
    /** @} */ // end of geqrt2

    /** @defgroup geqrt3 geqrt3
     * @ingroup QR
     * @{
     */
    /*! @brief GEQRT3 computes a blocked QR factorization of a M-by-N matrix A \n
     using the compact WY representation
 * @details
 * \b Purpose:
    \verbatim
     GEQRT3 computes a QR factorization of a real M-by-N matrix A,
     using the compact WY representation of Q.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= N. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the real M-by-N matrix A.  On exit, the elements on and
          above the diagonal contain the N-by-N upper triangular matrix R; the
          elements below the diagonal are the columns of V.  See below for
          further details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] T
          T is REAL array, dimension (LDT,N) \n
          The N-by-N upper triangular factor of the block reflector.
          The elements on and above the diagonal contain the block
          reflector T; the elements below the diagonal are not used.
          See below for further details. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void geqrt3(integer *m, integer *n, T *a, integer *lda, T *t, integer *ldt, integer *info)
    {
        geqrt3(m, n, a, lda, t, ldt, info);
    }
    /** @} */ // end of geqrt3

    /** @defgroup gemqrt gemqrt
     * @ingroup QR
     * @{
     */
    /*! @brief Multiplies a general matrix by the orthogonal/unitary matrix Q of the QR
factorization formed by ?geqrt.
 * @details
 * \b Purpose:
    \verbatim
     GEMQRT overwrites the general real M-by-N matrix C with

                    SIDE = 'L'     SIDE = 'R'
     TRANS = 'N':      Q C            C Q
     TRANS = 'T':   Q**T C            C Q**T

     where Q is a real orthogonal matrix defined as the product of K
     elementary reflectors:

          Q = H(1) H(2) . . . H(K) = I - V T V**T

     generated using the compact WY representation as returned by SGEQRT.

     Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q**T from the Left; \n
          = 'R': apply Q or Q**T from the Right. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N':  No transpose, apply Q; \n
          = 'T':  Transpose, apply Q**T. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines
          the matrix Q. \n
          If SIDE = 'L', M >= K >= 0; \n
          if SIDE = 'R', N >= K >= 0. \n
 * @param[in] NB
          NB is INTEGER \n
          The block size used for the storage of T.  K >= NB >= 1. \n
          This must be the same value of NB used to generate T
          in CGEQRT. \n
 * @param[in] V
          V is REAL array, dimension (LDV,K) \n
          The i-th column must contain the vector which defines the
          elementary reflector H(i), for i = 1,2,...,k, as returned by
          CGEQRT in the first K columns of its array argument A. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V. \n
          If SIDE = 'L', LDA >= fla_max(1,M); \n
          if SIDE = 'R', LDA >= fla_max(1,N). \n
 * @param[in] T
          T is REAL array, dimension (LDT,K) \n
          The upper triangular factors of the block reflectors
          as returned by CGEQRT, stored as a NB-by-N matrix. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= NB. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the M-by-N matrix C. \n
          On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out]	WORK
          WORK is COMPLEX array. The dimension of WORK is
           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

* * */
    template <typename T>
    void gemqrt(char *side, char *trans, integer *m, integer *n, integer *k, integer *nb, T *v,
                integer *ldv, T *t, integer *ldt, T *c, integer *ldc, T *work, integer *info)
    {
        gemqrt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info);
    }
    /** @} */ // end of gemqrt

    /** @defgroup geqrfp geqrfp
     * @ingroup QR
     * @{
     */
    /*! @brief GEQRP computes a QR factorization of a M-by-N matrix
 * @details
 * \b Purpose:
    \verbatim
     GEQRP computes a QR factorization of a M-by-N matrix A:
        A = Q * ( R),
                ( 0)

     where:
        Q is a M-by-M orthogonal matrix;
        R is an upper-triangular N-by-N matrix with nonnegative diagonal
        entries;
        0 is a (M-N)-by-N zero matrix, if M > N.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit, the elements on and above the diagonal of the array
          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
          upper triangular if m >= n). The diagonal entries of R
          are nonnegative; the elements below the diagonal,
          with the array TAU, represent the orthogonal matrix Q as a
          product of min(m,n) elementary reflectors (see Further
          Details). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] TAU
          TAU is array, dimension (min(M,N)) \n
          The scalar factors of the elementary reflectors (see Further
          Details). \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,N).
          For optimum performance LWORK >= N*NB, where NB is
          the optimal blocksize. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void geqrfp(integer *m, integer *n, T *a, integer *lda, T *tau, T *work, integer *lwork,
                integer *info)
    {
        geqrfp(m, n, a, lda, tau, work, lwork, info);
    }
    /** @} */ // end of geqrfp

    /** @defgroup geqr2p geqr2p
     * @ingroup QR
     * @{
     */
    /*! @brief GEQR2P computes the QR factorization of a general rectangular  \n
     matrix with non-negative diagonal elements using an unblocked algorithm

 * @details
 * \b Purpose:
    \verbatim
    GEQR2P computes a QR factorization of a real m-by-n matrix A:

       A = Q * (R),
               (0)

    where:

       Q is a m-by-m orthogonal matrix;
       R is an upper-triangular n-by-n matrix with nonnegative diagonal
       entries;
       0 is a (m-n)-by-n zero matrix, if m > n.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the m by n matrix A. \n
          On exit, the elements on and above the diagonal of the array
          contain the min(m,n) by n upper trapezoidal matrix R (R is
          upper triangular if m >= n). The diagonal entries of R
          are nonnegative; the elements below the diagonal,
          with the array TAU, represent the orthogonal matrix Q as a
          product of elementary reflectors (see Further Details). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] TAU
          TAU is REAL array, dimension (min(M,N)) \n
          The scalar factors of the elementary reflectors (see Further
          Details). \n
 * @param[out] WORK
          WORK is REAL array, dimension (N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void geqr2p(integer *m, integer *n, T *a, integer *lda, T *tau, T *work, integer *info)
    {
        geqr2p(m, n, a, lda, tau, work, info);
    }
    /** @} */ // end of geqr2p

    /** @defgroup geqp3 geqp3
     * @ingroup QR
     * @{
     */
    /*! @brief QR factorization of a real m-by-n matrix a
    *
    * @details
    * \b Purpose:
    * \verbatim
        QR factorization with column pivoting of a real m-by-n matrix a:
            A*P = Q*R
        using Level 3 BLAS.
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows of the matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the m-by-n matrix a. \n
              On exit, the upper triangle of the array contains the
              min(m,n)-by-n upper triangular matrix R; the elements
              below the diagonal, together with the array tau,
              represent the orthogonal matrix Q as a product of
              min(m,n) elementary reflectors. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,m) \n
    * @param[in,out] jpvt
              jpvt is integer array, dimension (n) \n
              On entry, if jpvt(J).ne.0, the J-th column of A is permuted
              to the front of A*P (a leading column); if jpvt(J) = 0,
              the J-th column of A is a free column. \n
              On exit, if jpvt(J) = K, then the J-th column of A*P
              was the K-th column of A. \n
    * @param[out] tau
              tau is float/double array, dimension (min(m,n)) \n
              The scalar factors of the elementary reflectors.
              *
              * \n
              * **Further Details**
              * \verbatim
                      The matrix Q is represented as a product of elementary reflectors

                      Q = H(1) H(2) . . . H(k), where k = min(m,n).

                      Each H(i) has the form

                      H = I - tau * V * V**T

                      where tau is a real scalar, and V is a real/complex vector with V(1:i-1) = 0
 and V(i) = 1; V(i+1:M) is stored on exit in A(i+1:M,i), and tau in tau(i) \endverbatim
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO=0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. LWORK >= N+1.
              For optimal performance LWORK >= ( N+1 )*NB, where NB
              is the optimal blocksize. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	RWORK
              RWORK is REAL array, dimension (2*N) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit. \n
              < 0: if INFO = -i, the i-th argument had an illegal value. \n

    *     *  */
    template <typename T>
    void geqp3(integer *m, integer *n, T *a, integer *lda, integer *jpvt, T *tau, T *work,
               integer *lwork, integer *info)
    {
        geqp3(m, n, a, lda, jpvt, tau, work, lwork, info);
    }
    template <typename T, typename Ta>
    void geqp3(integer *m, integer *n, T *a, integer *lda, integer *jpvt, T *tau, T *work,
               integer *lwork, Ta *rwork, integer *info)
    {
        geqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
    }
    /** @} */ // end of geqp3

    /** @defgroup laqp2 laqp2
     * @ingroup QR
     * @{
     */
    /*! @brief LAQP2 computes a QR factorization with column pivoting of the matrix block

 * @details
 * \b Purpose:
    \verbatim
    LAQP2 computes a QR factorization with column pivoting of
    the block A(OFFSET+1:M,1:N).
    The block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A. N >= 0. \n
 * @param[in] OFFSET
          OFFSET is INTEGER \n
          The number of rows of the matrix A that must be pivoted
          but no factorized. OFFSET >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit, the upper triangle of block A(OFFSET+1:M,1:N) is
          the triangular factor obtained; the elements in block
          A(OFFSET+1:M,1:N) below the diagonal, together with the
          array TAU, represent the orthogonal matrix Q as a product of
          elementary reflectors. Block A(1:OFFSET,1:N) has been
          accordingly pivoted, but no factorized. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[in,out] JPVT
          JPVT is INTEGER array, dimension (N) \n
          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
          to the front of A*P (a leading column); if JPVT(i) = 0,
          the i-th column of A is a free column. \n
          On exit, if JPVT(i) = k, then the i-th column of A*P
          was the k-th column of A. \n
 * @param[out] TAU
          TAU is REAL array, dimension (min(M,N)) \n
          The scalar factors of the elementary reflectors. \n
 * @param[in,out] VN1
          VN1 is REAL array, dimension (N) \n
          The vector with the partial column norms. \n
 * @param[in,out] VN2
          VN2 is REAL array, dimension (N) \n
          The vector with the exact column norms. \n
 * @param[out] WORK
          WORK is REAL array, dimension (N)  \n

 * */
    template <typename T>
    void laqp2(integer *m, integer *n, integer *offset, T *a, integer *lda, integer *jpvt, T *tau,
               T *vn1, T *vn2, T *work)
    {
        laqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
    }
    template <typename T, typename Ta>
    void laqp2(integer *m, integer *n, integer *offset, T *a, integer *lda, integer *jpvt, T *tau,
               Ta *vn1, Ta *vn2, T *work)
    {
        laqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
    }
    /** @} */ // end of laqp3

    /** @defgroup laqps laqps
     * @ingroup QR
     * @{
     */
    /*! @brief LAQPS computes a step of QR factorization with column pivoting of a real m-by-n
 matrix A by using BLAS level 3

 * @details
 * \b Purpose:
    \verbatim
    LAQPS computes a step of QR factorization with column pivoting
    of a real M-by-N matrix A by using Blas-3.  It tries to factorize
    NB columns from A starting from the row OFFSET+1, and updates all
    of the matrix with Blas-3 xGEMM.

    In some cases, due to catastrophic cancellations, it cannot
    factorize NB columns.  Hence, the actual number of factorized
    columns is   returned in KB.

    Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A. N >= 0 \n
 * @param[in] OFFSET
          OFFSET is INTEGER \n
          The number of rows of A that have been factorized in
          previous steps. \n
 * @param[in] NB
          NB is INTEGER \n
          The number of columns to factorize. \n
 * @param[out] KB
          KB is INTEGER \n
          The number of columns actually factorized. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit, block A(OFFSET+1:M,1:KB) is the triangular
          factor obtained and block A(1:OFFSET,1:N) has been
          accordingly pivoted, but no factorized.
          The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has
          been updated. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[in,out] JPVT
          JPVT is INTEGER array, dimension (N) \n
          JPVT(I) = K <==> Column K of the full matrix A has been
          permuted into position I in AP. \n
 * @param[out] TAU
          TAU is REAL array, dimension (KB) \n
          The scalar factors of the elementary reflectors. \n
 * @param[in,out] VN1
          VN1 is REAL array, dimension (N) \n
          The vector with the partial column norms. \n
 * @param[in,out] VN2
          VN2 is REAL array, dimension (N) \n
          The vector with the exact column norms. \n
 * @param[in,out] AUXV
          AUXV is REAL array, dimension (NB) \n
          Auxiliary vector. \n
 * @param[in,out] F
          F is REAL array, dimension (LDF,NB) \n
          Matrix F**T = L*Y**T*A. \n
 * @param[in] LDF
          LDF is INTEGER \n
          The leading dimension of the array F. LDF >= fla_max(1,N).  \n

 * */
    template <typename T>
    void laqps(integer *m, integer *n, integer *offset, integer *nb, integer *kb, T *a,
               integer *lda, integer *jpvt, T *tau, T *vn1, T *vn2, T *auxv, T *f, integer *ldf)
    {
        laqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
    }
    template <typename T, typename Ta>
    void laqps(integer *m, integer *n, integer *offset, integer *nb, integer *kb, T *a,
               integer *lda, integer *jpvt, T *tau, Ta *vn1, Ta *vn2, T *auxv, T *f, integer *ldf)
    {
        laqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
    }
    /** @} */ // end of laqps

    /** @defgroup latsqr latsqr
     * @ingroup QR
     * @{
     */
    /*! @brief LATSQR computes a blocked Tall-Skinny QR factorization of a real M-by-N matrix A for
 M >= N

 * @details
 * \b Purpose:
    \verbatim
    LATSQR computes a blocked Tall-Skinny QR factorization of
    a real M-by-N matrix A for M >= N:

       A = Q * (R),
               (0)

    where:

       Q is a M-by-M orthogonal matrix, stored on exit in an implicit
       form in the elements below the digonal of the array A and in
       the elemenst of the array T;

       R is an upper-triangular N-by-N matrix, stored on exit in
       the elements on and above the diagonal of the array A.

       0 is a (M-N)-by-N zero matrix, and is not stored.
    \endverbatim

  * @param[in] M
           M is INTEGER \n
           The number of rows of the matrix A.  M >= 0. \n
  * @param[in] N
           N is INTEGER \n
           The number of columns of the matrix A. M >= N >= 0. \n
  * @param[in] MB
           MB is INTEGER \n
           The row block size to be used in the blocked QR.
           MB > N. \n
  * @param[in] NB
           NB is INTEGER \n
           The column block size to be used in the blocked QR.
           N >= NB >= 1. \n
  * @param[in,out] A
           A is REAL array, dimension (LDA,N) \n
           On entry, the M-by-N matrix A. \n
           On exit, the elements on and above the diagonal
           of the array contain the N-by-N upper triangular matrix R;
           the elements below the diagonal represent Q by the columns
           of blocked V (see Further Details). \n
  * @param[in] LDA
           LDA is INTEGER \n
           The leading dimension of the array A.  LDA >= fla_max(1,M). \n
  * @param[out] T
           T is REAL array,
           dimension (LDT, N * Number_of_row_blocks) \n
           where Number_of_row_blocks = CEIL((M-N)/(MB-N)) \n
           The blocked upper triangular block reflectors stored in compact form
           as a sequence of upper triangular blocks.
           See Further Details below. \n
  * @param[in] LDT
           LDT is INTEGER \n
           The leading dimension of the array T.  LDT >= NB. \n
  * @param[out] WORK
          (workspace) REAL array, dimension (MAX(1,LWORK)) \n
  * @param[in] LWORK
           The dimension of the array WORK.  LWORK >= NB*N. \n
           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array, returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA. \n
  * @param[out] INFO
           INFO is INTEGER \n
           = 0:  successful exit \n
           < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void latsqr(integer *m, integer *n, integer *mb, integer *nb, T *a, integer *lda, T *t,
                integer *ldt, T *work, integer *lwork, integer *info)
    {
        latsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    }
    /** @} */ // end of latsqr

    /** @defgroup orgtsqr {un,or}gtsqr
     * @ingroup QR
     * @{
     */
    /*! @brief ORGTSQR generates an M-by-N real matrix Q_out with orthonormal columns

 * @details
 * \b Purpose:
    \verbatim
    ORGTSQR generates an M-by-N real matrix Q_out with orthonormal columns,
    which are the first N columns of a product of real orthogonal
    matrices of order M which are   returned by SLATSQR

         Q_out = first_N_columns_of(Q(1)_in * Q(2)_in * ... * Q(k)_in).

    See the documentation for LATSQR.
    \endverbatim

  * @param[in] M
           M is INTEGER \n
           The number of rows of the matrix A.  M >= 0. \n
  * @param[in] N
           N is INTEGER \n
           The number of columns of the matrix A. M >= N >= 0. \n
  * @param[in] MB
           MB is INTEGER \n
           The row block size used by SLATSQR to   return
           arrays A and T. MB > N. \n
           (Note that if MB > M, then M is used instead of MB
           as the row block size). \n
  * @param[in] NB
           NB is INTEGER \n
           The column block size used by SLATSQR to   return
           arrays A and T. NB >= 1. \n
           (Note that if NB > N, then N is used instead of NB
           as the column block size). \n
  * @param[in,out] A
           A is REAL array, dimension (LDA,N) \n
  \n
           On entry:
  \n
              The elements on and above the diagonal are not accessed.
              The elements below the diagonal represent the unit
              lower-trapezoidal blocked matrix V computed by SLATSQR
              that defines the input matrices Q_in(k) (ones on the
              diagonal are not stored) (same format as the output A
              below the diagonal in SLATSQR). \n
  \n
           On exit:
  \n
              The array A contains an M-by-N orthonormal matrix Q_out,
              i.e the columns of A are orthogonal unit vectors. \n
  * @param[in] LDA
           LDA is INTEGER \n
           The leading dimension of the array A.  LDA >= fla_max(1,M). \n
  * @param[in] T
           T is REAL array,
           dimension (LDT, N * NIRB) \n
           where NIRB = Number_of_input_row_blocks \n
                      = MAX(1, CEIL((M-N)/(MB-N))) \n
           Let NICB = Number_of_input_col_blocks \n
                    = CEIL(N/NB) \n
  \n
           The upper-triangular block reflectors used to define the
           input matrices Q_in(k), k=(1:NIRB*NICB). The block
           reflectors are stored in compact form in NIRB block
           reflector sequences. Each of NIRB block reflector sequences
           is stored in a larger NB-by-N column block of T and consists
           of NICB smaller NB-by-NB upper-triangular column blocks.
           (same format as the output T in SLATSQR). \n
  * @param[in] LDT
           LDT is INTEGER \n
           The leading dimension of the array T.
           LDT >= fla_max(1,min(NB1,N)). \n
  * @param[out] WORK
           (workspace) REAL array, dimension (MAX(2,LWORK)) \n
           On exit, if INFO = 0, WORK(1)   returns the optimal LWORK. \n
  * @param[in] LWORK
           The dimension of the array WORK.  LWORK >= (M+NB)*N. \n
           If LWORK = -1, then a workspace query is assumed. \n
           The routine only calculates the optimal size of the WORK
           array,   returns this value as the first entry of the WORK
           array, and no error message related to LWORK is issued
           by XERBLA. \n
  * @param[out] INFO
           INFO is INTEGER \n
           = 0:  successful exit \n
           < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void orgtsqr(integer *m, integer *n, integer *mb, integer *nb, T *a, integer *lda, T *t,
                 integer *ldt, T *work, integer *lwork, integer *info)
    {
        orgtsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    }
    template <typename T>
    void ungtsqr(integer *m, integer *n, integer *mb, integer *nb, T *a, integer *lda, T *t,
                 integer *ldt, T *work, integer *lwork, integer *info)
    {
        ungtsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    }
    /** @} */ // end of orgtsqr

    /** @defgroup  orgtqr_ orstqr_row
     * @ingroup QR
     * @{
     */
    /*! @brief ORGTSQR_ROW generates an M-by-N real matrix Q_out with
           orthonormal columns from the output of LATSQR.

 * @details
 * \b Purpose:
    \verbatim
    ORGTSQR_ROW generates an M-by-N real matrix Q_out with
    orthonormal columns from the output of LATSQR. These N orthonormal
    columns are the first N columns of a product of complex unitary
    matrices Q(k)_in of order M, which are returned by SLATSQR in
    a special format.

        Q_out = first_N_columns_of( Q(1)_in * Q(2)_in * ... * Q(k)_in ).

    The input matrices Q(k)_in are stored in row and column blocks in A.
    See the documentation of SLATSQR for more details on the format of
    Q(k)_in, where each Q(k)_in is represented by block Householder
    transformations. This routine calls an auxiliary routine SLARFB_GETT,
    where the computation is performed on each individual block. The
    algorithm first sweeps NB-sized column blocks from the right to left
    starting in the bottom row block and continues to the top row block
    (hence _ROW in the routine name). This sweep is in reverse order of
    the order in which SLATSQR generates the output blocks.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A. M >= N >= 0. \n
 * @param[in] MB
          MB is INTEGER \n
          The row block size used by SLATSQR to return
          arrays A and T. MB > N.
          (Note that if MB > M, then M is used instead of MB
          as the row block size). \n
 * @param[in] NB
          NB is INTEGER \n
          The column block size used by SLATSQR to return
          arrays A and T. NB >= 1. \n
          (Note that if NB > N, then N is used instead of NB
          as the column block size). \n
 * @param[in, out] A
          A is REAL array, dimension (LDA,N) \n
 \n
          On entry: \n
 \n
             The elements on and above the diagonal are not used as
             input. The elements below the diagonal represent the unit
             lower-trapezoidal blocked matrix V computed by SLATSQR
             that defines the input matrices Q_in(k) (ones on the
             diagonal are not stored). See SLATSQR for more details.
 \n
          On exit: \n
 \n
             The array A contains an M-by-N orthonormal matrix Q_out,
             i.e the columns of A are orthogonal unit vectors. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[in] T
          T is REAL array,
          dimension (LDT, N * NIRB) \n
          where NIRB = Number_of_input_row_blocks
                     = MAX( 1, CEIL((M-N)/(MB-N)) ) \n
          Let NICB = Number_of_input_col_blocks \n
                   = CEIL(N/NB) \n
 \n
          The upper-triangular block reflectors used to define the
          input matrices Q_in(k), k=(1:NIRB*NICB). The block
          reflectors are stored in compact form in NIRB block
          reflector sequences. Each of the NIRB block reflector
          sequences is stored in a larger NB-by-N column block of T
          and consists of NICB smaller NB-by-NB upper-triangular
          column blocks. See SLATSQR for more details on the format
          of T. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= fla_max(1,min(NB,N)). \n
 * @param[out] WORK
          (workspace) REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in] LWORK
          The dimension of the array WORK. \n
          LWORK >= NBLOCAL * MAX(NBLOCAL,(N-NBLOCAL)),
          where NBLOCAL=MIN(NB,N). \n
          If LWORK = -1, then a workspace query is assumed. \n
          The routine only calculates the optimal size of the WORK
          array, returns this value as the first entry of the WORK
          array, and no error message related to LWORK is issued
          by XERBLA. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void gtsqr_row(integer *m, integer *n, integer *mb, integer *nb, T *a, integer *lda, T *t,
                   integer *ldt, T *work, integer *lwork, integer *info)
    {
        gtsqr_row(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    }
    /** @} */ // end of

    /** @defgroup larfb_ larfb_gett
     * @ingroup QR
     * @{
     */
    /*! @brief LARFB_GETT applies a real Householder block reflector H from the
    left to a real (K+M)-by-N  "triangular-pentagonal" matrix.

 * @details
 * \b Purpose:
    \verbatim
    LARFB_GETT applies a real Householder block reflector H from the
    left to a real (K+M)-by-N  "triangular-pentagonal" matrix
    composed of two block matrices: an upper trapezoidal K-by-N matrix A
    stored in the array A, and a rectangular M-by-(N-K) matrix B, stored
    in the array B. The block reflector H is stored in a compact
    WY-representation, where the elementary reflectors are in the
    arrays A, B and T. See Further Details section.
    \endverbatim

 * @param[in] IDENT
          IDENT is CHARACTER*1 \n
          If IDENT = not 'I', or not 'i', then V1 is unit
             lower-triangular and stored in the left K-by-K block of
             the input matrix A, \n
          If IDENT = 'I' or 'i', then  V1 is an identity matrix and
             not stored. \n
          See Further Details section. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix B.
          M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrices A and B.
          N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number or rows of the matrix A. \n
          K is also order of the matrix T, i.e. the number of
          elementary reflectors whose product defines the block
          reflector. 0 <= K <= N. \n
 * @param[in] T
          T is REAL array, dimension (LDT,K) \n
          The upper-triangular K-by-K matrix T in the representation
          of the block reflector. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= K. \n
 * @param[in, out] A
          A is REAL array, dimension (LDA,N) \n
 \n
          On entry: \n
           a) In the K-by-N upper-trapezoidal part A: input matrix A. \n
           b) In the columns below the diagonal: columns of V1
              (ones are not stored on the diagonal). \n
 \n
          On exit: \n
            A is overwritten by rectangular K-by-N product H*A. \n
 \n
          See Further Details section. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,K). \n
 * @param[in, out] B
          B is REAL array, dimension (LDB,N) \n
 \n
          On entry: \n
            a) In the M-by-(N-K) right block: input matrix B. \n
            b) In the M-by-N left block: columns of V2. \n
 \n
          On exit: \n
            B is overwritten by rectangular M-by-N product H*B. \n
 \n
          See Further Details section. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,M). \n
 * @param[out] WORK
          WORK is REAL array, \n
          dimension (LDWORK,max(K,N-K)) \n
 * @param[in] LDWORK
          LDWORK is INTEGER \n
          The leading dimension of the array WORK. LDWORK>=max(1,K). \n

 *  * */
    template <typename T>
    void larfb_gett(char *ident, integer *m, integer *n, integer *k, T *t, integer *ldt, T *a,
                    integer *lda, T *b, integer *ldb, T *work, integer *ldwork)
    {
        larfb_gett(ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork);
    }
    /** @} */ // end of larf_gett

    /** @defgroup lamtsqr lamtsqr
     * @ingroup QR
     * @{
     */
    /*! @brief LAMTSQR overwrites the general real M-by-N matrix C

 * @details
 * \b Purpose:
    \verbatim
     LAMTSQR overwrites the general real M-by-N matrix C with

                    SIDE = 'L'     SIDE = 'R'
    TRANS = 'N':      Q * C          C * Q
    TRANS = 'T':      Q**T * C       C * Q**T
    where Q is a real orthogonal matrix defined as the product
    of blocked elementary reflectors computed by tall skinny
    QR factorization (DLATSQR)
    \endverbatim

  * @param[in] SIDE
           SIDE is CHARACTER*1 \n
           = 'L': apply Q or Q**T from the Left; \n
           = 'R': apply Q or Q**T from the Right. \n
  * @param[in] TRANS
           TRANS is CHARACTER*1 \n
           = 'N':  No transpose, apply Q; \n
           = 'T':  Transpose, apply Q**T. \n
  * @param[in] M
           M is INTEGER \n
           The number of rows of the matrix A.  M >=0. \n
  * @param[in] N
           N is INTEGER \n
           The number of columns of the matrix C. M >= N >= 0. \n
  * @param[in] K
           K is INTEGER \n
           The number of elementary reflectors whose product defines
           the matrix Q.
           N >= K >= 0; \n
  * @param[in] MB
           MB is INTEGER \n
           The block size to be used in the blocked QR.
           MB > N. (must be the same as DLATSQR) \n
  * @param[in] NB
           NB is INTEGER \n
           The column block size to be used in the blocked QR.
           N >= NB >= 1. \n
  * @param[in] A
           A is REAL array, dimension (LDA,K) \n
           The i-th column must contain the vector which defines the
           blockedelementary reflector H(i), for i = 1,2,...,k, as
           returned by DLATSQR in the first k columns of
           its array argument A. \n
  * @param[in] LDA
           LDA is INTEGER \n
           The leading dimension of the array A. \n
           If SIDE = 'L', LDA >= fla_max(1,M); \n
           if SIDE = 'R', LDA >= fla_max(1,N). \n
  * @param[in] T
           T is REAL array, dimension \n
           (N * Number of blocks(CEIL(M-K/MB-K)), \n
           The blocked upper triangular block reflectors stored in compact form
           as a sequence of upper triangular blocks.  See below
           for further details. \n
  * @param[in] LDT
           LDT is INTEGER \n
           The leading dimension of the array T.  LDT >= NB. \n
  * @param[in,out] C
           C is REAL array, dimension (LDC,N) \n
           On entry, the M-by-N matrix C.
           On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
  * @param[in] LDC
           LDC is INTEGER \n
           The leading dimension of the array C. LDC >= fla_max(1,M). \n
  * @param[out] WORK
          (workspace) REAL array, dimension (MAX(1,LWORK)) \n
  * @param[in] LWORK
           LWORK is INTEGER \n
           The dimension of the array WORK. \n
  \n
           If SIDE = 'L', LWORK >= fla_max(1,N)*NB; \n
           if SIDE = 'R', LWORK >= fla_max(1,MB)*NB. \n
           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array,   returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA. \n
  * @param[out] INFO
           INFO is INTEGER \n
           = 0:  successful exit \n
           < 0:  if INFO = -i, the i-th argument had an illegal value \n

 * */
    template <typename T>
    void lamtsqr(char *side, char *trans, integer *m, integer *n, integer *k, integer *mb,
                 integer *nb, T *a, integer *lda, T *t, integer *ldt, T *c, integer *ldc, T *work,
                 integer *lwork, integer *info)
    {
        lamtsqr(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
    }
    /** @} */ // end of lamtsqr

    /** @defgroup getsqrhrt getsqrhrt
     * @ingroup QR
     * @{
     */
    /*! @brief Computes a NB2-sized column blocked QR-factorization.

 * @details
 * \b Purpose:
    \verbatim
    The routine computes a NB2-sized column blocked QR-factorization of a
    complex M-by-N matrix A with M >= N,

    A = Q * R.

    The routine uses internally a NB1-sized column blocked and MB1-sized
    row blocked TSQR-factorization and perfors the reconstruction
    of the Householder vectors from the TSQR output. The routine also
    converts the R_tsqr factor from the TSQR-factorization output into
    the R factor that corresponds to the Householder QR-factorization,

      A = Q_tsqr * R_tsqr = Q * R.

    The output Q and R factors are stored in the same format as in CGEQRT
    (Q is in blocked compact WY-representation). See the documentation
    of CGEQRT for more details on the format.
    \endverbatim

 * @param[in] m
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] n
          N is INTEGER \n
          The number of columns of the matrix A. M >= N >= 0. \n
 * @param[in] mb1
          MB1 is INTEGER \n
          The row block size to be used in the blocked TSQR. \n
          MB1 > N. \n
 * @param[in, out] nb1
          NB1 is INTEGER \n
          The column block size to be used in the blocked TSQR. \n
          N >= NB1 >= 1. \n
 * @param[in, out] nb2
          NB2 is INTEGER \n
          The block size to be used in the blocked QR that is \n
          output. NB2 >= 1. \n
 * @param[in, out] a
          A is REAL/DOUBLE/COMPLEX/COMPLEX*16 array, dimension (LDA,N) \n

          On entry: an M-by-N matrix A.\n

          On exit:\n
           a) the elements on and above the diagonal
              of the array contain the N-by-N upper-triangular
              matrix R corresponding to the Householder QR; \n
           b) the elements below the diagonal represent Q by
              the columns of blocked V (compact WY-representation). \n
 * @param[in] lda
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[in, out] t
          T is REAL/DOUBLE/COMPLEX/COMPLEX*16 array, dimension (LDT,N)) \n
          The upper triangular block reflectors stored in compact form
          as a sequence of upper triangular blocks. \n
 * @param[in] ldt
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= NB2. \n
 * @param[out]	work
          WORK is REAL/DOUBLE/COMPLEX/COMPLEX*16 array, dimension (N,NRHS) \n
          This array is used to hold the residual vectors. \n
 * @param[out]	lwork
          The dimension of the array WORK.\n
          LWORK >= MAX( LWT + LW1, MAX( LWT+N*N+LW2, LWT+N*N+N ) ), \n
          where \n
             NUM_ALL_ROW_BLOCKS = CEIL((M-N)/(MB1-N)), \n
             NB1LOCAL = MIN(NB1,N). \n
             LWT = NUM_ALL_ROW_BLOCKS * N * NB1LOCAL, \n
             LW1 = NB1LOCAL * N, \n
             LW2 = NB1LOCAL * MAX( NB1LOCAL, ( N - NB1LOCAL ) ), \n
          If LWORK = -1, then a workspace query is assumed. \n
          The routine only calculates the optimal size of the WORK
          array, returns this value as the first entry of the WORK
          array, and no error message related to LWORK is issued
          by XERBLA. \n
 * @param[in] info
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void getsqrhrt(integer *m, integer *n, integer *mb1, integer *nb1, integer *nb2, T *a,
                   integer *lda, T *t, integer *ldt, T *work, integer *lwork, integer *info)
    {
        getsqrhrt(m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info);
    }
    /** @} */ // end of getsqrht

    /** @defgroup  orhr_col {un,or}hr_col
     * @ingroup QR
     * @{
     */
    /*! @brief ORHR_COL takes an M-by-N real matrix Q_in with orthonormal columns \n
     as input and performs Householder Reconstruction

 * @details
 * \b Purpose:
    \verbatim
    ORHR_COL takes an M-by-N real matrix Q_in with orthonormal columns
    as input, stored in A, and performs Householder Reconstruction (HR),
    i.e. reconstructs Householder vectors V(i) implicitly representing
    another M-by-N matrix Q_out, with the property that Q_in = Q_out*S,
    where S is an N-by-N diagonal matrix with diagonal entries
    equal to +1 or -1. The Householder vectors (columns V(i) of V) are
    stored in A on output, and the diagonal entries of S are stored in D.
    Block reflectors are also   returned in T
    (same output format as SGEQRT).
    \endverbatim

  * @param[in] M
           M is INTEGER \n
           The number of rows of the matrix A. M >= 0. \n
  * @param[in] N
           N is INTEGER \n
           The number of columns of the matrix A. M >= N >= 0. \n
  * @param[in] NB
           NB is INTEGER \n
           The column block size to be used in the reconstruction
           of Householder column vector blocks in the array A and
           corresponding block reflectors in the array T. NB >= 1.
           (Note that if NB > N, then N is used instead of NB
           as the column block size.) \n
  * @param[in,out] A
           A is REAL array, dimension (LDA,N) \n
  \n
           On entry:
  \n
              The array A contains an M-by-N orthonormal matrix Q_in,
              i.e the columns of A are orthogonal unit vectors. \n
  \n
           On exit:
  \n
              The elements below the diagonal of A represent the unit
              lower-trapezoidal matrix V of Householder column vectors
              V(i). The unit diagonal entries of V are not stored
              (same format as the output below the diagonal in A from
              SGEQRT). The matrix T and the matrix V stored on output
              in A implicitly define Q_out.
  \n
              The elements above the diagonal contain the factor U
              of the "modified" LU-decomposition: \n
                 Q_in - (S) = V * U \n
                        (0) \n
              where 0 is a (M-N)-by-(M-N) zero matrix. \n
  * @param[in] LDA
           LDA is INTEGER \n
           The leading dimension of the array A.  LDA >= fla_max(1,M). \n
  * @param[out] T
           T is REAL array,
           dimension (LDT, N) \n
  \n
           Let NOCB = Number_of_output_col_blocks \n
                    = CEIL(N/NB) \n
  \n
           On exit, T(1:NB, 1:N) contains NOCB upper-triangular
           block reflectors used to define Q_out stored in compact
           form as a sequence of upper-triangular NB-by-NB column
           blocks (same format as the output T in SGEQRT).
           The matrix T and the matrix V stored on output in A
           implicitly define Q_out. NOTE: The lower triangles
           below the upper-triangular blcoks will be filled with
           zeros. See Further Details. \n
  * @param[in] LDT
           LDT is INTEGER \n
           The leading dimension of the array T.
           LDT >= fla_max(1,min(NB,N)). \n
  * @param[out] D
           D is REAL array, dimension min(M,N). \n
           The elements can be only plus or minus one.
  \n
           D(i) is constructed as D(i) = -SIGN(Q_in_i(i,i)), where
           1 <= i <= min(M,N), and Q_in_i is Q_in after performing
           i-1 steps of â€œmodifiedâ€ Gaussian elimination.
           See Further Details. \n
  * @param[out] INFO
           INFO is INTEGER \n
           = 0:  successful exit \n
           < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void orhr_col(integer *m, integer *n, integer *nb, T *a, integer *lda, T *t, integer *ldt, T *d,
                  integer *info)
    {
        orhr_col(m, n, nb, a, lda, t, ldt, d, info);
    }
    template <typename T>
    void unhr_col(integer *m, integer *n, integer *nb, T *a, integer *lda, T *t, integer *ldt, T *d,
                  integer *info)
    {
        unhr_col(m, n, nb, a, lda, t, ldt, d, info);
    }
    /** @} */ // end of orhr_col

    /** @defgroup launhr_col_getrfnp launhr_col_getrfnp
     * @ingroup QR
     * @{
     */
    /*! @brief LAUNHR_COL_GETRFNP computes the modified LU factorization without pivoting of a
 complex general M-by-N matrix A

 * @details
 * \b Purpose:
    \verbatim
    LAUNHR_COL_GETRFNP computes the modified LU factorization without
    pivoting of a complex general M-by-N matrix A. The factorization has
    the form:

        A - S = L * U,

    where:
       S is a m-by-n diagonal sign matrix with the diagonal D, so that
       D(i) = S(i,i), 1 <= i <= min(M,N). The diagonal D is constructed
       as D(i)=-SIGN(A(i,i)), where A(i,i) is the value after performing
       i-1 steps of Gaussian elimination. This means that the diagonal
       element at each step of "modified" Gaussian elimination is
       at least one in absolute value (so that division-by-zero not
       not possible during the division by the diagonal element);

       L is a M-by-N lower triangular matrix with unit diagonal elements
       (lower trapezoidal if M > N);

       and U is a M-by-N upper triangular matrix
       (upper trapezoidal if M < N).

    This routine is an auxiliary routine used in the Householder
    reconstruction routine CUNHR_COL. In CUNHR_COL, this routine is
    applied to an M-by-N matrix A with orthonormal columns, where each
    element is bounded by one in absolute value. With the choice of
    the matrix S above, one can show that the diagonal element at each
    step of Gaussian elimination is the largest (in absolute value) in
    the column on or below the diagonal, so that no pivoting is required
    for numerical stability [1].

    For more details on the Householder reconstruction algorithm,
    including the modified LU factorization, see [1].

    This is the blocked right-looking version of the algorithm,
    calling Level 3 BLAS to update the submatrix. To factorize a block,
    this routine calls the recursive routine CLAUNHR_COL_GETRFNP2.

    [1] "Reconstructing Householder vectors from tall-skinny QR",
        G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen,
        E. Solomonik, J. Parallel Distrib. Comput.,
        vol. 85, pp. 3-31, 2015.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is COMPLEX array, dimension (LDA,N) \n
          On entry, the M-by-N matrix to be factored.
          On exit, the factors L and U from the factorization
          A-S=L*U; the unit diagonal elements of L are not stored. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] D
          D is COMPLEX array, dimension min(M,N) \n
          The diagonal elements of the diagonal M-by-N sign matrix S,
          D(i) = S(i,i), where 1 <= i <= min(M,N). The elements can be
          only (+1.0, 0.0) or (-1.0, 0.0). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 * */
    template <typename T>
    void launhr_col_getrfnp(integer *m, integer *n, T *a, integer *lda, T *d, integer *info)
    {
        launhr_col_getrfnp(m, n, a, lda, d, info);
    }
    /** @} */ // end of launhr_col_getrfnp

    /** @defgroup  laorhr_col_getrfnp laorhr_col_getrfnp
     * @ingroup QR
     * @{
     */
    /*! @brief LAORHR_COL_GETRFNP computes the modified LU factorization without pivoting of a real
 general M-by-N matrix A

 * @details
 * \b Purpose:
    \verbatim
    LAORHR_COL_GETRFNP computes the modified LU factorization without
    pivoting of a real general M-by-N matrix A. The factorization has
    the form:

        A - S = L * U,

    where:
       S is a m-by-n diagonal sign matrix with the diagonal D, so that
       D(i) = S(i,i), 1 <= i <= min(M,N). The diagonal D is constructed
       as D(i)=-SIGN(A(i,i)), where A(i,i) is the value after performing
       i-1 steps of Gaussian elimination. This means that the diagonal
       element at each step of "modified" Gaussian elimination is
       at least one in absolute value (so that division-by-zero not
       not possible during the division by the diagonal element);

       L is a M-by-N lower triangular matrix with unit diagonal elements
       (lower trapezoidal if M > N);

       and U is a M-by-N upper triangular matrix
       (upper trapezoidal if M < N).

    This routine is an auxiliary routine used in the Householder
    reconstruction routine SORHR_COL. In SORHR_COL, this routine is
    applied to an M-by-N matrix A with orthonormal columns, where each
    element is bounded by one in absolute value. With the choice of
    the matrix S above, one can show that the diagonal element at each
    step of Gaussian elimination is the largest (in absolute value) in
    the column on or below the diagonal, so that no pivoting is required
    for numerical stability [1].

    For more details on the Householder reconstruction algorithm,
    including the modified LU factorization, see [1].

    This is the blocked right-looking version of the algorithm,
    calling Level 3 BLAS to update the submatrix. To factorize a block,
    this routine calls the recursive routine SLAORHR_COL_GETRFNP2.

    [1] "Reconstructing Householder vectors from tall-skinny QR",
        G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen,
        E. Solomonik, J. Parallel Distrib. Comput.,
        vol. 85, pp. 3-31, 2015.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix to be factored.
          On exit, the factors L and U from the factorization
          A-S=L*U; the unit diagonal elements of L are not stored. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] D
          D is REAL array, dimension min(M,N) \n
          The diagonal elements of the diagonal M-by-N sign matrix S,
          D(i) = S(i,i), where 1 <= i <= min(M,N). The elements can
          be only plus or minus one. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value  \n

 * */
    template <typename T>
    void laorhr_col_getrfnp(integer *m, integer *n, T *a, integer *lda, T *d, integer *info)
    {
        laorhr_col_getrfnp(m, n, a, lda, d, info);
    }
    /** @} */ // end of laorhr_col_getrfnp

    /** @defgroup launhr_col_getrfnp2 launhr_col_getrfnp2
     * @ingroup QR
     * @{
     */
    /*! @brief LAUNHR_COL_GETRFNP2 computes the modified LU factorization without pivoting of a real
 general M-by-N matrix A

 * @details
 * \b Purpose:
    \verbatim
    LAUNHR_COL_GETRFNP2 computes the modified LU factorization without
    pivoting of a complex general M-by-N matrix A. The factorization has
    the form:

        A - S = L * U,

    where:
       S is a m-by-n diagonal sign matrix with the diagonal D, so that
       D(i) = S(i,i), 1 <= i <= min(M,N). The diagonal D is constructed
       as D(i)=-SIGN(A(i,i)), where A(i,i) is the value after performing
       i-1 steps of Gaussian elimination. This means that the diagonal
       element at each step of "modified" Gaussian elimination is at
       least one in absolute value (so that division-by-zero not
       possible during the division by the diagonal element);

       L is a M-by-N lower triangular matrix with unit diagonal elements
       (lower trapezoidal if M > N);

       and U is a M-by-N upper triangular matrix
       (upper trapezoidal if M < N).

    This routine is an auxiliary routine used in the Householder
    reconstruction routine CUNHR_COL. In CUNHR_COL, this routine is
    applied to an M-by-N matrix A with orthonormal columns, where each
    element is bounded by one in absolute value. With the choice of
    the matrix S above, one can show that the diagonal element at each
    step of Gaussian elimination is the largest (in absolute value) in
    the column on or below the diagonal, so that no pivoting is required
    for numerical stability [1].

    For more details on the Householder reconstruction algorithm,
    including the modified LU factorization, see [1].

    This is the recursive version of the LU factorization algorithm.
    Denote A - S by B. The algorithm divides the matrix B into four
    submatrices:

           [  B11 | B12  ]  where B11 is n1 by n1,
       B = [ -----|----- ]        B21 is (m-n1) by n1,
           [  B21 | B22  ]        B12 is n1 by n2,
                                  B22 is (m-n1) by n2,
                                  with n1 = min(m,n)/2, n2 = n-n1.


    The subroutine calls itself to factor B11, solves for B21,
    solves for B12, updates B22, then calls itself to factor B22.

    For more details on the recursive LU algorithm, see [2].

    CLAUNHR_COL_GETRFNP2 is called to factorize a block by the blocked
    routine CLAUNHR_COL_GETRFNP, which uses blocked code calling
   *. Level 3 BLAS to update the submatrix. However, CLAUNHR_COL_GETRFNP2
    is self-sufficient and can be used without CLAUNHR_COL_GETRFNP.

    [1] "Reconstructing Householder vectors from tall-skinny QR",
        G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen,
        E. Solomonik, J. Parallel Distrib. Comput.,
        vol. 85, pp. 3-31, 2015.

    [2] "Recursion leads to automatic variable blocking for dense linear
        algebra algorithms", F. Gustavson, IBM J. of Res. and Dev.,
        vol. 41, no. 6, pp. 737-755, 1997.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is COMPLEX array, dimension (LDA,N) \n
          On entry, the M-by-N matrix to be factored.
          On exit, the factors L and U from the factorization
          A-S=L*U; the unit diagonal elements of L are not stored. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] D
          D is COMPLEX array, dimension min(M,N) \n
          The diagonal elements of the diagonal M-by-N sign matrix S,
          D(i) = S(i,i), where 1 <= i <= min(M,N). The elements can be
          only (+1.0, 0.0) or (-1.0, 0.0). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 * */
    template <typename T>
    void launhr_col_getrfnp2(integer *m, integer *n, T *a, integer *lda, T *d, integer *info)
    {
        launhr_col_getrfnp2(m, n, a, lda, d, info);
    }
    /** @} */ // end of launhr_col_getrfnp2

    /** @defgroup laorhr_col_getrfnp2 laorhr_col_getrfnp2
     * @ingroup QR
     * @{
     */
    /*! @brief LAORHR_COL_GETRFNP2 computes the modified LU factorization without pivoting of a real
 general M-by-N matrix A

 * @details
 * \b Purpose:
    \verbatim
    LAORHR_COL_GETRFNP2 computes the modified LU factorization without
    pivoting of a real general M-by-N matrix A. The factorization has
    the form:

        A - S = L * U,

    where:
       S is a m-by-n diagonal sign matrix with the diagonal D, so that
       D(i) = S(i,i), 1 <= i <= min(M,N). The diagonal D is constructed
       as D(i)=-SIGN(A(i,i)), where A(i,i) is the value after performing
       i-1 steps of Gaussian elimination. This means that the diagonal
       element at each step of "modified" Gaussian elimination is at
       least one in absolute value (so that division-by-zero not
       possible during the division by the diagonal element);

       L is a M-by-N lower triangular matrix with unit diagonal elements
       (lower trapezoidal if M > N);

       and U is a M-by-N upper triangular matrix
       (upper trapezoidal if M < N).

    This routine is an auxiliary routine used in the Householder
    reconstruction routine SORHR_COL. In SORHR_COL, this routine is
    applied to an M-by-N matrix A with orthonormal columns, where each
    element is bounded by one in absolute value. With the choice of
    the matrix S above, one can show that the diagonal element at each
    step of Gaussian elimination is the largest (in absolute value) in
    the column on or below the diagonal, so that no pivoting is required
    for numerical stability [1].

    For more details on the Householder reconstruction algorithm,
    including the modified LU factorization, see [1].

    This is the recursive version of the LU factorization algorithm.
    Denote A - S by B. The algorithm divides the matrix B into four
    submatrices:

           [  B11 | B12  ]  where B11 is n1 by n1,
       B = [ -----|----- ]        B21 is (m-n1) by n1,
           [  B21 | B22  ]        B12 is n1 by n2,
                                  B22 is (m-n1) by n2,
                                  with n1 = min(m,n)/2, n2 = n-n1.


    The subroutine calls itself to factor B11, solves for B21,
    solves for B12, updates B22, then calls itself to factor B22.

    For more details on the recursive LU algorithm, see [2].

    SLAORHR_COL_GETRFNP2 is called to factorize a block by the blocked
    routine SLAORHR_COL_GETRFNP, which uses blocked code calling
   *. Level 3 BLAS to update the submatrix. However, SLAORHR_COL_GETRFNP2
    is self-sufficient and can be used without SLAORHR_COL_GETRFNP.

    [1] "Reconstructing Householder vectors from tall-skinny QR",
        G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen,
        E. Solomonik, J. Parallel Distrib. Comput.,
        vol. 85, pp. 3-31, 2015.

    [2] "Recursion leads to automatic variable blocking for dense linear
        algebra algorithms", F. Gustavson, IBM J. of Res. and Dev.,
        vol. 41, no. 6, pp. 737-755, 1997.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix to be factored.
          On exit, the factors L and U from the factorization
          A-S=L*U; the unit diagonal elements of L are not stored. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] D
          D is REAL array, dimension min(M,N) \n
          The diagonal elements of the diagonal M-by-N sign matrix S,
          D(i) = S(i,i), where 1 <= i <= min(M,N). The elements can
          be only plus or minus one. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 * */
    template <typename T>
    void laorhr_col_getrfnp2(integer *m, integer *n, T *a, integer *lda, T *d, integer *info)
    {
        laorhr_col_getrfnp2(m, n, a, lda, d, info);
    }
    /** @} */ // end of laorhr_col_getrfnp2

    /** @defgroup tpqrt tpqrt
     * @ingroup QR
     * @{
     */
    /*! @brief TPQRT computes a blocked QR factorization of a real
     triangular-pentagonal matrix
 * @details
 * \b Purpose:
    \verbatim
     TPQRT computes a blocked QR factorization of a real
     "triangular-pentagonal" matrix C, which is composed of a
     triangular block A and pentagonal block B, using the compact
     WY representation for Q.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix B. \n
          M >= 0.
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix B, and the order of the
          triangular matrix A.
          N >= 0. \n
 * @param[in] L
          L is INTEGER \n
          The number of rows of the upper trapezoidal part of B.
          MIN(M,N) >= L >= 0.  See Further Details. \n
 * @param[in] NB
          NB is INTEGER \n
          The block size to be used in the blocked QR.  N >= NB >= 1. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the upper triangular N-by-N matrix A.
          On exit, the elements on and above the diagonal of the array
          contain the upper triangular matrix R. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
          are rectangular, and the last L rows are upper trapezoidal. \n
          On exit, B contains the pentagonal matrix V.  See Further Details. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,M). \n
 * @param[out] T
          T is REAL array, dimension (LDT,N) \n
          The upper triangular block reflectors stored in compact form
          as a sequence of upper triangular blocks.  See Further Details. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= NB. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (NB*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tpqrt(integer *m, integer *n, integer *l, integer *nb, T *a, integer *lda, T *b,
               integer *ldb, T *t, integer *ldt, T *work, integer *info)
    {
        tpqrt(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info);
    }
    /** @} */ // end of tpqrt

    /** @defgroup tpqrt2 tpqrt2
     * @ingroup QR
     * @{
     */
    /*! @brief TPQRT2 computes a QR factorization of a real or complex    \n
     "triangular-pentagonal" matrix, which is composed of a triangular block \n
     and a pentagonal block, using the compact WY representation for Q
 * @details
 * \b Purpose:
    \verbatim
     TPQRT2 computes a QR factorization of a real "triangular-pentagonal"
     matrix C, which is composed of a triangular block A and pentagonal block B,
     using the compact WY representation for Q.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The total number of rows of the matrix B.
          M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix B, and the order of
          the triangular matrix A.
          N >= 0. \n
 * @param[in] L
          L is INTEGER \n
          The number of rows of the upper trapezoidal part of B.
          MIN(M,N) >= L >= 0.  See Further Details. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the upper triangular N-by-N matrix A.
          On exit, the elements on and above the diagonal of the array
          contain the upper triangular matrix R. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
          are rectangular, and the last L rows are upper trapezoidal.
          On exit, B contains the pentagonal matrix V.  See Further Details. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,M). \n
 * @param[out] T
          T is REAL array, dimension (LDT,N) \n
          The N-by-N upper triangular factor T of the block reflector.
          See Further Details. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= fla_max(1,N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tpqrt2(integer *m, integer *n, integer *l, T *a, integer *lda, T *b, integer *ldb, T *t,
                integer *ldt, integer *info)
    {
        tpqrt2(m, n, l, a, lda, b, ldb, t, ldt, info);
    }
    /** @} */ // end of tpqrt2

    /** @defgroup tpmqrt tpmqrt
     * @ingroup QR
     * @{
     */
    /*! @brief TPMQRT applies a real orthogonal matrix Q obtained from a   \n
     triangular-pentagonal
 * @details
 * \b Purpose:
    \verbatim
    TPMQRT applies a real orthogonal matrix Q obtained from a
    "triangular-pentagonal" real block reflector H to a general
    real matrix C, which consists of two blocks A and B.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q^T from the Left; \n
          = 'R': apply Q or Q^T from the Right. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N':  No transpose, apply Q; \n
          = 'T':  Transpose, apply Q^T. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix B. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix B. N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines
          the matrix Q. \n
 * @param[in] L
          L is INTEGER \n
          The order of the trapezoidal part of V.
          K >= L >= 0.  See Further Details. \n
 * @param[in] NB
          NB is INTEGER \n
          The block size used for the storage of T.  K >= NB >= 1.
          This must be the same value of NB used to generate T
          in CTPQRT. \n
 * @param[in] V
          V is REAL array, dimension (LDV,K) \n
          The i-th column must contain the vector which defines the
          elementary reflector H(i), for i = 1,2,...,k, as returned by
          CTPQRT in B.  See Further Details. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V. \n
          If SIDE = 'L', LDV >= fla_max(1,M); \n
          if SIDE = 'R', LDV >= fla_max(1,N). \n
 * @param[in] T
          T is REAL array, dimension (LDT,K) \n
          The upper triangular factors of the block reflectors
          as returned by CTPQRT, stored as a NB-by-K matrix. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= NB. \n
 * @param[in,out] A
          A is REAL array, dimension \n
          (LDA,N) if SIDE = 'L' or \n
          (LDA,K) if SIDE = 'R' \n
          On entry, the K-by-N or M-by-K matrix A.
          On exit, A is overwritten by the corresponding block of
          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A \n
          If SIDE = 'L', LDC >= fla_max(1,K); \n
          If SIDE = 'R', LDC >= fla_max(1,M). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the M-by-N matrix B. \n
          On exit, B is overwritten by the corresponding block of
          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.
          LDB >= fla_max(1,M). \n
 * @param[out]	WORK
          WORK is REAL array. The dimension of WORK is
           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tpmqrt(char *side, char *trans, integer *m, integer *n, integer *k, integer *l,
                integer *nb, T *v, integer *ldv, T *t, integer *ldt, T *a, integer *lda, T *b,
                integer *ldb, T *work, integer *info)
    {
        tpmqrt(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
    }
    /** @} */ // end of tpmqrt

    /** @defgroup tprfb tprfb
     * @ingroup QR
     * @{
     */
    /*! @brief TPRFB applies a real or complex "triangular-pentagonal" blocked  \n
     reflector to a real or complex matrix, which is composed of two blocks.
 * @details
 * \b Purpose:
    \verbatim
     TPRFB applies a real "triangular-pentagonal" block reflector H or its
     conjugate transpose H^H to a real matrix C, which is composed of two
     blocks A and B, either from the left or right.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply H or H^H from the Left \n
          = 'R': apply H or H^H from the Right \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N': apply H (No transpose) \n
          = 'C': apply H^H (Conjugate transpose) \n
 * @param[in] DIRECT
          DIRECT is CHARACTER*1 \n
          Indicates how H is formed from a product of elementary
          reflectors \n
          = 'F': H = H(1) H(2) . . . H(k) (Forward) \n
          = 'B': H = H(k) . . . H(2) H(1) (Backward) \n
 * @param[in] STOREV
          STOREV is CHARACTER*1 \n
          Indicates how the vectors which define the elementary
          reflectors are stored: \n
          = 'C': Columns \n
          = 'R': Rows \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix B.
          M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix B.
          N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The order of the matrix T, i.e. the number of elementary
          reflectors whose product defines the block reflector.
          K >= 0. \n
 * @param[in] L
          L is INTEGER \n
          The order of the trapezoidal part of V.
          K >= L >= 0.  See Further Details. \n
 * @param[in] V
          V is REAL array, dimension \n
                                (LDV,K) if STOREV = 'C' \n
                                (LDV,M) if STOREV = 'R' and SIDE = 'L' \n
                                (LDV,N) if STOREV = 'R' and SIDE = 'R' \n
          The pentagonal matrix V, which contains the elementary reflectors
          H(1), H(2), ..., H(K).  See Further Details. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V. \n
          If STOREV = 'C' and SIDE = 'L', LDV >= fla_max(1,M); \n
          if STOREV = 'C' and SIDE = 'R', LDV >= fla_max(1,N); \n
          if STOREV = 'R', LDV >= K. \n
 * @param[in] T
          T is REAL array, dimension (LDT,K) \n
          The triangular K-by-K matrix T in the representation of the
          block reflector. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.
          LDT >= K. \n
 * @param[in,out] A
          A is REAL array, dimension
          (LDA,N) if SIDE = 'L' or (LDA,K) if SIDE = 'R' \n
          On entry, the K-by-N or M-by-K matrix A. \n
          On exit, A is overwritten by the corresponding block of
          H*C or H^H*C or C*H or C*H^H.  See Further Details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. \n
          If SIDE = 'L', LDA >= fla_max(1,K); \n
          If SIDE = 'R', LDA >= fla_max(1,M). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the M-by-N matrix B.
          On exit, B is overwritten by the corresponding block of
          H*C or H^H*C or C*H or C*H^H.  See Further Details. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.
          LDB >= fla_max(1,M). \n
 * @param[out]	WORK
          WORK is REAL array, dimension \n
          (LDWORK,N) if SIDE = 'L', \n
          (LDWORK,K) if SIDE = 'R'. \n
 * @param[in]	LDWORK
          LDWORK is INTEGER \n
          The leading dimension of the array WORK. \n
          If SIDE = 'L', LDWORK >= K; \n
          if SIDE = 'R', LDWORK >= M. \n

 *  * */
    template <typename T>
    void tprfb(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
               integer *k, integer *l, T *v, integer *ldv, T *t, integer *ldt, T *a, integer *lda,
               T *b, integer *ldb, T *work, integer *ldwork)
    {
        tprfb(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work,
              ldwork);
    }

    /** @} */ // end of tprfb

    /** @defgroup ggqrf ggqrf
     * @ingroup QR
     * @{
     */
    /*! @brief GGQRF computes a generalized QR factorization of an N-by-M matrix A and an N-by-P
 matrix B
 * @details
 * \b Purpose:
    \verbatim
     GGQRF computes a generalized QR factorization of an N-by-M matrix A
     and an N-by-P matrix B:
                 A = Q*R,        B = Q*T*Z,

     where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal
     matrix, and R and T assume one of the forms:
     if N >= M,  R = ( R11) M  ,   or if N < M,  R = ( R11  R12) N,
                     (  0 ) N-M                         N   M-N
                        M
     where R11 is upper triangular, and
     if N <= P,  T = ( 0  T12) N,   or if N > P,  T = ( T11) N-P,
                      P-N  N                           ( T21) P
                                                          P
     where T12 or T21 is upper triangular.
     In particular, if B is square and nonsingular, the GQR factorization
     of A and B implicitly gives the QR factorization of inv(B)*A:
                  inv(B)*A = Z**T*(inv(T)*R)
     where inv(B) denotes the inverse of the matrix B, and Z**T denotes the
     transpose of the matrix Z.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of rows of the matrices A and B. N >= 0. \n
 * @param[in] M
          M is INTEGER \n
          The number of columns of the matrix A.  M >= 0. \n
 * @param[in] P
          P is INTEGER \n
          The number of columns of the matrix B.  P >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,M) \n
          On entry, the N-by-M matrix A. \n
          On exit, the elements on and above the diagonal of the array
          contain the min(N,M)-by-M upper trapezoidal matrix R (R is
          upper triangular if N >= M); the elements below the diagonal,
          with the array TAUA, represent the orthogonal matrix Q as a
          product of min(N,M) elementary reflectors (see Further
          Details). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,N). \n
  * @param[out] TAUA
          TAUA is REAL array, dimension (min(N,M)) \n
          The scalar factors of the elementary reflectors which
          represent the orthogonal matrix Q (see Further Details). \n
  * @param[in,out] B
          B is REAL array, dimension (LDB,P) \n
          On entry, the N-by-P matrix B. \n
          On exit, if N <= P, the upper triangle of the subarray
          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T;
          if N > P, the elements on and above the (N-P)-th subdiagonal
          contain the N-by-P upper trapezoidal matrix T; the remaining
          elements, with the array TAUB, represent the orthogonal
          matrix Z as a product of elementary reflectors (see Further
          Details). \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,N). \n
 * @param[out] TAUB
          TAUB is REAL array, dimension (min(N,P)) \n
          The scalar factors of the elementary reflectors which
          represent the orthogonal matrix Z (see Further Details). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= fla_max(1,N,M,P).
          For optimum performance LWORK >= fla_max(N,M,P)*max(NB1,NB2,NB3),
          where NB1 is the optimal blocksize for the QR factorization
          of an N-by-M matrix, NB2 is the optimal blocksize for the
          RQ factorization of an N-by-P matrix, and NB3 is the optimal
          blocksize for a call of SORMQR. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void ggqrf(integer *n, integer *m, integer *p, T *a, integer *lda, T *taua, T *b, integer *ldb,
               T *taub, T *work, integer *lwork, integer *info)
    {
        ggqrf(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
    }
    /** @} */ // end of ggqrf
    /**@} */ // end of QR

    /** @defgroup LQ LQ
     * @ingroup Orthogonal
     * @{
     */
    /** @defgroup gelq gelq
     * @ingroup LQ
     * @{
     */
    /*! @brief Computes LQ factorization of a real matrix
* @details
* \b Purpose:
* \verbatim
   GELQ computes an LQ factorization of a real M-by-N matrix A:

     A = ( L 0) *  Q

   where:

   Q is a N-by-N orthogonal matrix;
   L is an lower-triangular M-by-M matrix;
   0 is a M-by-(N-M) zero matrix, if M < N.
  \endverbatim

* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix A.  M >= 0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix A.  N >= 0. \n
* @param[in,out] A
         A is REAL array, dimension (LDA,N) \n
         On entry, the M-by-N matrix A. \n
         On exit, the elements on and below the diagonal of the array
         contain the M-by-min(M,N) lower trapezoidal matrix L
         (L is lower triangular if M <= N);
         the elements above the diagonal are used to store part of the
         data structure to represent -
         computes row and column scaling to reduce condition number of matrix \n
* @param[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A.  LDA >= fla_max(1,M). \n
* @param[out] T
         T is REAL array, dimension (MAX(5,TSIZE)) \n
         On exit, if INFO = 0, T(1) returns optimal (or either minimal
         or optimal, if query is assumed) TSIZE. See TSIZE for details.
         Remaining T contains part of the data structure used to represent Q.
         If one wants to apply or construct Q, then one needs to keep T
         (in addition to A) and pass it to further subroutines. \n
* @param[in] TSIZE
         TSIZE is INTEGER \n
         If TSIZE >= 5, the dimension of the array T. \n
         If TSIZE = -1 or -2, then a workspace query is assumed. The routine
         only calculates the sizes of the T and WORK arrays, returns these
         values as the first entries of the T and WORK arrays, and no error
         message related to T or WORK is issued by XERBLA. \n
         If TSIZE = -1, the routine calculates optimal size of T for the
         optimum performance and returns this value in T(1). \n
         If TSIZE = -2, the routine calculates minimal size of T and
         returns this value in T(1). \n
* @param[out]	WORK
         (workspace) REAL array, dimension (MAX(1,LWORK)) \n
         On exit, if INFO = 0, WORK(1) contains optimal (or either minimal
         or optimal, if query was assumed) LWORK.
         See LWORK for details. \n
* @param[in]	LWORK
         LWORK is INTEGER \n
         The dimension of the array WORK. \n
         If LWORK = -1 or -2, then a workspace query is assumed. The routine
         only calculates the sizes of the T and WORK arrays, returns these
         values as the first entries of the T and WORK arrays, and no error
         message related to T or WORK is issued by XERBLA. \n
         If LWORK = -1, the routine calculates optimal size of WORK for the
         optimal performance and returns this value in WORK(1). \n
         If LWORK = -2, the routine calculates minimal size of WORK and
         returns this value in WORK(1). \n
* @param[out]	INFO
         INFO is INTEGER \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument had an illegal value \n

* * */
    template <typename T>
    void gelq(integer *m, integer *n, T *a, integer *lda, T *t, integer *tsize, T *work,
              integer *lwork, integer *info)
    {
        gelq(m, n, a, lda, t, tsize, work, lwork, info);
    }
    /**@} */ // end of gelq

    /** @defgroup gemlq gemlq
     * @ingroup LQ
     * @{
     */
    /*! @brief Overwrites general matrix with a form compatible with orthogonal matrix
* @details
* \b Purpose:
   \verbatim
    GEMLQ overwrites the general real M-by-N matrix C with

                   SIDE = 'L'     SIDE = 'R'
   TRANS = 'N':      Q * C          C * Q
   TRANS = 'T':      Q**T * C       C * Q**T
   where Q is a real orthogonal matrix defined as the product
   of blocked elementary reflectors computed by short wide LQ
   factorization (GELQ)
   \endverbatim

* @param[in] SIDE
         SIDE is CHARACTER*1 \n
         = 'L': apply Q or Q**T from the Left; \n
         = 'R': apply Q or Q**T from the Right. \n
* @param[in] TRANS
         TRANS is CHARACTER*1 \n
         = 'N':  No transpose, apply Q; \n
         = 'T':  Transpose, apply Q**T. \n
* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix A.  M >=0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix C. N >= 0. \n
* @param[in] K
         K is INTEGER \n
         The number of elementary reflectors whose product defines
         the matrix Q. \n
         If SIDE = 'L', M >= K >= 0; \n
         if SIDE = 'R', N >= K >= 0. \n
* @param[in] A
         A is REAL array, dimension \n
                              (LDA,M) if SIDE = 'L', \n
                              (LDA,N) if SIDE = 'R' \n
         Part of the data structure to represent Q as returned by DGELQ. \n
 * @para@[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A. LDA >= fla_max(1,K). \n
* @param[in] T
         T is REAL array, dimension (MAX(5,TSIZE)). \n
         Part of the data structure to represent Q as returned by SGELQ. \n
* @param[in] TSIZE
         TSIZE is INTEGER \n
         The dimension of the array T. TSIZE >= 5. \n
* @param[in,out] C
         C is REAL array, dimension (LDC,N) \n
         On entry, the M-by-N matrix C. \n
         On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
* @param[in] LDC
         LDC is INTEGER \n
         The leading dimension of the array C. LDC >= fla_max(1,M). \n
* @param[out]	WORK
        (workspace) REAL array, dimension (MAX(1,LWORK)) \n
* @param[in]	LWORK
         LWORK is INTEGER \n
         The dimension of the array WORK. \n
         If LWORK = -1, then a workspace query is assumed. The routine
         only calculates the size of the WORK array, returns this
         value as WORK(1), and no error message related to WORK
         is issued by XERBLA. \n
* @param[out]	INFO
         INFO is INTEGER \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument had an illegal value \n
* */
    template <typename T>
    void gemlq(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *t, integer *tsize, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        gemlq(side, trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info);
    }
    /**@} */ // end of gemlq

    /** @defgroup gelqf gelqf
     * @ingroup LQ
     * @{
     */
    /*! @brief LQ factorization of a real m-by-n matrix a
   *
   * @details
   * \b Purpose:
   * \verbatim
       LQ factorization of a real m-by-n matrix a
       The factorization has the form
           A = L * Q
   \endverbatim

   * @param[in] m
             m is integer* \n
             The number of rows of the matrix a.  m >= 0. \n
   * @param[in] n
             n is integer* \n
             The number of columns of the matrix a.  n >= 0. \n
   * @param[in,out] a
             a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
             On entry, the m-by-n matrix. \n
             On exit, the elements on and below the diagonal of the array
             contain the m-by-min(m,n) lower trapezoidal matrix L (L is
             lower triangular if m <= n); the elements above the diagonal,
             with the array tau, represent the orthogonal matrix Q as a
             product of elementary reflectors (see Further Details). \n
   * @param[in] lda
             lda is integer* \n
             The leading dimension of the matrix a, lda >= fla_max(1,m) \n
   * @param[out] tau
             tau is float/double/COMPLEX/COMPLEX*16 array, dimension (min(m,n)) \n
             The scalar factors of the elementary reflectors (see Further Details).
             *
             * \n
             * **Further Details**
             * \verbatim
                     The matrix Q is represented as a product of elementary reflectors
                     Q = H(k) . . . H(2) H(1), where k = min(m,n).

                     Each H(i) has the form

                     H(i) = I - tau * V * V**T

                     where tau is a real scalar, and V is a real vector with V(1:i-1) = 0 and V(i) =
1; V(i+1:N) is stored on exit in A(i,i+1:N), and tau in tau(i). \endverbatim
   * @param[out]	WORK
             WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
             On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
   * @param[in]	LWORK
             LWORK is INTEGER \n
             The dimension of the array WORK.  LWORK >= fla_max(1,M).
             For optimum performance LWORK >= M*NB, where NB is the
             optimal blocksize. \n
\n
             If LWORK = -1, then a workspace query is assumed; the routine
             only calculates the optimal size of the WORK array, returns
             this value as the first entry of the WORK array, and no error
             message related to LWORK is issued by XERBLA. \n
   * @param[out]	INFO
             INFO is INTEGER \n
             = 0:  successful exit \n
             < 0:  if INFO = -i, the i-th argument had an illegal value \n

   *     *  */
    template <typename T>
    void gelqf(integer *m, integer *n, T *a, integer *lda, T *tau, T *work, integer *lwork,
               integer *info)
    {
        gelqf(m, n, a, lda, tau, work, lwork, info);
    }
    /**@} */ // end of gelqf

    /** @defgroup gelq2 gelq2
     * @ingroup LQ
     * @{
     */
    /*! @brief LQ factorization of a real m-by-n matrix a
   *
   * @details
   * \b Purpose:
   * \verbatim
       LQ factorization of a real m-by-n matrix a
       The factorization has the form
       A = Q * R
   \endverbatim

   * @param[in] m
             m is integer* \n
             The number of rows of the matrix a.  m >= 0. \n
   * @param[in] n
             n is integer* \n
             The number of columns of the matrix a.  n >= 0. \n
   * @param[in,out] a
             a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
             On entry, the m-by-n matrix. \n
             On exit, the elements on and below the diagonal of the array
             contain the m by min(m,n) lower trapezoidal matrix L (L is
             lower triangular if m <= n); the elements above the diagonal,
             with the array tau, represent the orthogonal matrix Q as a
             product of elementary reflectors (see Further Details). \n
   * @param[in] lda
             lda is integer* \n
             The leading dimension of the matrix a, lda >= fla_max(1,m) \n
   * @param[out] tau
             tau is float/double/COMPLEX/COMPLEX*16 array, dimension (min(m,n)) \n
             The scalar factors of the elementary reflectors (see Further
             Details).
             *
             * \n
             * **Further Details**
             * \verbatim
                     The matrix Q is represented as a product of elementary reflectors
                     Q = H(k) . . . H(2) H(1), where k = min(m,n).

                     Each H(i) has the form
                     H(i) = I - tau * V * V**T

                     where tau is a real scalar, and V is a real vector with V(1:i-1) = 0 and V(i) =
   1; V(i+1:N) is stored on exit in A(i,i+1:N), and tau in tau(i). \endverbatim
   * @param[out]	WORK
             WORK is COMPLEX array, dimension (M) \n
   * @param[out]	INFO
             INFO is INTEGER \n
             = 0: successful exit \n
             < 0: if INFO = -i, the i-th argument had an illegal value \n

   *     *  */
    template <typename T>
    void gelq2(integer *m, integer *n, T *a, integer *lda, T *tau, T *work, integer *info)
    {
        gelq2(m, n, a, lda, tau, work, info);
    }
    /**@} */ // end of gelqf

    /** @defgroup unglq unglq
     * @ingroup LQ
     * @{
     */
    /*! @brief Form Q from LQ factorization
   *
   * @details
   * \b Purpose:
   * \verbatim
       Generate an m-by-n complex matrix Q with orthonormal rows, which is defined as the first
       M rows of a product of K elementary reflectors of order N

       Q  =  H(k) . . . H(2) H(1)

       as returned by CGELQF.
   \endverbatim

   * @param[in] m
             m is integer* \n
             The number of rows of the matrix Q. m >= 0. \n
   * @param[in] n
             n is integer* \n
             The number of columns of the matrix Q. n >= m. \n
   * @param[in] k
             k is integer* \n
             The number of elementary reflectors whose product defines the
             matrix Q. m >= k >= 0. \n
   * @param[in,out] a
             a is COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
             On entry, the i-th row must contain the vector which defines
             the elementary reflector H(i), for i = 1,2,...,k, as returned
             by CGELQF in the first k rows of its array argument a. \n
             On exit, the m-by-n matrix Q. \n
   * @param[in] lda
             lda is integer* \n
             The first dimension of the array a. lda >= fla_max(1,m). \n
   * @param[in] tau
             tau is COMPLEX/COMPLEX*16 array, dimension (k) \n
             tau(i) must contain the scalar factor of the elementary
             reflector H(i), as returned by CGELQF. \n
   * @param[out]	WORK
             WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
             On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
   * @param[in]	LWORK
             LWORK is INTEGER \n
             The dimension of the array WORK. LWORK >= fla_max(1,M).
             For optimum performance LWORK >= M*NB, where NB is
             the optimal blocksize. \n
\n
             If LWORK = -1, then a workspace query is assumed; the routine
             only calculates the optimal size of the WORK array, returns
             this value as the first entry of the WORK array, and no error
             message related to LWORK is issued by XERBLA. \n
   * @param[out]	INFO
             INFO is INTEGER \n
             = 0:  successful exit; \n
             < 0:  if INFO = -i, the i-th argument has an illegal value \n

   *     *  */
    template <typename T>
    void unglq(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        unglq(m, n, k, a, lda, tau, work, lwork, info);
    }
    /**@} */ // end of unglq

    /** @defgroup orglq orglq
     * @ingroup LQ
     * @{
     */
    /*! @brief Form Q from LQ factorization
   *
   * @details
   * \b Purpose:
   * \verbatim
       Generate an m-by-n real matrix Q with orthonormal rows, which is defined as the first M
       rows of a product of K elementary reflectors of order N

       Q  =  H(k) . . . H(2) H(1)

       as returned by SGELQF.
   \endverbatim

   * @param[in] m
             m is integer* \n
             The number of rows of the matrix Q. m >= 0. \n
   * @param[in] n
             n is integer* \n
             The number of columns of the matrix Q. n >= m. \n
   * @param[in] k
             k is integer* \n
             The number of elementary reflectors whose product defines the
             matrix Q. m >= k >= 0. \n
   * @param[in,out] a
             a is float/doublearray, dimension (lda,n) \n
             On entry, the i-th row must contain the vector which defines
             the elementary reflector H(i), for i = 1,2,...,k, as returned
             by SGELQF in the first k rows of its array argument a.
             On exit, the m-by-n matrix Q. \n
   * @param[in] lda
             lda is integer* \n
             The first dimension of the array a. lda >= fla_max(1,m). \n
   * @param[in] tau
             tau is float/double array, dimension (k) \n
             tau(i) must contain the scalar factor of the elementary
             reflector H(i), as returned by SGELQF. \n
   * @param[out]	WORK
             WORK is REAL array, dimension (MAX(1,LWORK)) \n
             On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
   * @param[in]	LWORK
             LWORK is INTEGER \n
             The dimension of the array WORK. LWORK >= fla_max(1,M).
             For optimum performance LWORK >= M*NB, where NB is
             the optimal blocksize. \n
\n
             If LWORK = -1, then a workspace query is assumed; the routine
             only calculates the optimal size of the WORK array, returns
             this value as the first entry of the WORK array, and no error
             message related to LWORK is issued by XERBLA. \n
   * @param[out]	INFO
             INFO is INTEGER \n
             = 0:  successful exit \n
             < 0:  if INFO = -i, the i-th argument has an illegal value \n

   *     *  */
    template <typename T>
    void orglq(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        orglq(m, n, k, a, lda, tau, work, lwork, info);
    }
    /**@} */ // end of orglq

    /** @defgroup orgl2 orgl2
     * @ingroup LQ
     * @{
     */
    /*! @brief ORGL2 generates an m by n real matrix Q with orthonormal rows

* @details
* \b Purpose:
   \verbatim
   ORGL2 generates an m by n real matrix Q with orthonormal rows,
   which is defined as the first m rows of a product of k elementary
   reflectors of order n

         Q  =  H(k) . . . H(2) H(1)

   as   returned by GELQF.
   \endverbatim

* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix Q. M >= 0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix Q. N >= M. \n
* @param[in] K
         K is INTEGER \n
         The number of elementary reflectors whose product defines the
         matrix Q. M >= K >= 0. \n
* @param[in,out] A
         A is REAL array, dimension (LDA,N) \n
         On entry, the i-th row must contain the vector which defines
         the elementary reflector H(i), for i = 1,2,...,k, as   returned
         by SGELQF in the first k rows of its array argument A. \n
         On exit, the m-by-n matrix Q. \n
* @param[in] LDA
         LDA is INTEGER \n
         The first dimension of the array A. LDA >= fla_max(1,M). \n
* @param[in] TAU
         TAU is REAL array, dimension (K) \n
         TAU(i) must contain the scalar factor of the elementary
         reflector H(i), as   returned by SGELQF. \n
* @param[out] WORK
         WORK is REAL array, dimension (M) \n
* @param[out] INFO
         INFO is INTEGER \n
         = 0: successful exit \n
         < 0: if INFO = -i, the i-th argument has an illegal value \n

*  * */
    template <typename T>
    void orgl2(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *info)
    {
        orgl2(m, n, k, a, lda, tau, work, info);
    }
    template <typename T>
    void ungl2(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *info)
    {
        ungl2(m, n, k, a, lda, tau, work, info);
    }
    /**@} */ // end of orgl2

    /** @defgroup ormlq ormlq
     * @ingroup LQ
     * @{
     */
    /*! @brief Apply Q or Q' from LQ factorization
   *
   * @details
   * \b Purpose:
   * \verbatim
       Apply Q or Q' from LQ factorization. Overwrite the general real m-by-n matrix c with

       side = 'L'  side = 'R'
       trans = 'N':   Q * C C * Q
       trans = 'T':   Q**T * C C * Q**T

       where Q is a real orthogonal matrix defined as the product of k elementary reflectors

       Q = H(k) . . . H(2) H(1)

       as returned by SGELQF. Q is of order M if side = 'L' and of order N if side = 'R'.
   \endverbatim

   * @param[in] side
             side is char* \n
             = 'L': apply Q or Q**T from the Left; \n
             = 'R': apply Q or Q**T from the Right. \n
   * @param[in] trans
             trans is char* \n
             = 'N':  No transpose, apply Q; \n
             = 'T':  Transpose, apply Q**T. \n
   * @param[in] m
             m is integer* \n
             The number of rows of the matrix c. m >= 0. \n
   * @param[in] n
             n is integer* \n
             The number of columns of the matrix c. n >= 0. \n
   * @param[in] k
             k is integer* \n
             The number of elementary reflectors whose product defines
             the matrix Q. \n
             If side = 'L', m >= k >= 0; \n
             if side = 'R', n >= k >= 0. \n
   * @param[in] a
             a is float/double/ array, dimension \n
             (lda,m) if side = 'L', \n
             (lda,n) if side = 'R' \n
             The i-th row must contain the vector which defines the
             elementary reflector H(i), for i = 1,2,...,k, as returned by
             SGELQF in the first k rows of its array argument a. \n
   * @param[in] lda
             lda is integer* \n
             The leading dimension of the array a. lda >= fla_max(1,k). \n
   * @param[in] tau
             tau is float/double array, dimension (k) \n
             tau(i) must contain the scalar factor of the elementary
             reflector H(i), as returned by SGELQF. \n
   * @param[in,out] c
             c is float/double array, dimension (ldc,n) \n
             On entry, the m-by-n matrix c. \n
             On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
   * @param[in] ldc
             ldc is integer* \n
             The leading dimension of the array c. ldc >= fla_max(1,m). \n
   * @param[out]	WORK
             WORK is REAL array, dimension (MAX(1,LWORK)) \n
             On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
   * @param[in]	LWORK
             LWORK is INTEGER \n
             The dimension of the array WORK. \n
             If SIDE = 'L', LWORK >= fla_max(1,N); \n
             if SIDE = 'R', LWORK >= fla_max(1,M). \n
             For good performance, LWORK should generally be larger. \n
\n
             If LWORK = -1, then a workspace query is assumed; the routine
             only calculates the optimal size of the WORK array, returns
             this value as the first entry of the WORK array, and no error
             message related to LWORK is issued by XERBLA. \n
   * @param[out]	INFO
             INFO is INTEGER \n
             = 0:  successful exit \n
             < 0:  if INFO = -i, the i-th argument had an illegal value \n

   *     *  */
    template <typename T>
    void ormlq(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        ormlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
    /**@} */ // end of ormlq

    /** @defgroup unmlq unmlq
     * @ingroup LQ
     * @{
     */
    /*! @brief Apply Q or Q' from LQ factorization
   *
   * @details
   * \b Purpose:
   * \verbatim
       Apply Q or Q' from LQ factorization. Overwrite the general complex m-by-n matrix c with

       side = 'L'  side = 'R'
       trans = 'N':   Q * C C * Q
       trans = 'C':   Q**H * C C * Q**H

       where Q is a complex unitary matrix defined as the product of k elementary reflectors

       Q = H(k)**H . . . H(2)**H H(1)**H

       as returned by CGELQF. Q is of order M if side = 'L' and of order N if side = 'R'.
   \endverbatim

   * @param[in] side
             side is char* \n
             = 'L': apply Q or Q**H from the Left; \n
             = 'R': apply Q or Q**H from the Right. \n
   * @param[in] trans
             trans is char* \n
             = 'N':  No transpose, apply Q; \n
             = 'C':  Conjugate transpose, apply Q**H. \n
   * @param[in] m
             m is integer* \n
             The number of rows of the matrix c. m >= 0. \n
   * @param[in] n
             n is integer* \n
             The number of columns of the matrix c. n >= 0. \n
   * @param[in] k
             k is integer* \n
             The number of elementary reflectors whose product defines
             the matrix Q. \n
             If side = 'L', m >= k >= 0; \n
             if side = 'R', n >= k >= 0. \n
   * @param[in] a
             a is COMPLEX/COMPLEX*16 array, dimension \n
             (lda,m) if side = 'L', \n
             (lda,n) if side = 'R' \n
             The i-th row must contain the vector which defines the
             elementary reflector H(i), for i = 1,2,...,k, as returned by
             CGELQF in the first k rows of its array argument a. \n
   * @param[in] lda
             lda is integer* \n
             The leading dimension of the array a. lda >= fla_max(1,k). \n
   * @param[in] tau
             tau is COMPLEX/COMPLEX*16 array, dimension (k) \n
             tau(i) must contain the scalar factor of the elementary
             reflector H(i), as returned by CGELQF. \n
   * @param[in,out] c
             c is COMPLEX/COMPLEX*16 array, dimension (ldc,n) \n
             On entry, the m-by-n matrix c. \n
             On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
   * @param[in] ldc
             ldc is integer* \n
             The leading dimension of the array c. ldc >= fla_max(1,m). \n
   * @param[out]	WORK
             WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
             On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
   * @param[in]	LWORK
             LWORK is INTEGER \n
             The dimension of the array WORK. \n
             If SIDE = 'L', LWORK >= fla_max(1,N); \n
             if SIDE = 'R', LWORK >= fla_max(1,M). \n
             For good performance, LWORK should generally be larger. \n
\n
             If LWORK = -1, then a workspace query is assumed; the routine
             only calculates the optimal size of the WORK array, returns
             this value as the first entry of the WORK array, and no error
             message related to LWORK is issued by XERBLA. \n
   * @param[out]	INFO
             INFO is INTEGER \n
             = 0:  successful exit \n
             < 0:  if INFO = -i, the i-th argument had an illegal value \n

   *     *  */
    template <typename T>
    void unmlq(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        unmlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
    /**@} */ // end of unmlq

    /** @defgroup orml2 orml2
     * @ingroup LQ
     * @{
     */
    /*! @brief Apply Q or Q' from LQ factorization
   *
   * @details
   * \b Purpose:
   * \verbatim
       Apply Q or Q' from LQ factorization. Overwrite the general real m-by-n matrix c with

       Q * C  if side = 'L' and trans = 'N', or

       Q**T* C  if side = 'L' and trans = 'T', or

       C * Q  if side = 'R' and trans = 'N', or

       C * Q**T if side = 'R' and trans = 'T',

       where Q is a real orthogonal matrix defined as the product of k elementary reflectors

       Q = H(k) . . . H(2) H(1)

       as returned by SGELQF. Q is of order m if side = 'L' and of order n
   \endverbatim

   * @param[in] side
             side is char* \n
             = 'L': apply Q or Q**T from the Left; \n
             = 'R': apply Q or Q**T from the Right. \n
   * @param[in] trans
             trans is char* \n
             = 'N':  No transpose, apply Q; \n
             = 'T':  Transpose, apply Q**T. \n
   * @param[in] m
             m is integer* \n
             The number of rows of the matrix c. m >= 0. \n
   * @param[in] n
             n is integer* \n
             The number of columns of the matrix c. n >= 0. \n
   * @param[in] k
             k is integer* \n
             The number of elementary reflectors whose product defines
             the matrix Q. \n
             If side = 'L', m >= k >= 0; \n
             if side = 'R', n >= k >= 0. \n
   * @param[in] a
             a is float/double array, dimension \n
             (lda,m) if side = 'L', \n
             (lda,n) if side = 'R' \n
             The i-th row must contain the vector which defines the
             elementary reflector H(i), for i = 1,2,...,k, as returned by
             SGELQF in the first k rows of its array argument a. \n
             a is modified by the routine but restored on exit. \n
   * @param[in] lda
             lda is integer* \n
             The leading dimension of the array a. lda >= fla_max(1,k). \n
   * @param[in] tau
             tau is float/double array, dimension (k) \n
             tau(i) must contain the scalar factor of the elementary
             reflector H(i), as returned by SGELQF. \n
   * @param[in,out] c
             c is float/double array, dimension (ldc,n) \n
             On entry, the m-by-n matrix c. \n
             On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
   * @param[in] ldc
             ldc is integer* \n
             The leading dimension of the array c. ldc >= fla_max(1,m). \n
   * @param[out]	WORK
             WORK is REAL array, dimension \n
                                      (N) if SIDE = 'L', \n
                                      (M) if SIDE = 'R' \n
   * @param[out]	INFO
             INFO is INTEGER \n
             = 0: successful exit \n
             < 0: if INFO = -i, the i-th argument had an illegal value \n

   *     *  */
    template <typename T>
    void orml2(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *info)
    {
        orml2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    }
    /**@} */ // end of orml2

    /** @defgroup unml2 unml2
     * @ingroup LQ
     * @{
     */
    /*! @brief Apply Q or Q' from LQ factorization
   *
   * @details
   * \b Purpose:
   * \verbatim
       Apply Q or Q' from LQ factorization. Overwrite the general complex m-by-n matrix c with

       Q * C  if side = 'L' and trans = 'N', or

       Q**H* C  if side = 'L' and trans = 'C', or

       C * Q  if side = 'R' and trans = 'N', or

       C * Q**H if side = 'R' and trans = 'C',

       where Q is a complex unitary matrix defined as the product of k elementary reflectors

       Q = H(k)**H . . . H(2)**H H(1)**H

       as returned by CGELQF. Q is of order m if side = 'L' and of order n
   \endverbatim

   * @param[in] side
             side is char* \n
             = 'L': apply Q or Q**T from the Left; \n
             = 'R': apply Q or Q**T from the Right. \n
   * @param[in] trans
             trans is char* \n
             = 'N':  No transpose, apply Q; \n
             = 'T':  Transpose, apply Q**T. \n
   * @param[in] m
             m is integer* \n
             The number of rows of the matrix c. m >= 0. \n
   * @param[in] n
             n is integer* \n
             The number of columns of the matrix c. n >= 0. \n
   * @param[in] k
             k is integer* \n
             The number of elementary reflectors whose product defines
             the matrix Q. \n
             If side = 'L', m >= k >= 0; \n
             if side = 'R', n >= k >= 0. \n
   * @param[in] a
             a is COMPLEX/COMPLEX*16 array, dimension \n
             (lda,m) if side = 'L', \n
             (lda,n) if side = 'R' \n
             The i-th row must contain the vector which defines the
             elementary reflector H(i), for i = 1,2,...,k, as returned by
             CGELQF in the first k rows of its array argument a. \n
             a is modified by the routine but restored on exit. \n
   * @param[in] lda
             lda is integer* \n
             The leading dimension of the array a. lda >= fla_max(1,k). \n
   * @param[in] tau
             tau is COMPLEX/COMPLEX*16 array, dimension (k) \n
             tau(i) must contain the scalar factor of the elementary
             reflector H(i), as returned by CGELQF. \n
   * @param[in,out] c
             c is COMPLEX/COMPLEX*16 array, dimension (ldc,n) \n
             On entry, the m-by-n matrix c. \n
             On exit, c is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. \n
   * @param[in] ldc
             ldc is integer* \n
             The leading dimension of the array c. ldc >= fla_max(1,m). \n
   * @param[out]	WORK
             WORK is COMPLEX array, dimension \n
                                      (N) if SIDE = 'L', \n
                                      (M) if SIDE = 'R' \n
   * @param[out]	INFO
             INFO is INTEGER \n
             = 0: successful exit \n
             < 0: if INFO = -i, the i-th argument had an illegal value \n

   *     *  */
    template <typename T>
    void unml2(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *info)
    {
        unml2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    }
    /**@} */ // end of unml2

    /** @defgroup gelqt gelqt
     * @ingroup LQ
     * @{
     */
    /*! @brief GELQT computes a blocked LQ factorization of a real M-by-N matrix A
   using the compact WY representation of Q.

* @details
* \b Purpose:
   \verbatim
   GELQT computes a blocked LQ factorization of a real M-by-N matrix A
   using the compact WY representation of Q.
   \endverbatim

* @param[in]	M
         M is INTEGER \n
         The number of rows of the matrix A.  M >= 0. \n
* @param[in]	N
         N is INTEGER \n
         The number of columns of the matrix A.  N >= 0. \n
* @param[in]	MB
         MB is INTEGER \n
         The block size to be used in the blocked QR.  MIN(M,N) >= MB >= 1. \n
* @param[in,out]	A
         A is REAL array, dimension (LDA,N) \n
         On entry, the M-by-N matrix A. \n
         On exit, the elements on and below the diagonal of the array
         contain the M-by-MIN(M,N) lower trapezoidal matrix L (L is
         lower triangular if M <= N); the elements above the diagonal
         are the rows of V. \n
* @param[in]	LDA
         LDA is INTEGER \n
         The leading dimension of the array A.  LDA >= fla_max(1,M). \n
* @param[out]	T
         T is REAL array, dimension (LDT,MIN(M,N)) \n
         The upper triangular block reflectors stored in compact form
         as a sequence of upper triangular blocks.  See below
         for further details. \n
* @param[in]	LDT
         LDT is INTEGER \n
         The leading dimension of the array T.  LDT >= MB. \n
* @param[out]	WORK
         WORK is REAL array, dimension (MB*N) \n
* @param[out]	INFO
         INFO is INTEGER \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void gelqt(integer *m, integer *n, integer *mb, T *a, integer *lda, T *t, integer *ldt, T *work,
               integer *info)
    {
        gelqt(m, n, mb, a, lda, t, ldt, work, info);
    }
    /**@} */ // end of gelqt

    /** @defgroup gelqt3 gelqt3
     * @ingroup LQ
     * @{
     */
    /*! @brief GELQT3 recursively computes a LQ factorization of a real M-by-N matrix A

* @details
* \b Purpose:
   \verbatim
   GELQT3 recursively computes a LQ factorization of a real M-by-N
   matrix A, using the compact WY representation of Q.

   Based on the algorithm of Elmroth and Gustavson,
   IBM J. Res. Develop. Vol 44 No. 4 July 2000.
   \endverbatim

* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix A.  M =< N. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix A.  N >= 0. \n
* @param[in,out] A
         A is REAL array, dimension (LDA,N) \n
         On entry, the real M-by-N matrix A.  On exit, the elements on and
         below the diagonal contain the N-by-N lower triangular matrix L; the
         elements above the diagonal are the rows of V.  See below for
         further details. \n
* @param[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A.  LDA >= fla_max(1,M). \n
* @param[out] T
         T is REAL array, dimension (LDT,N) \n
         The N-by-N upper triangular factor of the block reflector.
         The elements on and above the diagonal contain the block
         reflector T; the elements below the diagonal are not used.
         See below for further details. \n
* @param[in] LDT
         LDT is INTEGER \n
         The leading dimension of the array T.  LDT >= fla_max(1,N). \n
* @param[out] INFO
         INFO is INTEGER \n
         = 0: successful exit \n
         < 0: if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void gelqt3(integer *m, integer *n, T *a, integer *lda, T *t, integer *ldt, integer *info)
    {
        gelqt3(m, n, a, lda, t, ldt, info);
    }
    /**@} */ // end of gelqt3

    /** @defgroup gemlqr gemlqt
     * @ingroup LQ
     * @{
     */
    /*! @brief GEMLQT overwrites the general real M-by-N matrix C

* @details
* \b Purpose:
   \verbatim
   GEMLQT overwrites the general real M-by-N matrix C with

                  SIDE = 'L'     SIDE = 'R'
   TRANS = 'N':      Q C            C Q
   TRANS = 'T':   Q**T C            C Q**T

   where Q is a real orthogonal matrix defined as the product of K
   elementary reflectors:

        Q = H(1) H(2) . . . H(K) = I - V T V**T

   generated using the compact WY representation as returned by DGELQT.

   Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.
   \endverbatim

* @param[in]	SIDE
         SIDE is CHARACTER*1 \n
         = 'L': apply Q or Q**T from the Left; \n
         = 'R': apply Q or Q**T from the Right. \n
* @param[in]	TRANS
         TRANS is CHARACTER*1 \n
         = 'N':  No transpose, apply Q; \n
         = 'C':  Transpose, apply Q**T. \n
* @param[in]	M
         M is INTEGER \n
         The number of rows of the matrix C. M >= 0. \n
* @param[in]	N
         N is INTEGER \n
         The number of columns of the matrix C. N >= 0. \n
* @param[in]	K
         K is INTEGER \n
         The number of elementary reflectors whose product defines
         the matrix Q. \n
         If SIDE = 'L', M >= K >= 0; \n
         if SIDE = 'R', N >= K >= 0. \n
* @param[in]	MB
         MB is INTEGER \n
         The block size used for the storage of T.  K >= MB >= 1. \n
         This must be the same value of MB used to generate T
         in DGELQT. \n
* @param[in]	V
         V is REAL array, dimension \n
                              (LDV,M) if SIDE = 'L', \n
                              (LDV,N) if SIDE = 'R' \n
         The i-th row must contain the vector which defines the
         elementary reflector H(i), for i = 1,2,...,k, as returned by
         DGELQT in the first K rows of its array argument A. \n
* @param[in]	LDV
         LDV is INTEGER \n
         The leading dimension of the array V. LDV >= fla_max(1,K). \n
* @param[in]	T
         T is REAL array, dimension (LDT,K) \n
         The upper triangular factors of the block reflectors
         as returned by DGELQT, stored as a MB-by-K matrix. \n
* @param[in]	LDT
         LDT is INTEGER \n
         The leading dimension of the array T.  LDT >= MB. \n
* @param[in,out]	C
         C is REAL array, dimension (LDC,N) \n
         On entry, the M-by-N matrix C. \n
         On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q. \n
* @param[in]	LDC
         LDC is INTEGER \n
         The leading dimension of the array C. LDC >= fla_max(1,M). \n
* @param[out]	WORK
         WORK is REAL array. The dimension of \n
         WORK is N*MB if SIDE = 'L', or  M*MB if SIDE = 'R'. \n
* @param[out]	INFO
         INFO is INTEGER \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void gemlqt(char *side, char *trans, integer *m, integer *n, integer *k, integer *mb, T *v,
                integer *ldv, T *t, integer *ldt, T *c, integer *ldc, T *work, integer *info)
    {
        gemlqt(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info);
    }
    /**@} */ // end of gemlqt

    /** @defgroup lamswlq lamswql
     * @ingroup LQ
     * @{
     */
    /*! @brief LAMSWLQ overwrites the general real M-by-N matrix C

* @details
* \b Purpose:
   \verbatim
   LAMSWLQ overwrites the general real M-by-N matrix C with

                   SIDE = 'L'     SIDE = 'R'
   TRANS = 'N':      Q * C          C * Q
   TRANS = 'T':      Q**T * C       C * Q**T
   where Q is a real orthogonal matrix defined as the product of blocked
   elementary reflectors computed by short wide LQ
   factorization (LASWLQ)
   \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q**T from the Left; \n
          = 'R': apply Q or Q**T from the Right. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N':  No transpose, apply Q; \n
          = 'T':  Transpose, apply Q**T. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C.  M >=0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. N >= M. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines
          the matrix Q.
          M >= K >= 0; \n
 * @param[in] MB
          MB is INTEGER \n
          The row block size to be used in the blocked QR.
          M >= MB >= 1 \n
 * @param[in] NB
          NB is INTEGER \n
          The column block size to be used in the blocked QR.
          NB > M. \n
 * @param[in] NB
          NB is INTEGER \n
          The block size to be used in the blocked QR.
                MB > M. \n
 * @param[in] A
          A is REAL array, dimension \n
                               (LDA,M) if SIDE = 'L', \n
                               (LDA,N) if SIDE = 'R' \n
          The i-th row must contain the vector which defines the blocked
          elementary reflector H(i), for i = 1,2,...,k, as   returned by
          SLASWLQ in the first k rows of its array argument A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. \n
          If SIDE = 'L', LDA >= fla_max(1,M); \n
          if SIDE = 'R', LDA >= fla_max(1,N). \n
 * @param[in] T
          T is REAL array, dimension \n
          (M * Number of blocks(CEIL(N-K/NB-K)), \n
          The blocked upper triangular block reflectors stored in compact form
          as a sequence of upper triangular blocks.  See below
          for further details. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= MB. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the M-by-N matrix C.
          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out] WORK
         (workspace) REAL array, dimension (MAX(1,LWORK)) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If SIDE = 'L', LWORK >= fla_max(1,NB) * MB; \n
          if SIDE = 'R', LWORK >= fla_max(1,M) * MB. \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array,   returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

* */
    template <typename T>
    void lamswlq(char *side, char *trans, integer *m, integer *n, integer *k, integer *mb,
                 integer *nb, T *a, integer *lda, T *t, integer *ldt, T *c, integer *ldc, T *work,
                 integer *lwork, integer *info)
    {
        lamswlq(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
    }
    /**@} */ // end of lamswql

    /** @defgroup tplqt tplqt
     * @ingroup LQ
     * @{
     */
    /*! @brief TPLQT computes a blocked LQ factorization of a real "triangular-pentagonal" matrix C

* @details
* \b Purpose:
   \verbatim
   TPLQT computes a blocked LQ factorization of a real
   "triangular-pentagonal" matrix C, which is composed of a
   triangular block A and pentagonal block B, using the compact
   WY representation for Q.
   \endverbatim

* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix B, and the order of the
         triangular matrix A. \n
         M >= 0.
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix B.
         N >= 0. \n
* @param[in] L
         L is INTEGER \n
         The number of rows of the lower trapezoidal part of B.
         MIN(M,N) >= L >= 0.  See Further Details. \n
* @param[in] MB
         MB is INTEGER \n
         The block size to be used in the blocked QR.  M >= MB >= 1. \n
* @param[in,out] A
         A is REAL array, dimension (LDA,M) \n
         On entry, the lower triangular M-by-M matrix A.
         On exit, the elements on and below the diagonal of the array
         contain the lower triangular matrix L. \n
* @param[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A.  LDA >= fla_max(1,M). \n
* @param[in,out] B
         B is REAL array, dimension (LDB,N) \n
         On entry, the pentagonal M-by-N matrix B.  The first N-L columns
         are rectangular, and the last L columns are lower trapezoidal.
         On exit, B contains the pentagonal matrix V.  See Further Details. \n
* @param[in] LDB
         LDB is INTEGER \n
         The leading dimension of the array B.  LDB >= fla_max(1,M). \n
* @param[out] T
         T is REAL array, dimension (LDT,N) \n
         The lower triangular block reflectors stored in compact form
         as a sequence of upper triangular blocks.  See Further Details. \n
* @param[in] LDT
         LDT is INTEGER \n
         The leading dimension of the array T.  LDT >= MB. \n
* @param[out] WORK
         WORK is REAL array, dimension (MB*M) \n
* @param[out] INFO
         INFO is INTEGER \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument had an illegal value  \n

*  * */
    template <typename T>
    void tplqt(integer *m, integer *n, integer *l, integer *mb, T *a, integer *lda, T *b,
               integer *ldb, T *t, integer *ldt, T *work, integer *info)
    {
        tplqt(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info);
    }
    /**@} */ // end of tplqt

    /** @defgroup tplqt2 tplqt2
     * @ingroup LQ
     * @{
     */
    /*! @brief TPLQT2 computes a LQ factorization of a real or complex "triangular-pentagonal"
matrix

* @details
* \b Purpose:
   \verbatim
   TPLQT2 computes a LQ a factorization of a real "triangular-pentagonal"
   matrix C, which is composed of a triangular block A and pentagonal block B,
   using the compact WY representation for Q.
   \endverbatim

* @param[in] M
         M is INTEGER \n
         The total number of rows of the matrix B.
         M >= 0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix B, and the order of
         the triangular matrix A.
         N >= 0. \n
* @param[in] L
         L is INTEGER \n
         The number of rows of the lower trapezoidal part of B.
         MIN(M,N) >= L >= 0.  See Further Details. \n
* @param[in,out] A
         A is REAL array, dimension (LDA,M) \n
         On entry, the lower triangular M-by-M matrix A. \n
         On exit, the elements on and below the diagonal of the array
         contain the lower triangular matrix L. \n
* @param[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A.  LDA >= fla_max(1,M). \n
* @param[in,out] B
         B is REAL array, dimension (LDB,N) \n
         On entry, the pentagonal M-by-N matrix B.  The first N-L columns
         are rectangular, and the last L columns are lower trapezoidal. \n
         On exit, B contains the pentagonal matrix V.  See Further Details. \n
* @param[in] LDB
         LDB is INTEGER \n
         The leading dimension of the array B.  LDB >= fla_max(1,M). \n
* @param[out] T
         T is REAL array, dimension (LDT,M) \n
         The N-by-N upper triangular factor T of the block reflector.
         See Further Details. \n
* @param[in] LDT
         LDT is INTEGER \n
         The leading dimension of the array T.  LDT >= fla_max(1,M) \n
* @param[out] INFO
         INFO is INTEGER \n
         = 0: successful exit \n
         < 0: if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void tplqt2(integer *m, integer *n, integer *l, T *a, integer *lda, T *b, integer *ldb, T *t,
                integer *ldt, integer *info)
    {
        tplqt2(m, n, l, a, lda, b, ldb, t, ldt, info);
    }
    /**@} */ // end of tplqt2

    /** @defgroup tpmlqt tpmlqt
     * @ingroup LQ
     * @{
     */
    /*! @brief TPMLQT applies a real orthogonal matrix Q obtained from a "triangular-pentagonal"
real block reflector

* @details
* \b Purpose:
   \verbatim
   TPMLQT applies a real orthogonal matrix Q obtained from a
   "triangular-pentagonal" real block reflector H to a general
   real matrix C, which consists of two blocks A and B.
   \endverbatim

* @param[in] SIDE
         SIDE is CHARACTER*1 \n
         = 'L': apply Q or Q**T from the Left; \n
         = 'R': apply Q or Q**T from the Right. \n
* @param[in] TRANS
         TRANS is CHARACTER*1 \n
         = 'N':  No transpose, apply Q; \n
         = 'T':  Transpose, apply Q**T. \n
* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix B. M >= 0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix B. N >= 0. \n
* @param[in] K
         K is INTEGER \n
         The number of elementary reflectors whose product defines
         the matrix Q. \n
* @param[in] L
         L is INTEGER \n
         The order of the trapezoidal part of V.
         K >= L >= 0.  See Further Details. \n
* @param[in] MB
         MB is INTEGER \n
         The block size used for the storage of T.  K >= MB >= 1.
         This must be the same value of MB used to generate T
         in DTPLQT. \n
* @param[in] V
         V is REAL array, dimension (LDV,K) \n
         The i-th row must contain the vector which defines the
         elementary reflector H(i), for i = 1,2,...,k, as   returned by
         DTPLQT in B.  See Further Details. \n
* @param[in] LDV
         LDV is INTEGER \n
         The leading dimension of the array V. \n
         If SIDE = 'L', LDV >= fla_max(1,M); \n
         if SIDE = 'R', LDV >= fla_max(1,N). \n
* @param[in] T
         T is REAL array, dimension (LDT,K) \n
         The upper triangular factors of the block reflectors
         as   returned by DTPLQT, stored as a MB-by-K matrix. \n
* @param[in] LDT
         LDT is INTEGER \n
         The leading dimension of the array T.  LDT >= MB. \n
* @param[in,out] A
         A is REAL array, dimension
         (LDA,N) if SIDE = 'L' or
         (LDA,K) if SIDE = 'R' \n
         On entry, the K-by-N or M-by-K matrix A. \n
         On exit, A is overwritten by the corresponding block of
         Q*C or Q**T*C or C*Q or C*Q**T.  See Further Details. \n
* @param[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A. \n
         If SIDE = 'L', LDC >= fla_max(1,K); \n
         If SIDE = 'R', LDC >= fla_max(1,M). \n
* @param[in,out] B
         B is REAL array, dimension (LDB,N) \n
         On entry, the M-by-N matrix B. \n
         On exit, B is overwritten by the corresponding block of
         Q*C or Q**T*C or C*Q or C*Q**T.  See Further Details. \n
* @param[in] LDB
         LDB is INTEGER \n
         The leading dimension of the array B.
         LDB >= fla_max(1,M). \n
* @param[out] WORK
         WORK is REAL array. The dimension of WORK is
          N*MB if SIDE = 'L', or  M*MB if SIDE = 'R'. \n
* @param[out] INFO
         INFO is INTEGER \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void tpmlqt(char *side, char *trans, integer *m, integer *n, integer *k, integer *l,
                integer *mb, T *v, integer *ldv, T *t, integer *ldt, T *a, integer *lda, T *b,
                integer *ldb, T *work, integer *info)
    {
        tpmlqt(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
    }
    /**@} */ // end of tpmlqt

    /** @defgroup geqlf geqlf
     * @ingroup LQ
     * @{
     */
    /*! @brief GEQLF computes a QL factorization of M-by-N matrix
* @details
* \b Purpose:
   \verbatim
   GEQLF computes a QL factorization of M-by-N matrix A:
   A = Q * L.
   \endverbatim

* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix A.  M >= 0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix A.  N >= 0. \n
* @param[in,out] A
         A is REAL array, dimension (LDA,N) \n
         On entry, the M-by-N matrix A. \n
         On exit,
         if m >= n, the lower triangle of the subarray
         A(m-n+1:m,1:n) contains the N-by-N lower triangular matrix L;
         if m <= n, the elements on and below the (n-m)-th
         superdiagonal contain the M-by-N lower trapezoidal matrix L;
         the remaining elements, with the array TAU, represent the
         orthogonal matrix Q as a product of elementary reflectors
         (see Further Details). \n
* @param[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A.  LDA >= fla_max(1,M). \n
* @param[out] TAU
         TAU is array, dimension (min(M,N)) \n
* @param[out]	WORK
         WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
         On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
* @param[in]	LWORK
         LWORK is INTEGER \n
         The dimension of the array WORK.  LWORK >= fla_max(1,N).
         For optimum performance LWORK >= N*NB, where NB is
         the optimal blocksize. \n
\n
         If LWORK = -1, then a workspace query is assumed; the routine
         only calculates the optimal size of the WORK array, returns
         this value as the first entry of the WORK array, and no error
         message related to LWORK is issued by XERBLA. \n
* @param[out]	INFO
         INFO is INTEGER \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void geqlf(integer *m, integer *n, T *a, integer *lda, T *tau, T *work, integer *lwork,
               integer *info)
    {
        geqlf(m, n, a, lda, tau, work, lwork, info);
    }
    /**@} */ // end of  geqlf

    /** @defgroup gerq2 gerq2
     * @ingroup LQ
     * @{
     */
    /*! @brief GERQ2 computes the RQ factorization of a general rectangular matrix using an
unblocked algorithm

* @details
* \b Purpose:
   \verbatim
   GERQ2 computes an RQ factorization of a real m by n matrix A:
   A = R * Q.
   \endverbatim

* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix A.  M >= 0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix A.  N >= 0. \n
* @param[in,out] A
         A is REAL array, dimension (LDA,N) \n
         On entry, the m by n matrix A. \n
         On exit, if m <= n, the upper triangle of the subarray
         A(1:m,n-m+1:n) contains the m by m upper triangular matrix R;
         if m >= n, the elements on and above the (m-n)-th subdiagonal
         contain the m by n upper trapezoidal matrix R; the remaining
         elements, with the array TAU, represent the orthogonal matrix
         Q as a product of elementary reflectors (see Further
         Details). \n
* @param[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A.  LDA >= fla_max(1,M). \n
* @param[out] TAU
         TAU is REAL array, dimension (min(M,N)) \n
         The scalar factors of the elementary reflectors (see Further
         Details). \n
* @param[out] WORK
         WORK is REAL array, dimension (M) \n
* @param[out] INFO
         INFO is INTEGER \n
         = 0: successful exit \n
         < 0: if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void gerq2(integer *m, integer *n, T *a, integer *lda, T *tau, T *work, integer *info)
    {
        gerq2(m, n, a, lda, tau, work, info);
    }
    /**@} */ // end of gerq2

    /** @defgroup ungrq ungrq
     * @ingroup LQ
     * @{
     */
    /*! @brief ORGRQ generates an M-by-N real matrix Q with orthonormal rows

* @details
* \b Purpose:
   \verbatim
    ORGRQ generates an M-by-N real matrix Q with orthonormal rows,
    which is defined as the last M rows of a product of K elementary
    reflectors of order N

          Q  =  H(1) H(2) . . . H(k)

    as returned by SGERQF.
   \endverbatim

* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix Q. M >= 0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix Q. N >= M. \n
* @param[in] K
         K is INTEGER \n
         The number of elementary reflectors whose product defines the
         matrix Q. M >= K >= 0. \n
* @param[in,out] A
         A is REAL array, dimension (LDA,N) \n
         On entry, the (m-k+i)-th row must contain the vector which
         defines the elementary reflector H(i), for i = 1,2,...,k, as
         returned by SGERQF in the last k rows of its array argument
         A. \n
         On exit, the M-by-N matrix Q. \n
* @param[in] LDA
         LDA is INTEGER \n
         The first dimension of the array A. LDA >= fla_max(1,M). \n
* @param[in] TAU
         TAU is REAL array, dimension (K) \n
         TAU(i) must contain the scalar factor of the elementary
         reflector H(i), as returned by SGERQF. \n
* @param[out]	WORK
         WORK is REAL array, dimension (MAX(1,LWORK)) \n
         On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
* @param[in]	LWORK
         LWORK is INTEGER \n
         The dimension of the array WORK. LWORK >= fla_max(1,M).
         For optimum performance LWORK >= M*NB, where NB is the
         optimal blocksize. \n
\n
         If LWORK = -1, then a workspace query is assumed; the routine
         only calculates the optimal size of the WORK array, returns
         this value as the first entry of the WORK array, and no error
         message related to LWORK is issued by XERBLA. \n
* @param[out]	INFO
         INFO is INTEGER \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument has an illegal value \n

*  * */
    template <typename T>
    void orgrq(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        orgrq(m, n, k, a, lda, tau, work, lwork, info);
    }
    template <typename T>
    void ungrq(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        ungrq(m, n, k, a, lda, tau, work, lwork, info);
    }
    /**@} */ // end of ungrq

    /** @defgroup unmrq unmrq
     * @ingroup LQ
     * @{
     */
    /*! @brief ORMRQ overwrites the general real M-by-N matrix

* @details
* \b Purpose:
   \verbatim
    ORMRQ overwrites the general real M-by-N matrix C with

                    SIDE = 'L'     SIDE = 'R'
    TRANS = 'N':      Q * C          C * Q
    TRANS = 'T':      Q**T * C       C * Q**T

    where Q is a real orthogonal matrix defined as the product of k
    elementary reflectors

          Q = H(1) H(2) . . . H(k)

    as returned by SGERQF. Q is of order M if SIDE = 'L' and of order N
    if SIDE = 'R'.
   \endverbatim

* @param[in] SIDE
         SIDE is CHARACTER*1 \n
         = 'L': apply Q or Q**T from the Left; \n
         = 'R': apply Q or Q**T from the Right. \n
* @param[in] TRANS
         TRANS is CHARACTER*1 \n
         = 'N':  No transpose, apply Q; \n
         = 'T':  Transpose, apply Q**T. \n
* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix C. M >= 0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix C. N >= 0. \n
* @param[in] K
         K is INTEGER \n
         The number of elementary reflectors whose product defines
         the matrix Q. \n
         If SIDE = 'L', M >= K >= 0; \n
         if SIDE = 'R', N >= K >= 0. \n
* @param[in] A
         A is REAL array, dimension \n
                              (LDA,M) if SIDE = 'L', \n
                              (LDA,N) if SIDE = 'R' \n
         The i-th row must contain the vector which defines the
         elementary reflector H(i), for i = 1,2,...,k, as returned by
         SGERQF in the last k rows of its array argument A. \n
* @param[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A. LDA >= fla_max(1,K). \n
* @param[in] TAU
         TAU is REAL array, dimension (K) \n
         TAU(i) must contain the scalar factor of the elementary
         reflector H(i), as returned by SGERQF. \n
* @param[in,out] C
         C is REAL array, dimension (LDC,N) \n
         On entry, the M-by-N matrix C. \n
         On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
* @param[in] LDC
         LDC is INTEGER \n
         The leading dimension of the array C. LDC >= fla_max(1,M). \n
* @param[out]	WORK
         WORK is REAL array, dimension (MAX(1,LWORK)) \n
         On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
* @param[in]	LWORK
         LWORK is INTEGER \n
         The dimension of the array WORK. \n
         If SIDE = 'L', LWORK >= fla_max(1,N); \n
         if SIDE = 'R', LWORK >= fla_max(1,M). \n
         For good performance, LWORK should generally be larger. \n
\n
         If LWORK = -1, then a workspace query is assumed; the routine
         only calculates the optimal size of the WORK array, returns
         this value as the first entry of the WORK array, and no error
         message related to LWORK is issued by XERBLA. \n
* @param[out]	INFO
         INFO is INTEGER \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void ormrq(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        ormrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
    template <typename T>
    void unmrq(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        unmrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
    /**@} */ // end of unmrq

    /** @defgroup ormr2 ormr2
     * @ingroup LQ
     * @{
     */
    /*! @brief ORMR2 multiplies a general matrix by the orthogonal matrix from a RQ factorization
determined by sgerqf

* @details
* \b Purpose:
   \verbatim
   ORMR2 overwrites the general real m by n matrix C with

         Q * C  if SIDE = 'L' and TRANS = 'N', or

         Q**T* C  if SIDE = 'L' and TRANS = 'T', or

         C * Q  if SIDE = 'R' and TRANS = 'N', or

         C * Q**T if SIDE = 'R' and TRANS = 'T',

   where Q is a real orthogonal matrix defined as the product of k
   elementary reflectors

         Q = H(1) H(2) . . . H(k)

   as   returned by SGERQF. Q is of order m if SIDE = 'L' and of order n
   if SIDE = 'R'.
   \endverbatim

* @param[in] SIDE
         SIDE is CHARACTER*1 \n
         = 'L': apply Q or Q**T from the Left \n
         = 'R': apply Q or Q**T from the Right \n
* @param[in] TRANS
         TRANS is CHARACTER*1 \n
         = 'N': apply Q  (No transpose) \n
         = 'T': apply Q' (Transpose) \n
* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix C. M >= 0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix C. N >= 0. \n
* @param[in] K
         K is INTEGER \n
         The number of elementary reflectors whose product defines
         the matrix Q. \n
         If SIDE = 'L', M >= K >= 0; \n
         if SIDE = 'R', N >= K >= 0. \n
* @param[in] A
         A is REAL array, dimension \n
                              (LDA,M) if SIDE = 'L', \n
                              (LDA,N) if SIDE = 'R' \n
         The i-th row must contain the vector which defines the
         elementary reflector H(i), for i = 1,2,...,k, as   returned by
         SGERQF in the last k rows of its array argument A.
         A is modified by the routine but restored on exit. \n
* @param[in] LDA
         LDA is INTEGER \n
         The leading dimension of the array A. LDA >= fla_max(1,K). \n
* @param[in] TAU
         TAU is REAL array, dimension (K) \n
         TAU(i) must contain the scalar factor of the elementary
         reflector H(i), as   returned by SGERQF. \n
* @param[in,out] C
         C is REAL array, dimension (LDC,N) \n
         On entry, the m by n matrix C. \n
         On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
* @param[in] LDC
         LDC is INTEGER \n
         The leading dimension of the array C. LDC >= fla_max(1,M). \n
* @param[out] WORK
         WORK is REAL array, dimension \n
                                  (N) if SIDE = 'L', \n
                                  (M) if SIDE = 'R' \n
* @param[out] INFO
         INFO is INTEGER \n
         = 0: successful exit \n
         < 0: if INFO = -i, the i-th argument had an illegal value  \n

*  * */
    template <typename T>
    void ormr2(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *info)
    {
        ormr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    }
    template <typename T>
    void unmr2(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *info)
    {
        unmr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    }
    /**@} */ // end of ormr2

    /** @defgroup orgr2 orgr2
     * @ingroup LQ
     * @{
     */
    /*! @brief ORGR2 generates all or part of the orthogonal matrix Q from an RQ factorization
determined by sgerqf

* @details
* \b Purpose:
   \verbatim
   ORGR2 generates an m by n real matrix Q with orthonormal rows,
   which is defined as the last m rows of a product of k elementary
   reflectors of order n

         Q  =  H(1) H(2) . . . H(k)

   as   returned by SGERQF.
   \endverbatim

* @param[in] M
         M is INTEGER \n
         The number of rows of the matrix Q. M >= 0. \n
* @param[in] N
         N is INTEGER \n
         The number of columns of the matrix Q. N >= M. \n
* @param[in] K
         K is INTEGER \n
         The number of elementary reflectors whose product defines the
         matrix Q. M >= K >= 0. \n
* @param[in,out] A
         A is REAL array, dimension (LDA,N) \n
         On entry, the (m-k+i)-th row must contain the vector which
         defines the elementary reflector H(i), for i = 1,2,...,k, as
         returned by SGERQF in the last k rows of its array argument
         A. \n
         On exit, the m by n matrix Q. \n
* @param[in] LDA
         LDA is INTEGER \n
         The first dimension of the array A. LDA >= fla_max(1,M). \n
* @param[in] TAU
         TAU is REAL array, dimension (K) \n
         TAU(i) must contain the scalar factor of the elementary
         reflector H(i), as   returned by SGERQF. \n
* @param[out] WORK
         WORK is REAL array, dimension (M) \n
* @param[out] INFO
         INFO is INTEGER \n
         = 0: successful exit \n
         < 0: if INFO = -i, the i-th argument has an illegal value \n

*  * */
    template <typename T>
    void orgr2(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *info)
    {
        orgr2(m, n, k, a, lda, tau, work, info);
    }
    template <typename T>
    void ungr2(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *info)
    {
        ungr2(m, n, k, a, lda, tau, work, info);
    }
    /**@} */ // end of orgr2

    /** @defgroup geql2 geql2
     * @ingroup LQ
     * @{
     */
    /*! @brief GEQL2 computes the QL factorization of a general rectangular matrix using an
 unblocked algorithm

 * @details
 * \b Purpose:
    \verbatim
    GEQL2 computes a QL factorization of a real m by n matrix A:
    A = Q * L.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the m by n matrix A. \n
          On exit, if m >= n, the lower triangle of the subarray
          A(m-n+1:m,1:n) contains the n by n lower triangular matrix L;
          if m <= n, the elements on and below the (n-m)-th
          superdiagonal contain the m by n lower trapezoidal matrix L;
          the remaining elements, with the array TAU, represent the
          orthogonal matrix Q as a product of elementary reflectors
          (see Further Details). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] TAU
          TAU is REAL array, dimension (min(M,N)) \n
          The scalar factors of the elementary reflectors (see Further
          Details). \n
 * @param[out] WORK
          WORK is REAL array, dimension (N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void geql2(integer *m, integer *n, T *a, integer *lda, T *tau, T *work, integer *info)
    {
        geql2(m, n, a, lda, tau, work, info);
    }
    /** @}*/ // end of

    /** @defgroup orgql {un,or}gql
     * @ingroup LQ
     * @{
     */
    /*! @brief ORGQL generates an M-by-N real matrix Q with orthonormal columns

 * @details
 * \b Purpose:
    \verbatim
    ORGQL generates an M-by-N real matrix Q with orthonormal columns,
    which is defined as the last N columns of a product of K elementary
    reflectors of order M
          Q  =  H(k) . . . H(2) H(1)
    as returned by GEQLF.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix Q. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix Q. M >= N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines the
          matrix Q. N >= K >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the (n-k+i)-th column must contain the vector which
          defines the elementary reflector H(i), for i = 1,2,...,k, as
          returned by SGEQLF in the last k columns of its array
          argument A. \n
          On exit, the M-by-N matrix Q. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The first dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[in] TAU
          TAU is REAL array, dimension (K) \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as returned by GEQLF. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= fla_max(1,N).
          For optimum performance LWORK >= N*NB, where NB is the
          optimal blocksize. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument has an illegal value \n

 *  * */
    template <typename T>
    void orgql(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        orgql(m, n, k, a, lda, tau, work, lwork, info);
    }
    template <typename T>
    void ungql(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        ungql(m, n, k, a, lda, tau, work, lwork, info);
    }
    /** @}*/ // end of orgql

    /** @defgroup ormql {un,or}mql
     * @ingroup LQ
     * @{
     */
    /*! @brief ORMQL overwrites the general real M-by-N matrix

 * @details
 * \b Purpose:
    \verbatim
     ORMQL overwrites the general real M-by-N matrix C with

                     SIDE = 'L'     SIDE = 'R'
     TRANS = 'N':      Q * C          C * Q
     TRANS = 'T':      Q**T * C       C * Q**T

     where Q is a real orthogonal matrix defined as the product of k
     elementary reflectors

           Q = H(k) . . . H(2) H(1)

     as returned by SGEQLF. Q is of order M if SIDE = 'L' and of order N
     if SIDE = 'R'.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q**T from the Left; \n
          = 'R': apply Q or Q**T from the Right. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N':  No transpose, apply Q; \n
          = 'T':  Transpose, apply Q**T. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines
          the matrix Q. \n
          If SIDE = 'L', M >= K >= 0; \n
          if SIDE = 'R', N >= K >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,K) \n
          The i-th column must contain the vector which defines the
          elementary reflector H(i), for i = 1,2,...,k, as returned by
          SGEQLF in the last k columns of its array argument A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. \n
          If SIDE = 'L', LDA >= fla_max(1,M); \n
          if SIDE = 'R', LDA >= fla_max(1,N). \n
 * @param[in] TAU
          TAU is REAL array, dimension (K) \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as returned by SGEQLF. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the M-by-N matrix C. \n
          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If SIDE = 'L', LWORK >= fla_max(1,N); \n
          if SIDE = 'R', LWORK >= fla_max(1,M). \n
          For good performance, LWORK should generally be larger. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void ormql(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        ormql(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
    template <typename T>
    void unmql(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        unmql(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
    /** @}*/ // end of ormql

    /** @defgroup org2l {un,or}g2l
     * @ingroup LQ
     * @{
     */
    /*! @brief ORG2L generates all or part of the orthogonal matrix Q from a QL factorization
 determined by sgeqlf.

 * @details
 * \b Purpose:
    \verbatim
    ORG2L generates an m by n real matrix Q with orthonormal columns,
    which is defined as the last n columns of a product of k elementary
    reflectors of order m

          Q  =  H(k) . . . H(2) H(1)

    as   returned by GEQLF.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix Q. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix Q. M >= N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines the
          matrix Q. N >= K >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the (n-k+i)-th column must contain the vector which
          defines the elementary reflector H(i), for i = 1,2,...,k, as
          returned by SGEQLF in the last k columns of its array
          argument A. \n
          On exit, the m by n matrix Q. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The first dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[in] TAU
          TAU is REAL array, dimension (K) \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as   returned by SGEQLF. \n
 * @param[out] WORK
          WORK is REAL array, dimension (N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument has an illegal value \n

 *  * */
    template <typename T>
    void org2l(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *info)
    {
        org2l(m, n, k, a, lda, tau, work, info);
    }
    template <typename T>
    void ung2l(integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *info)
    {
        ung2l(m, n, k, a, lda, tau, work, info);
    }
    /** @}*/ // end of org2l

    /** @defgroup orm2l orm2l
     * @ingroup LQ
     * @{
     */
    /*! @brief ORM2L multiplies a general matrix by the orthogonal matrix from a QL factorization
 determined by sgeqlf

 * @details
 * \b Purpose:
    \verbatim
    ORM2L overwrites the general real m by n matrix C with

          Q * C  if SIDE = 'L' and TRANS = 'N', or

          Q**T * C  if SIDE = 'L' and TRANS = 'T', or

          C * Q  if SIDE = 'R' and TRANS = 'N', or

          C * Q**T if SIDE = 'R' and TRANS = 'T',

    where Q is a real orthogonal matrix defined as the product of k
    elementary reflectors

          Q = H(k) . . . H(2) H(1)

    as   returned by SGEQLF. Q is of order m if SIDE = 'L' and of order n
    if SIDE = 'R'.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q**T from the Left \n
          = 'R': apply Q or Q**T from the Right \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N': apply Q  (No transpose) \n
          = 'T': apply Q**T (Transpose) \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines
          the matrix Q. \n
          If SIDE = 'L', M >= K >= 0; \n
          if SIDE = 'R', N >= K >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,K) \n
          The i-th column must contain the vector which defines the
          elementary reflector H(i), for i = 1,2,...,k, as   returned by
          SGEQLF in the last k columns of its array argument A.
          A is modified by the routine but restored on exit. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. \n
          If SIDE = 'L', LDA >= fla_max(1,M); \n
          if SIDE = 'R', LDA >= fla_max(1,N). \n
 * @param[in] TAU
          TAU is REAL array, dimension (K) \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as   returned by SGEQLF. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the m by n matrix C. \n
          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out] WORK
          WORK is REAL array, dimension \n
                                   (N) if SIDE = 'L', \n
                                   (M) if SIDE = 'R' \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void orm2l(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *info)
    {
        orm2l(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    }
    template <typename T>
    void unm2l(char *side, char *trans, integer *m, integer *n, integer *k, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *info)
    {
        unm2l(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    }
    /** @}*/ // end of orm2l

    /** @defgroup ggrqf ggrqf
     * @ingroup LQ
     * @{
     */
    /*! @brief GGRQF computes a generalized RQ factorization of an M-by-N matrix A and a P-by-N
 matrix B
 * @details
 * \b Purpose:
    \verbatim
     GGRQF computes a generalized RQ factorization of an M-by-N matrix A
     and a P-by-N matrix B:
                 A = R*Q,        B = Z*T*Q,

     where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal
     matrix, and R and T assume one of the forms:
     if M <= N,  R = ( 0  R12) M,   or if M > N,  R = ( R11) M-N,
                      N-M  M                           ( R21) N
                                                          N
     where R12 or R21 is upper triangular, and
     if P >= N,  T = ( T11) N  ,   or if P < N,  T = ( T11  T12) P,
                     (  0 ) P-N                         P   N-P
                        N
     where T11 is upper triangular.
     In particular, if B is square and nonsingular, the GRQ factorization
     of A and B implicitly gives the RQ factorization of A*inv(B):
                  A*inv(B) = (R*inv(T))*Z**T
     where inv(B) denotes the inverse of the matrix B, and Z**T denotes the
     transpose of the matrix Z.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows of the matrix B.  P >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrices A and B. N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit, if M <= N, the upper triangle of the subarray
          A(1:M,N-M+1:N) contains the M-by-M upper triangular matrix R;
          if M > N, the elements on and above the (M-N)-th subdiagonal
          contain the M-by-N upper trapezoidal matrix R; the remaining
          elements, with the array TAUA, represent the orthogonal
          matrix Q as a product of elementary reflectors (see Further
          Details). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[out] TAUA
          TAUA is REAL array, dimension (min(M,N)) \n
          The scalar factors of the elementary reflectors which
          represent the orthogonal matrix Q (see Further Details). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the P-by-N matrix B. \n
          On exit, the elements on and above the diagonal of the array
          contain the min(P,N)-by-N upper trapezoidal matrix T (T is
          upper triangular if P >= N); the elements below the diagonal,
          with the array TAUB, represent the orthogonal matrix Z as a
          product of elementary reflectors (see Further Details). \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,P). \n
 * @param[out] TAUB
          TAUB is REAL array, dimension (min(P,N)) \n
          The scalar factors of the elementary reflectors which
          represent the orthogonal matrix Z (see Further Details). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= fla_max(1,N,M,P). \n
          For optimum performance LWORK >= fla_max(N,M,P)*max(NB1,NB2,NB3),
          where NB1 is the optimal blocksize for the RQ factorization
          of an M-by-N matrix, NB2 is the optimal blocksize for the
          QR factorization of a P-by-N matrix, and NB3 is the optimal
          blocksize for a call of SORMRQ. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INF0= -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void ggrqf(integer *m, integer *p, integer *n, T *a, integer *lda, T *taua, T *b, integer *ldb,
               T *taub, T *work, integer *lwork, integer *info)
    {
        ggrqf(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info);
    }
    /** @}*/ // end of ggrqf

    /** @defgroup tzrzf tzrzf
     * @ingroup LQ
     * @{
     */
    /*! @brief TZRZF reduces the M-by-N ( M<=N) real upper trapezoidal matrix A \n
     to upper triangular form by means of orthogonal transformations.
 * @details
 * \b Purpose:
    \verbatim
    TZRZF reduces the M-by-N ( M<=N) real upper trapezoidal matrix A
    to upper triangular form by means of orthogonal transformations.

    The upper trapezoidal matrix A is factored as

       A = ( R  0) * Z,

    where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
    triangular matrix.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= M. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the leading M-by-N upper trapezoidal part of the
          array A must contain the matrix to be factorized.
          On exit, the leading M-by-M upper triangular part of A
          contains the upper triangular matrix R, and elements M+1 to
          N of the first M rows of A, with the array TAU, represent the
          orthogonal matrix Z as a product of M elementary reflectors. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] TAU
          TAU is REAL array, dimension (M) \n
          The scalar factors of the elementary reflectors. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,M).
          For optimum performance LWORK >= M*NB, where NB is
          the optimal blocksize. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tzrzf(integer *m, integer *n, T *a, integer *lda, T *tau, T *work, integer *lwork,
               integer *info)
    {
        tzrzf(m, n, a, lda, tau, work, lwork, info);
    }
    /** @}*/ // end of tzrzf

    /** @defgroup latrz latrz
     * @ingroup LQ
     * @{
     */
    /*! @brief LATRZ factors an upper trapezoidal matrix by means of orthogonal transformations

 * @details
 * \b Purpose:
    \verbatim
    LATRZ factors the M-by-(M+L) real upper trapezoidal matrix
    [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as (R  0) * Z, by means
    of orthogonal transformations.  Z is an (M+L)-by-(M+L) orthogonal
    matrix and, R and A1 are M-by-M upper triangular matrices.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] L
          L is INTEGER \n
          The number of columns of the matrix A containing the
          meaningful part of the Householder vectors. N-M >= L >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the leading M-by-N upper trapezoidal part of the
          array A must contain the matrix to be factorized.
          On exit, the leading M-by-M upper triangular part of A
          contains the upper triangular matrix R, and elements N-L+1 to
          N of the first M rows of A, with the array TAU, represent the
          orthogonal matrix Z as a product of M elementary reflectors. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] TAU
          TAU is REAL array, dimension (M) \n
          The scalar factors of the elementary reflectors. \n
 * @param[out] WORK
          WORK is REAL array, dimension (M)  \n

 *  * */
    template <typename T>
    void latrz(integer *m, integer *n, integer *l, T *a, integer *lda, T *tau, T *work)
    {
        latrz(m, n, l, a, lda, tau, work);
    }
    /** @}*/ // end of latrz

    /** @defgroup larz larz
     * @ingroup LQ
     * @{
     */
    /*! @brief LARZ applies an elementary reflector (as   returned by stzrzf) to a general matrix

 * @details
 * \b Purpose:
    \verbatim
    LARZ applies a real elementary reflector H to a real M-by-N
    matrix C, from either the left or the right. H is represented in the
    form

          H = I - tau * v * v**T

    where tau is a real scalar and v is a real vector.

    If tau = 0, then H is taken to be the unit matrix.


    H is a product of k elementary reflectors as   returned by STZRZF.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': form  H * C \n
          = 'R': form  C * H \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. \n
 * @param[in] L
          L is INTEGER \n
          The number of entries of the vector V containing
          the meaningful part of the Householder vectors. \n
          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. \n
 * @param[in] V
          V is REAL array, dimension (1+(L-1)*abs(INCV)) \n
          The vector v in the representation of H as   returned by
          STZRZF. V is not used if TAU = 0. \n
 * @param[in] INCV
          INCV is INTEGER \n
          The increment between elements of v. INCV <> 0. \n
 * @param[in] TAU
          TAU is REAL \n
          The value tau in the representation of H. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the M-by-N matrix C.
          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
          or C * H if SIDE = 'R'. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out] WORK
          WORK is REAL array, dimension \n
                         (N) if SIDE = 'L' \n
                      or (M) if SIDE = 'R'   \n

 *  * */
    template <typename T>
    void larz(char *side, integer *m, integer *n, integer *l, T *v, integer *incv, T *tau, T *c,
              integer *ldc, T *work)
    {
        larz(side, m, n, l, v, incv, tau, c, ldc, work);
    }
    /** @}*/ // end of larz

    /** @defgroup larzb larzb
     * @ingroup LQ
     * @{
     */
    /*! @brief LARZB applies a block reflector or its transpose to a general matrix

 * @details
 * \b Purpose:
    \verbatim
    LARZB applies a real block reflector H or its transpose H**T to
    a real distributed M-by-N  C from the left or the right.

    Currently, only STOREV = 'R' and DIRECT = 'B' are supported.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply H or H**T from the Left \n
          = 'R': apply H or H**T from the Right \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N': apply H (No transpose) \n
          = 'C': apply H**T (Transpose) \n
 * @param[in] DIRECT
          DIRECT is CHARACTER*1 \n
          Indicates how H is formed from a product of elementary
          reflectors \n
          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet) \n
          = 'B': H = H(k) . . . H(2) H(1) (Backward) \n
 * @param[in] STOREV
          STOREV is CHARACTER*1 \n
          Indicates how the vectors which define the elementary
          reflectors are stored: \n
          = 'C': Columnwise (not supported yet) \n
          = 'R': Rowwise \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. \n
 * @param[in] K
          K is INTEGER \n
          The order of the matrix T (= the number of elementary
          reflectors whose product defines the block reflector). \n
 * @param[in] L
          L is INTEGER \n
          The number of columns of the matrix V containing the
          meaningful part of the Householder reflectors. \n
          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. \n
 * @param[in] V
          V is REAL array, dimension (LDV,NV). \n
          If STOREV = 'C', NV = K; if STOREV = 'R', NV = L. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V. \n
          If STOREV = 'C', LDV >= L; if STOREV = 'R', LDV >= K. \n
 * @param[in] T
          T is REAL array, dimension (LDT,K) \n
          The triangular K-by-K matrix T in the representation of the
          block reflector. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= K. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the M-by-N matrix C.
          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out] WORK
          WORK is REAL array, dimension (LDWORK,K) \n
 * @param[in] LDWORK
          LDWORK is INTEGER \n
          The leading dimension of the array WORK. \n
          If SIDE = 'L', LDWORK >= fla_max(1,N); \n
          if SIDE = 'R', LDWORK >= fla_max(1,M). \n

 *  * */
    template <typename T>
    void larzb(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
               integer *k, integer *l, T *v, integer *ldv, T *t, integer *ldt, T *c, integer *ldc,
               T *work, integer *ldwork)
    {
        larzb(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work, ldwork);
    }
    /** @}*/ // end of larzb

    /** @defgroup larzt larzt
     * @ingroup LQ
     * @{
     */
    /*! @brief LARZT forms the triangular factor T of a block reflector H = I - vtvH.

 * @details
 * \b Purpose:
    \verbatim
    LARZT forms the triangular factor T of a real block reflector
    H of order > n, which is defined as a product of k elementary
    reflectors.

    If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;

    If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.

    If STOREV = 'C', the vector which defines the elementary reflector
    H(i) is stored in the i-th column of the array V, and

       H  =  I - V * T * V**T

    If STOREV = 'R', the vector which defines the elementary reflector
    H(i) is stored in the i-th row of the array V, and

       H  =  I - V**T * T * V

    Currently, only STOREV = 'R' and DIRECT = 'B' are supported.
    \endverbatim

 * @param[in] DIRECT
          DIRECT is CHARACTER*1 \n
          Specifies the order in which the elementary reflectors are
          multiplied to form the block reflector: \n
          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet) \n
          = 'B': H = H(k) . . . H(2) H(1) (Backward) \n
 * @param[in] STOREV
          STOREV is CHARACTER*1 \n
          Specifies how the vectors which define the elementary
          reflectors are stored (see also Further Details): \n
          = 'C': columnwise (not supported yet) \n
          = 'R': rowwise \n
 * @param[in] N
          N is INTEGER \n
          The order of the block reflector H. N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The order of the triangular factor T (= the number of
          elementary reflectors). K >= 1. \n
 * @param[in,out] V
          V is REAL array, dimension \n
                               (LDV,K) if STOREV = 'C' \n
                               (LDV,N) if STOREV = 'R' \n
          The matrix V. See further details. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V.
          If STOREV = 'C', LDV >= fla_max(1,N); if STOREV = 'R', LDV >= K. \n
 * @param[in] TAU
          TAU is REAL array, dimension (K) \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i). \n
 * @param[out] T
          T is REAL array, dimension (LDT,K) \n
          The k by k triangular factor T of the block reflector.
          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
          lower triangular. The rest of the array is not used. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= K. \n

 *  * */
    template <typename T>
    void larzt(char *direct, char *storev, integer *n, integer *k, T *v, integer *ldv, T *tau, T *t,
               integer *ldt)
    {
        larzt(direct, storev, n, k, v, ldv, tau, t, ldt);
    }
    /** @}*/ // end of larzt

    /** @defgroup ormzr {un,or}mzr
     * @ingroup LQ
     * @{
     */
    /*! @brief ORMRZ overwrites the general real M-by-N matrix

 * @details
 * \b Purpose:
    \verbatim
     ORMRZ overwrites the general real M-by-N matrix C with

                     SIDE = 'L'     SIDE = 'R'
     TRANS = 'N':      Q * C          C * Q
     TRANS = 'T':      Q**T * C       C * Q**T

     where Q is a real orthogonal matrix defined as the product of k
     elementary reflectors

           Q = H(1) H(2) . . . H(k)

     as returned by STZRZF. Q is of order M if SIDE = 'L' and of order N
     if SIDE = 'R'.
    \endverbatim
*
*  Arguments:
*  ==========
*
 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q**T from the Left; \n
          = 'R': apply Q or Q**T from the Right. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N':  No transpose, apply Q; \n
          = 'T':  Transpose, apply Q**T. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines
          the matrix Q. \n
          If SIDE = 'L', M >= K >= 0; \n
          if SIDE = 'R', N >= K >= 0. \n
 * @param[in] L
          L is INTEGER \n
          The number of columns of the matrix A containing
          the meaningful part of the Householder reflectors. \n
          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. \n
 * @param[in] A
          A is REAL array, dimension \n
                               (LDA,M) if SIDE = 'L', \n
                               (LDA,N) if SIDE = 'R' \n
          The i-th row must contain the vector which defines the
          elementary reflector H(i), for i = 1,2,...,k, as returned by
          STZRZF in the last k rows of its array argument A.
          A is modified by the routine but restored on exit. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,K). \n
 * @param[in] TAU
          TAU is REAL array, dimension (K) \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as returned by STZRZF. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the M-by-N matrix C. \n
          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If SIDE = 'L', LWORK >= fla_max(1,N); \n
          if SIDE = 'R', LWORK >= fla_max(1,M). \n
          For good performance, LWORK should generally be larger. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void ormrz(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, T *a,
               integer *lda, T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        ormrz(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info);
    }
    template <typename T>
    void unmrz(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, T *a,
               integer *lda, T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        unmrz(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info);
    }
    /** @}*/ // end of ormzr

    /** @defgroup ormr3 ormr3
     * @ingroup LQ
     * @{
     */
    /*! @brief ORMR3 multiplies a general matrix by the orthogonal matrix from a RZ factorization
 determined by stzrzf

 * @details
 * \b Purpose:
    \verbatim
    ORMR3 overwrites the general real m by n matrix C with

          Q * C  if SIDE = 'L' and TRANS = 'N', or

          Q**T* C  if SIDE = 'L' and TRANS = 'C', or

          C * Q  if SIDE = 'R' and TRANS = 'N', or

          C * Q**T if SIDE = 'R' and TRANS = 'C',

    where Q is a real orthogonal matrix defined as the product of k
    elementary reflectors

          Q = H(1) H(2) . . . H(k)

    as   returned by STZRZF. Q is of order m if SIDE = 'L' and of order n
    if SIDE = 'R'.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q**T from the Left \n
          = 'R': apply Q or Q**T from the Right \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N': apply Q  (No transpose) \n
          = 'T': apply Q**T (Transpose) \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. N >= 0. \n
 * @param[in] K
          K is INTEGER \n
          The number of elementary reflectors whose product defines
          the matrix Q. \n
          If SIDE = 'L', M >= K >= 0; \n
          if SIDE = 'R', N >= K >= 0. \n
 * @param[in] L
          L is INTEGER \n
          The number of columns of the matrix A containing
          the meaningful part of the Householder reflectors. \n
          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. \n
 * @param[in] A
          A is REAL array, dimension \n
                               (LDA,M) if SIDE = 'L', \n
                               (LDA,N) if SIDE = 'R' \n
          The i-th row must contain the vector which defines the
          elementary reflector H(i), for i = 1,2,...,k, as   returned by
          STZRZF in the last k rows of its array argument A. \n
          A is modified by the routine but restored on exit. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,K). \n
 * @param[in] TAU
          TAU is REAL array, dimension (K) \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as   returned by STZRZF. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the m-by-n matrix C. \n
          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out] WORK
          WORK is REAL array, dimension \n
                                   (N) if SIDE = 'L', \n
                                   (M) if SIDE = 'R' \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void ormr3(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, T *a,
               integer *lda, T *tau, T *c, integer *ldc, T *work, integer *info)
    {
        ormr3(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info);
    }
    template <typename T>
    void unmr3(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, T *a,
               integer *lda, T *tau, T *c, integer *ldc, T *work, integer *info)
    {
        unmr3(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info);
    }
    /** @}*/ // end of ormr3
    /** @} */ // end of LQ

    /** @defgroup CS Cosine-Sine (CS) decomposition
     * @ingroup Orthogonal
     * @{
     */

    /** @defgroup bbscd bbscd
     * @ingroup CS
     * @{
     */
    /*! @brief Bidiagonal block cs decomposition of orthogonal/unitary matrix

    @details
    \b Purpose:
    \verbatim
    bbcsd() is inline function to call respective LAPACK API using templates.
    Where T can be REAL or DOUBLE PRECISION.

    X is M-by-M, its top-left block is P-by-Q, and Q must be no larger
    than P, M-P, or M-Q.
    BBCSD computes the CS decomposition of an orthogonal matrix in
    bidiagonal-block form,

       [ B11 | B12 0  0 ]
       [  0  |  0 -I  0 ]
    X = [----------------]
       [ B21 | B22 0  0 ]
       [  0  |  0  0  I ]

                                 [  C | -S  0  0 ]
                     [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**T
                   = [---------] [---------------] [---------]   .
                     [    | U2 ] [  S |  C  0  0 ] [    | V2 ]
                                 [  0 |  0  0  I ]

    X is M-by-M, its top-left block is P-by-Q, and Q must be no larger
    than P, M-P, or M-Q. (If Q is not the smallest index, then X must be
    transposed and/or permuted. This can be done in constant time using
    the TRANS and SIGNS options. See SORCSD for details.)

    The bidiagonal matrices B11, B12, B21, and B22 are represented
    implicitly by angles THETA(1:Q) and PHI(1:Q-1).

    The orthogonal matrices U1, U2, V1T, and V2T are input/output.
    The input matrices are pre- or post-multiplied by the appropriate
    singular vector matrices.

    Real reference:
    http://www.netlib.org/lapack/explore-html/d1/df5/group__real_o_t_h_e_rcomputational_ga95bdd6e44aed23173e9a0c93c32dad78.html#ga95bdd6e44aed23173e9a0c93c32dad78
    Double reference:
    http://www.netlib.org/lapack/explore-html/da/dba/group__double_o_t_h_e_rcomputational_ga27a367582a76c7b48a8bf3eed068e216.html#ga27a367582a76c7b48a8bf3eed068e216
    Complex reference:
      http://www.netlib.org/lapack/explore-html/d3/db9/group__complex_o_t_h_e_rcomputational_gaa78ee3c0b2912f780622143726a5299e.html#gaa78ee3c0b2912f780622143726a5299e
      Complex double reference:
      http://www.netlib.org/lapack/explore-html/d0/da6/group__complex16_o_t_h_e_rcomputational_gab100b320bf854584daf3579ff6d96485.html#gab100b320bf854584daf3579ff6d96485
      \endverbatim

 * @param[in]	JOBU1
          JOBU1 is CHARACTER \n
          = 'Y':      U1 is updated; \n
          otherwise:  U1 is not updated. \n
 * @param[in]	JOBU2
          JOBU2 is CHARACTER \n
          = 'Y':      U2 is updated; \n
          otherwise:  U2 is not updated. \n
 * @param[in]	JOBV1T
          JOBV1T is CHARACTER \n
          = 'Y':      V1T is updated; \n
          otherwise:  V1T is not updated. \n
 * @param[in]	JOBV2T
          JOBV2T is CHARACTER \n
          = 'Y':      V2T is updated; \n
          otherwise:  V2T is not updated. \n
 * @param[in]	TRANS
          TRANS is CHARACTER \n
          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major
                      order; \n
          otherwise:  X, U1, U2, V1T, and V2T are stored in column-
                      major order. \n
 * @param[in]	M
          M is INTEGER \n
          The number of rows and columns in X, the orthogonal matrix in
          bidiagonal-block form. \n
 * @param[in]	P
          P is INTEGER \n
          The number of rows in the top-left block of X. 0 <= P <= M. \n
 * @param[in]	Q
          Q is INTEGER \n
          The number of columns in the top-left block of X.
          0 <= Q <= MIN(P,M-P,M-Q). \n
 * @param[in,out]	THETA
          THETA is REAL or DOUBLE PRECISION array, dimension (Q) \n
          On entry, the angles THETA(1),...,THETA(Q) that, along with
          PHI(1), ...,PHI(Q-1), define the matrix in bidiagonal-block
          form. On exit, the angles whose cosines and sines define the
          diagonal blocks in the CS decomposition. \n
 * @param[in,out]	PHI
          PHI is REAL or DOUBLE PRECISION array, dimension (Q-1) \n
          The angles PHI(1),...,PHI(Q-1) that, along with THETA(1),...,
          THETA(Q), define the matrix in bidiagonal-block form. \n
 * @param[in,out]	U1
          U1 is REAL or DOUBLE PRECISION or COMPLEX or COMPLEX*16 array, dimension (LDU1,P) \n
          On entry, a P-by-P matrix. On exit, U1 is postmultiplied
          by the left singular vector matrix common to [ B11 ; 0 ] and
          [ B12 0 0 ; 0 -I 0 0 ]. \n
 * @param[in]	LDU1
          LDU1 is INTEGER \n
          The leading dimension of the array U1, LDU1 >= MAX(1,P). \n
 * @param[in,out]	U2
          U2 is REAL or DOUBLE PRECISION or COMPLEX or COMPLEX*16 array, dimension (LDU2,M-P) \n
          On entry, an (M-P)-by-(M-P) matrix. On exit, U2 is
          postmultiplied by the left singular vector matrix common to
          [ B21 ; 0 ] and [ B22 0 0 ; 0 0 I ]. \n
 * @param[in]	LDU2
          LDU2 is INTEGER \n
          The leading dimension of the array U2, LDU2 >= MAX(1,M-P). \n
 * @param[in,out]	V1T
          V1T is REAL or DOUBLE PRECISION or COMPLEX or COMPLEX*16 array, dimension (LDV1T,Q) \n
          On entry, a Q-by-Q matrix. On exit, V1T is premultiplied
          by the transpose of the right singular vector
          matrix common to [ B11 ; 0 ] and [ B21 ; 0 ]. \n
 * @param[in]	LDV1T
          LDV1T is INTEGER \n
          The leading dimension of the array V1T, LDV1T >= MAX(1,Q). \n
 * @param[in,out]	V2T
          V2T is REAL or DOUBLE PRECISION or COMPLEX or COMPLEX*16 array, dimension (LDV2T,M-Q) \n
          On entry, an (M-Q)-by-(M-Q) matrix. On exit, V2T is
          premultiplied by the transpose of the right
          singular vector matrix common to [ B12 0 0 ; 0 -I 0 ] and
          [ B22 0 0 ; 0 0 I ]. \n
 * @param[in]	LDV2T
          LDV2T is INTEGER \n
          The leading dimension of the array V2T, LDV2T >= MAX(1,M-Q). \n
 * @param[out]	B11D
          B11D is REAL or DOUBLE PRECISION array, dimension (Q) \n
          When SBBCSD converges, B11D contains the cosines of THETA(1),
          ..., THETA(Q). If SBBCSD fails to converge, then B11D
          contains the diagonal of the partially reduced top-left
          block. \n
 * @param[out]	B11E
          B11E is REAL or DOUBLE PRECISION array, dimension (Q-1) \n
          When SBBCSD converges, B11E contains zeros. If SBBCSD fails
          to converge, then B11E contains the superdiagonal of the
          partially reduced top-left block. \n
 * @param[out]	B12D
          B12D is REAL or DOUBLE PRECISION array, dimension (Q) \n
          When SBBCSD converges, B12D contains the negative sines of
          THETA(1), ..., THETA(Q). If SBBCSD fails to converge, then
          B12D contains the diagonal of the partially reduced top-right
          block. \n
 * @param[out]	B12E
          B12E is REAL or DOUBLE PRECISION array, dimension (Q-1) \n
          When SBBCSD converges, B12E contains zeros. If SBBCSD fails
          to converge, then B12E contains the subdiagonal of the
          partially reduced top-right block. \n
 * @param[out]	B21D
          B21D is REAL or DOUBLE PRECISION array, dimension (Q) \n
          When SBBCSD converges, B21D contains the negative sines of
          THETA(1), ..., THETA(Q). If SBBCSD fails to converge, then
          B21D contains the diagonal of the partially reduced bottom-left
          block. \n
 * @param[out]	B21E
          B21E is REAL or DOUBLE PRECISION array, dimension (Q-1) \n
          When SBBCSD converges, B21E contains zeros. If SBBCSD fails
          to converge, then B21E contains the subdiagonal of the
          partially reduced bottom-left block. \n
 * @param[out]	B22D
          B22D is REAL or DOUBLE PRECISION array, dimension (Q) \n
          When SBBCSD converges, B22D contains the negative sines of
          THETA(1), ..., THETA(Q). If SBBCSD fails to converge, then
          B22D contains the diagonal of the partially reduced bottom-right
          block. \n
 * @param[out]	B22E
          B22E is REAL or DOUBLE PRECISION array, dimension (Q-1) \n
          When SBBCSD converges, B22E contains zeros. If SBBCSD fails
          to converge, then B22E contains the subdiagonal of the
          partially reduced bottom-right block. \n
 * @param[out]	WORK
          WORK is REAL or DOUBLE PRECISION array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= MAX(1,8*Q). \n
 \n
          If LWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of the WORK array,
          returns this value as the first entry of the work array, and
          no error message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if SBBCSD did not converge, INFO specifies the number
                of nonzero entries in PHI, and B11D, B11E, etc.,
                contain the partially reduced matrix. \n

 *  * */
    template <typename T>
    void bbcsd(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, integer *m,
               integer *p, integer *q, T *theta, T *phi, T *u1, integer *ldu1, T *u2, integer *ldu2,
               T *v1t, integer *ldv1t, T *v2t, integer *ldv2t, T *b11d, T *b11e, T *b12d, T *b12e,
               T *b21d, T *b21e, T *b22d, T *b22e, T *work, integer *lwork, integer *info)
    {
        bbcsd(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t,
              ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, work, lwork, info);
    }

    template <typename T, typename Ta>
    void bbcsd(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, integer *m,
               integer *p, integer *q, Ta *theta, Ta *phi, T *u1, integer *ldu1, T *u2,
               integer *ldu2, T *v1t, integer *ldv1t, T *v2t, integer *ldv2t, Ta *b11d, Ta *b11e,
               Ta *b12d, Ta *b12e, Ta *b21d, Ta *b21e, Ta *b22d, Ta *b22e, Ta *rwork,
               integer *lrwork, integer *info)
    {
        bbcsd(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t,
              ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, rwork, lrwork,
              info);
    }
    /** @}*/ // end of bbcsd

    /** @defgroup orcsd {un,or}csd
     * @ingroup CS
     * @{
     */
    /*! @brief ORCSD computes the CS decomposition of an M-by-M partitioned orthogonal/unitary
 matrix X

 * @details
 * \b Purpose:
    \verbatim
     ORCSD computes the CS decomposition of an M-by-M partitioned
     orthogonal/unitary  matrix X:

                                     [  I  0  0 |  0  0  0 ]
                                     [  0  C  0 |  0 -S  0 ]
         [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**T
     X = [-----------] = [---------] [---------------------] [---------]   .
         [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
                                     [  0  S  0 |  0  C  0 ]
                                     [  0  0  I |  0  0  0 ]

     X11 is P-by-Q. The orthogonal/unitary  matrices U1, U2, V1, and V2 are P-by-P,
     (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are
     R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
     which R = MIN(P,M-P,Q,M-Q).
    \endverbatim

 * @param[in] JOBU1
          JOBU1 is CHARACTER \n
          = 'Y':      U1 is computed; \n
          otherwise:  U1 is not computed. \n
 * @param[in] JOBU2
          JOBU2 is CHARACTER \n
          = 'Y':      U2 is computed; \n
          otherwise:  U2 is not computed. \n
 * @param[in] JOBV1T
          JOBV1T is CHARACTER \n
          = 'Y':      V1T is computed; \n
          otherwise:  V1T is not computed. \n
 * @param[in] JOBV2T
          JOBV2T is CHARACTER \n
          = 'Y':      V2T is computed; \n
          otherwise:  V2T is not computed. \n
 * @param[in] TRANS
          TRANS is CHARACTER \n
          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major
                      order; \n
          otherwise:  X, U1, U2, V1T, and V2T are stored in column-
                      major order. \n
 * @param[in] SIGNS
          SIGNS is CHARACTER \n
          = 'O':      The lower-left block is made nonpositive (the
                      "other" convention); \n
          otherwise:  The upper-right block is made nonpositive (the
                      "default" convention). \n
 * @param[in] M
          M is INTEGER \n
          The number of rows and columns in X. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows in X11 and X12. 0 <= P <= M. \n
 * @param[in] Q
          Q is INTEGER \n
          The number of columns in X11 and X21. 0 <= Q <= M. \n
 * @param[in,out] X11
          X11 is REAL array, dimension (LDX11,Q) \n
          On entry, part of the orthogonal matrix whose CSD is desired. \n
 * @param[in] LDX11
          LDX11 is INTEGER \n
          The leading dimension of X11. LDX11 >= MAX(1,P). \n
 * @param[in,out] X12
          X12 is REAL array, dimension (LDX12,M-Q) \n
          On entry, part of the orthogonal matrix whose CSD is desired. \n
 * @param[in] LDX12
          LDX12 is INTEGER \n
          The leading dimension of X12. LDX12 >= MAX(1,P). \n
 * @param[in,out] X21
          X21 is REAL array, dimension (LDX21,Q) \n
          On entry, part of the orthogonal matrix whose CSD is desired. \n
 * @param[in] LDX21
          LDX21 is INTEGER \n
          The leading dimension of X11. LDX21 >= MAX(1,M-P). \n
 * @param[in,out] X22
          X22 is REAL array, dimension (LDX22,M-Q) \n
          On entry, part of the orthogonal matrix whose CSD is desired. \n
 * @param[in] LDX22
          LDX22 is INTEGER \n
          The leading dimension of X11. LDX22 >= MAX(1,M-P). \n
 * @param[out] THETA
          THETA is REAL array, dimension (R), in which R =
          MIN(P,M-P,Q,M-Q). \n
          C = DIAG( COS(THETA(1)), ... , COS(THETA(R))) and \n
          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R))). \n
 * @param[out] U1
          U1 is REAL array, dimension (LDU1,P) \n
          If JOBU1 = 'Y', U1 contains the P-by-P orthogonal matrix U1. \n
 * @param[in] LDU1
          LDU1 is INTEGER \n
          The leading dimension of U1. If JOBU1 = 'Y', LDU1 >=
          MAX(1,P). \n
 * @param[out] U2
          U2 is REAL array, dimension (LDU2,M-P) \n
          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) orthogonal
          matrix U2. \n
 * @param[in] LDU2
          LDU2 is INTEGER \n
          The leading dimension of U2. If JOBU2 = 'Y', LDU2 >=
          MAX(1,M-P). \n
 * @param[out] V1T
          V1T is REAL array, dimension (LDV1T,Q) \n
          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix orthogonal
          matrix V1**T. \n
 * @param[in] LDV1T
          LDV1T is INTEGER \n
          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >=
          MAX(1,Q). \n
 * @param[out] V2T
          V2T is REAL array, dimension (LDV2T,M-Q) \n
          If JOBV2T = 'Y', V2T contains the (M-Q)-by-(M-Q) orthogonal
          matrix V2**T. \n
 * @param[in] LDV2T
          LDV2T is INTEGER \n
          The leading dimension of V2T. If JOBV2T = 'Y', LDV2T >=
          MAX(1,M-Q). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
          If INFO > 0 on exit, WORK(2:R) contains the values PHI(1),
          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R),
          define the matrix in intermediate bidiagonal-block form
          remaining after nonconvergence. INFO specifies the number
          of nonzero PHI's. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the work array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (M-MIN(P, M-P, Q, M-Q)) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  SBBCSD did not converge. See the description of WORK
                above for details. \n

 *  * */
    template <typename T>
    void orcsd(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
               integer *m, integer *p, integer *q, T *x11, integer *ldx11, T *x12, integer *ldx12,
               T *x21, integer *ldx21, T *x22, integer *ldx22, T *theta, T *u1, integer *ldu1,
               T *u2, integer *ldu2, T *v1t, integer *ldv1t, T *v2t, integer *ldv2t, T *work,
               integer *lwork, integer *iwork, integer *info)
    {
        orcsd(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21,
              ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork,
              iwork, info);
    }
    template <typename T, typename Ta>
    void uncsd(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
               integer *m, integer *p, integer *q, T *x11, integer *ldx11, T *x12, integer *ldx12,
               T *x21, integer *ldx21, T *x22, integer *ldx22, Ta *theta, T *u1, integer *ldu1,
               T *u2, integer *ldu2, T *v1t, integer *ldv1t, T *v2t, integer *ldv2t, T *work,
               integer *lwork, Ta *rwork, integer *lrwork, integer *iwork, integer *info)
    {
        uncsd(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21,
              ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork,
              rwork, lrwork, iwork, info);
    }
    /** @}*/ // end of orcsd

    /** @defgroup orcsd2by1 {un,or}csd2by1
     * @ingroup CS
     * @{
     */
    /*! @brief ORCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with orthonormal
 columns

 * @details
 * \b Purpose:
    \verbatim
     ORCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with
     orthonormal columns that has been partitioned into a 2-by-1 block
     structure:

                                    [  I1 0  0 ]
                                    [  0  C  0 ]
              [ X11 ]   [ U1 |    ] [  0  0  0 ]
          X = [-----] = [---------] [----------] V1**T .
              [ X21 ]   [    | U2 ] [  0  0  0 ]
                                    [  0  S  0 ]
                                    [  0  0  I2]

     X11 is P-by-Q. The orthogonal matrices U1, U2, and V1 are P-by-P,
     (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R
     nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which
     R = MIN(P,M-P,Q,M-Q). I1 is a K1-by-K1 identity matrix and I2 is a
     K2-by-K2 identity matrix, where K1 = MAX(Q+P-M,0), K2 = MAX(Q-P,0).
    \endverbatim

 * @param[in] JOBU1
          JOBU1 is CHARACTER \n
          = 'Y':      U1 is computed; \n
          otherwise:  U1 is not computed. \n
 * @param[in] JOBU2
          JOBU2 is CHARACTER \n
          = 'Y':      U2 is computed; \n
          otherwise:  U2 is not computed. \n
 * @param[in] JOBV1T
          JOBV1T is CHARACTER \n
          = 'Y':      V1T is computed; \n
          otherwise:  V1T is not computed. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows in X. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows in X11. 0 <= P <= M. \n
 * @param[in] Q
          Q is INTEGER \n
          The number of columns in X11 and X21. 0 <= Q <= M. \n
 * @param[in,out] X11
          X11 is REAL array, dimension (LDX11,Q) \n
          On entry, part of the orthogonal matrix whose CSD is desired. \n
 * @param[in] LDX11
          LDX11 is INTEGER \n
          The leading dimension of X11. LDX11 >= MAX(1,P). \n
 * @param[in,out] X21
          X21 is REAL array, dimension (LDX21,Q) \n
          On entry, part of the orthogonal matrix whose CSD is desired. \n
 * @param[in] LDX21
          LDX21 is INTEGER \n
          The leading dimension of X21. LDX21 >= MAX(1,M-P). \n
 * @param[out] THETA
          THETA is REAL array, dimension (R), in which R =
          MIN(P,M-P,Q,M-Q). \n
          C = DIAG( COS(THETA(1)), ... , COS(THETA(R))) and \n
          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R))). \n
 * @param[out] U1
          U1 is REAL array, dimension (P) \n
          If JOBU1 = 'Y', U1 contains the P-by-P orthogonal matrix U1. \n
 * @param[in] LDU1
          LDU1 is INTEGER \n
          The leading dimension of U1. If JOBU1 = 'Y', LDU1 >=
          MAX(1,P). \n
 * @param[out] U2
          U2 is REAL array, dimension (M-P) \n
          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) orthogonal
          matrix U2. \n
 * @param[in] LDU2
          LDU2 is INTEGER \n
          The leading dimension of U2. If JOBU2 = 'Y', LDU2 >=
          MAX(1,M-P). \n
 * @param[out] V1T
          V1T is REAL array, dimension (Q) \n
          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix orthogonal
          matrix V1**T. \n
 * @param[in] LDV1T
          LDV1T is INTEGER \n
          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >=
          MAX(1,Q). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
          If INFO > 0 on exit, WORK(2:R) contains the values PHI(1),
          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R),
          define the matrix in intermediate bidiagonal-block form
          remaining after nonconvergence. INFO specifies the number
          of nonzero PHI's. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the work array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (M-MIN(P,M-P,Q,M-Q)) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  SBBCSD did not converge. See the description of WORK
                above for details. \n

 *  * */
    template <typename T>
    void orcsd2by1(char *jobu1, char *jobu2, char *jobv1t, integer *m, integer *p, integer *q,
                   T *x11, integer *ldx11, T *x21, integer *ldx21, T *theta, T *u1, integer *ldu1,
                   T *u2, integer *ldu2, T *v1t, integer *ldv1t, T *work, integer *lwork,
                   integer *iwork, integer *info)
    {
        orcsd2by1(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2,
                  v1t, ldv1t, work, lwork, iwork, info);
    }
    template <typename T, typename Ta>
    void uncsd2by1(char *jobu1, char *jobu2, char *jobv1t, integer *m, integer *p, integer *q,
                   T *x11, integer *ldx11, T *x21, integer *ldx21, Ta *theta, T *u1, integer *ldu1,
                   T *u2, integer *ldu2, T *v1t, integer *ldv1t, T *work, integer *lwork, Ta *rwork,
                   integer *lrwork, integer *iwork, integer *info)
    {
        uncsd2by1(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2,
                  v1t, ldv1t, work, lwork, rwork, lrwork, iwork, info);
    }
    /** @}*/ // end of

    /** @defgroup orbdb {un,or}bdb
     * @ingroup CS
     * @{
     */
    /*! @brief ORBDB simultaneously bidiagonalizes the blocks of an M-by-M \n
 partitioned orthogonal matrix X
 * @details
 * \b Purpose:
    \verbatim
     ORBDB simultaneously bidiagonalizes the blocks of an M-by-M
     partitioned orthogonal matrix X:

                                     [ B11 | B12 0  0 ]
         [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**T
     X = [-----------] = [---------] [----------------] [---------]   .
         [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]
                                     [  0  |  0  0  I ]

     X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is
     not the case, then X must be transposed and/or permuted. This can be
     done in constant time using the TRANS and SIGNS options. See SORCSD
     for details.)

     The orthogonal matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-
     (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are
     represented implicitly by Householder vectors.

     B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented
     implicitly by angles THETA, PHI.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER \n
          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major
                      order; \n
          otherwise:  X, U1, U2, V1T, and V2T are stored in column-
                      major order. \n
 * @param[in] SIGNS
          SIGNS is CHARACTER \n
          = 'O':      The lower-left block is made nonpositive (the
                      "other" convention); \n
          otherwise:  The upper-right block is made nonpositive (the
                      "default" convention). \n
 * @param[in] M
          M is INTEGER \n
          The number of rows and columns in X. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows in X11 and X12. 0 <= P <= M. \n
 * @param[in] Q
          Q is INTEGER \n
          The number of columns in X11 and X21. 0 <= Q <=
          MIN(P,M-P,M-Q). \n
 * @param[in,out] X11
          X11 is REAL array, dimension (LDX11,Q) \n
          On entry, the top-left block of the orthogonal matrix to be
          reduced. On exit, the form depends on TRANS: \n
          If TRANS = 'N', then \n
             the columns of tril(X11) specify reflectors for P1,
             the rows of triu(X11,1) specify reflectors for Q1; \n
          else TRANS = 'T', and \n
             the rows of triu(X11) specify reflectors for P1,
             the columns of tril(X11,-1) specify reflectors for Q1. \n
 * @param[in] LDX11
          LDX11 is INTEGER \n
          The leading dimension of X11. If TRANS = 'N', then LDX11 >=
          P; else LDX11 >= Q. \n
 * @param[in,out] X12
          X12 is REAL array, dimension (LDX12,M-Q) \n
          On entry, the top-right block of the orthogonal matrix to
          be reduced. On exit, the form depends on TRANS: \n
          If TRANS = 'N', then \n
             the rows of triu(X12) specify the first P reflectors for
             Q2; \n
          else TRANS = 'T', and \n
             the columns of tril(X12) specify the first P reflectors
             for Q2. \n
 * @param[in] LDX12
          LDX12 is INTEGER \n
          The leading dimension of X12. If TRANS = 'N', then LDX12 >=
          P; else LDX11 >= M-Q. \n
 * @param[in,out] X21
          X21 is REAL array, dimension (LDX21,Q) \n
          On entry, the bottom-left block of the orthogonal matrix to
          be reduced. On exit, the form depends on TRANS: \n
          If TRANS = 'N', then \n
             the columns of tril(X21) specify reflectors for P2; \n
          else TRANS = 'T', and \n
             the rows of triu(X21) specify reflectors for P2. \n
 * @param[in] LDX21
          LDX21 is INTEGER \n
          The leading dimension of X21. If TRANS = 'N', then LDX21 >=
          M-P; else LDX21 >= Q. \n
 * @param[in,out] X22
          X22 is REAL array, dimension (LDX22,M-Q) \n
          On entry, the bottom-right block of the orthogonal matrix to
          be reduced. On exit, the form depends on TRANS: \n
          If TRANS = 'N', then \n
             the rows of triu(X22(Q+1:M-P,P+1:M-Q)) specify the last
             M-P-Q reflectors for Q2, \n
          else TRANS = 'T', and \n
             the columns of tril(X22(P+1:M-Q,Q+1:M-P)) specify the last
             M-P-Q reflectors for P2. \n
 * @param[in] LDX22
          LDX22 is INTEGER \n
          The leading dimension of X22. If TRANS = 'N', then LDX22 >=
          M-P; else LDX22 >= M-Q. \n
 * @param[out] THETA
          THETA is REAL array, dimension (Q) \n
          The entries of the bidiagonal blocks B11, B12, B21, B22 can
          be computed from the angles THETA and PHI. See Further
          Details. \n
 * @param[out] PHI
          PHI is REAL array, dimension (Q-1) \n
          The entries of the bidiagonal blocks B11, B12, B21, B22 can
          be computed from the angles THETA and PHI. See Further
          Details. \n
 * @param[out] TAUP1
          TAUP1 is REAL array, dimension (P) \n
          The scalar factors of the elementary reflectors that define
          P1. \n
 * @param[out] TAUP2
          TAUP2 is REAL array, dimension (M-P) \n
          The scalar factors of the elementary reflectors that define
          P2. \n
 * @param[out] TAUQ1
          TAUQ1 is REAL array, dimension (Q) \n
          The scalar factors of the elementary reflectors that define
          Q1. \n
 * @param[out] TAUQ2
          TAUQ2 is REAL array, dimension (M-Q) \n
          The scalar factors of the elementary reflectors that define
          Q2. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (LWORK) \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= M-Q. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n


 *  * */
    template <typename T>
    void orbdb(char *trans, char *signs, integer *m, integer *p, integer *q, T *x11, integer *ldx11,
               T *x12, integer *ldx12, T *x21, integer *ldx21, T *x22, integer *ldx22, T *theta,
               T *phi, T *taup1, T *taup2, T *tauq1, T *tauq2, T *work, integer *lwork,
               integer *info)
    {
        orbdb(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, phi,
              taup1, taup2, tauq1, tauq2, work, lwork, info);
    }

    template <typename T, typename Ta>
    void unbdb(char *trans, char *signs, integer *m, integer *p, integer *q, T *x11, integer *ldx11,
               T *x12, integer *ldx12, T *x21, integer *ldx21, T *x22, integer *ldx22, Ta *theta,
               Ta *phi, T *taup1, T *taup2, T *tauq1, T *tauq2, T *work, integer *lwork,
               integer *info)
    {
        unbdb(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, phi,
              taup1, taup2, tauq1, tauq2, work, lwork, info);
    }
    /** @}*/ // end of orbdb

    /** @defgroup orbdb1 {un,or}bdb1
     * @ingroup CS
     * @{
     */
    /*! @brief ORBDB1 simultaneously bidiagonalizes the blocks of a tall and skinny \n
     matrix X with orthonomal columns
 * @details
 * \b Purpose:
    \verbatim
    ORBDB1 simultaneously bidiagonalizes the blocks of a tall and skinny
    matrix X with orthonomal columns:

                               [ B11 ]
         [ X11 ]   [ P1 |    ] [  0  ]
         [-----] = [---------] [-----] Q1**T .
         [ X21 ]   [    | P2 ] [ B21 ]
                               [  0  ]

    X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P,
    M-P, or M-Q. Routines SORBDB2, SORBDB3, and SORBDB4 handle cases in
    which Q is not the minimum dimension.

    The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
    and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
    Householder vectors.

    B11 and B12 are Q-by-Q bidiagonal matrices represented implicitly by
    angles THETA, PHI.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows X11 plus the number of rows in X21. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows in X11. 0 <= P <= M. \n
 * @param[in] Q
          Q is INTEGER \n
          The number of columns in X11 and X21. 0 <= Q <=
          MIN(P,M-P,M-Q). \n
 * @param[in,out] X11
          X11 is REAL array, dimension (LDX11,Q) \n
          On entry, the top block of the matrix X to be reduced. On
          exit, the columns of tril(X11) specify reflectors for P1 and
          the rows of triu(X11,1) specify reflectors for Q1. \n
 * @param[in] LDX11
          LDX11 is INTEGER \n
          The leading dimension of X11. LDX11 >= P. \n
 * @param[in,out] X21
          X21 is REAL array, dimension (LDX21,Q) \n
          On entry, the bottom block of the matrix X to be reduced. On
          exit, the columns of tril(X21) specify reflectors for P2. \n
 * @param[in] LDX21
          LDX21 is INTEGER \n
          The leading dimension of X21. LDX21 >= M-P. \n
 * @param[out] THETA
          THETA is REAL array, dimension (Q) \n
          The entries of the bidiagonal blocks B11, B21 are defined by
          THETA and PHI. See Further Details. \n
 * @param[out] PHI
          PHI is REAL array, dimension (Q-1) \n
          The entries of the bidiagonal blocks B11, B21 are defined by
          THETA and PHI. See Further Details. \n
 * @param[out] TAUP1
          TAUP1 is REAL array, dimension (P) \n
          The scalar factors of the elementary reflectors that define
          P1. \n
 * @param[out] TAUP2
          TAUP2 is REAL array, dimension (M-P) \n
          The scalar factors of the elementary reflectors that define
          P2. \n
 * @param[out] TAUQ1
          TAUQ1 is REAL array, dimension (Q) \n
          The scalar factors of the elementary reflectors that define
          Q1. \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= M-Q.
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array,   returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void orbdb1(integer *m, integer *p, integer *q, T *x11, integer *ldx11, T *x21, integer *ldx21,
                T *theta, T *phi, T *taup1, T *taup2, T *tauq1, T *work, integer *lwork,
                integer *info)
    {
        orbdb1(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
    }
    template <typename T, typename Ta>
    void unbdb1(integer *m, integer *p, integer *q, T *x11, integer *ldx11, T *x21, integer *ldx21,
                Ta *theta, Ta *phi, T *taup1, T *taup2, T *tauq1, T *work, integer *lwork,
                integer *info)
    {
        unbdb1(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
    }
    /** @}*/ // end of orbdb1

    /** @defgroup orbdb2 {un,or}bdb2
     * @ingroup CS
     * @{
     */
    /*! @brief ORBDB2 simultaneously bidiagonalizes the blocks of a tall and skinny  \n
     matrix X with orthonomal columns
 * @details
 * \b Purpose:
    \verbatim
    ORBDB2 simultaneously bidiagonalizes the blocks of a tall and skinny
    matrix X with orthonomal columns:

                               [ B11 ]
         [ X11 ]   [ P1 |    ] [  0  ]
         [-----] = [---------] [-----] Q1**T .
         [ X21 ]   [    | P2 ] [ B21 ]
                               [  0  ]

    X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P,
    Q, or M-Q. Routines SORBDB1, SORBDB3, and SORBDB4 handle cases in
    which P is not the minimum dimension.

    The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
    and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
    Householder vectors.

    B11 and B12 are P-by-P bidiagonal matrices represented implicitly by
    angles THETA, PHI.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows X11 plus the number of rows in X21. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows in X11. 0 <= P <= min(M-P,Q,M-Q). \n
 * @param[in] Q
          Q is INTEGER \n
          The number of columns in X11 and X21. 0 <= Q <= M. \n
 * @param[in,out] X11
          X11 is REAL array, dimension (LDX11,Q) \n
          On entry, the top block of the matrix X to be reduced. On
          exit, the columns of tril(X11) specify reflectors for P1 and
          the rows of triu(X11,1) specify reflectors for Q1. \n
 * @param[in] LDX11
          LDX11 is INTEGER \n
          The leading dimension of X11. LDX11 >= P. \n
 * @param[in,out] X21
          X21 is REAL array, dimension (LDX21,Q) \n
          On entry, the bottom block of the matrix X to be reduced. On
          exit, the columns of tril(X21) specify reflectors for P2. \n
 * @param[in] LDX21
          LDX21 is INTEGER \n
          The leading dimension of X21. LDX21 >= M-P. \n
 * @param[out] THETA
          THETA is REAL array, dimension (Q) \n
          The entries of the bidiagonal blocks B11, B21 are defined by
          THETA and PHI. See Further Details. \n
 * @param[out] PHI
          PHI is REAL array, dimension (Q-1) \n
          The entries of the bidiagonal blocks B11, B21 are defined by
          THETA and PHI. See Further Details. \n
 * @param[out] TAUP1
          TAUP1 is REAL array, dimension (P) \n
          The scalar factors of the elementary reflectors that define
          P1. \n
 * @param[out] TAUP2
          TAUP2 is REAL array, dimension (M-P) \n
          The scalar factors of the elementary reflectors that define
          P2. \n
 * @param[out] TAUQ1
          TAUQ1 is REAL array, dimension (Q) \n
          The scalar factors of the elementary reflectors that define
          Q1. \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= M-Q.
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array,   returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value.  \n

 *  * */
    template <typename T>
    void orbdb2(integer *m, integer *p, integer *q, T *x11, integer *ldx11, T *x21, integer *ldx21,
                T *theta, T *phi, T *taup1, T *taup2, T *tauq1, T *work, integer *lwork,
                integer *info)
    {
        orbdb2(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
    }
    template <typename T, typename Ta>
    void unbdb2(integer *m, integer *p, integer *q, T *x11, integer *ldx11, T *x21, integer *ldx21,
                Ta *theta, Ta *phi, T *taup1, T *taup2, T *tauq1, T *work, integer *lwork,
                integer *info)
    {
        unbdb2(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
    }

    /** @}*/ // end of orbdb2

    /** @defgroup orbdb3 {un,or}bdb3
     * @ingroup CS
     * @{
     */
    /*! @brief ORBDB3 simultaneously bidiagonalizes the blocks of a tall and skinny  \n
     matrix X with orthonomal columns
 * @details
 * \b Purpose:
    \verbatim
    ORBDB3 simultaneously bidiagonalizes the blocks of a tall and skinny
    matrix X with orthonomal columns:

                               [ B11 ]
         [ X11 ]   [ P1 |    ] [  0  ]
         [-----] = [---------] [-----] Q1**T .
         [ X21 ]   [    | P2 ] [ B21 ]
                               [  0  ]

    X11 is P-by-Q, and X21 is (M-P)-by-Q. M-P must be no larger than P,
    Q, or M-Q. Routines SORBDB1, SORBDB2, and SORBDB4 handle cases in
    which M-P is not the minimum dimension.

    The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
    and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
    Householder vectors.

    B11 and B12 are (M-P)-by-(M-P) bidiagonal matrices represented
    implicitly by angles THETA, PHI.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows X11 plus the number of rows in X21. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows in X11. 0 <= P <= M. M-P <= min(P,Q,M-Q). \n
 * @param[in] Q
          Q is INTEGER \n
          The number of columns in X11 and X21. 0 <= Q <= M. \n
 * @param[in,out] X11
          X11 is REAL array, dimension (LDX11,Q) \n
          On entry, the top block of the matrix X to be reduced. On
          exit, the columns of tril(X11) specify reflectors for P1 and
          the rows of triu(X11,1) specify reflectors for Q1. \n
 * @param[in] LDX11
          LDX11 is INTEGER \n
          The leading dimension of X11. LDX11 >= P. \n
 * @param[in,out] X21
          X21 is REAL array, dimension (LDX21,Q) \n
          On entry, the bottom block of the matrix X to be reduced. On
          exit, the columns of tril(X21) specify reflectors for P2. \n
 * @param[in] LDX21
          LDX21 is INTEGER \n
          The leading dimension of X21. LDX21 >= M-P. \n
 * @param[out] THETA
          THETA is REAL array, dimension (Q) \n
          The entries of the bidiagonal blocks B11, B21 are defined by
          THETA and PHI. See Further Details. \n
 * @param[out] PHI
          PHI is REAL array, dimension (Q-1) \n
          The entries of the bidiagonal blocks B11, B21 are defined by
          THETA and PHI. See Further Details. \n
 * @param[out] TAUP1
          TAUP1 is REAL array, dimension (P) \n
          The scalar factors of the elementary reflectors that define
          P1. \n
 * @param[out] TAUP2
          TAUP2 is REAL array, dimension (M-P) \n
          The scalar factors of the elementary reflectors that define
          P2. \n
 * @param[out] TAUQ1
          TAUQ1 is REAL array, dimension (Q) \n
          The scalar factors of the elementary reflectors that define
          Q1. \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= M-Q.
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array,   returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void orbdb3(integer *m, integer *p, integer *q, T *x11, integer *ldx11, T *x21, integer *ldx21,
                T *theta, T *phi, T *taup1, T *taup2, T *tauq1, T *work, integer *lwork,
                integer *info)
    {
        orbdb3(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
    }
    template <typename T, typename Ta>
    void unbdb3(integer *m, integer *p, integer *q, T *x11, integer *ldx11, T *x21, integer *ldx21,
                Ta *theta, Ta *phi, T *taup1, T *taup2, T *tauq1, T *work, integer *lwork,
                integer *info)
    {
        unbdb3(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
    }
    /** @}*/ // end of orbdb3

    /** @defgroup orbdb4 {un,or}bdb4
     * @ingroup CS
     * @{
     */
    /*! @brief ORBDB4 simultaneously bidiagonalizes the blocks of a tall and skinny \n
     matrix X with orthonomal columns
 * @details
 * \b Purpose:
    \verbatim
    ORBDB4 simultaneously bidiagonalizes the blocks of a tall and skinny
    matrix X with orthonomal columns:

                               [ B11 ]
         [ X11 ]   [ P1 |    ] [  0  ]
         [-----] = [---------] [-----] Q1**T .
         [ X21 ]   [    | P2 ] [ B21 ]
                               [  0  ]

    X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P,
    M-P, or Q. Routines SORBDB1, SORBDB2, and SORBDB3 handle cases in
    which M-Q is not the minimum dimension.

    The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
    and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
    Householder vectors.

    B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented
    implicitly by angles THETA, PHI.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows X11 plus the number of rows in X21. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows in X11. 0 <= P <= M. \n
 * @param[in] Q
          Q is INTEGER \n
          The number of columns in X11 and X21. 0 <= Q <= M and
          M-Q <= min(P,M-P,Q). \n
 * @param[in,out] X11
          X11 is REAL array, dimension (LDX11,Q) \n
          On entry, the top block of the matrix X to be reduced. On
          exit, the columns of tril(X11) specify reflectors for P1 and
          the rows of triu(X11,1) specify reflectors for Q1. \n
 * @param[in] LDX11
          LDX11 is INTEGER \n
          The leading dimension of X11. LDX11 >= P. \n
 * @param[in,out] X21
          X21 is REAL array, dimension (LDX21,Q) \n
          On entry, the bottom block of the matrix X to be reduced. On
          exit, the columns of tril(X21) specify reflectors for P2. \n
 * @param[in] LDX21
          LDX21 is INTEGER \n
          The leading dimension of X21. LDX21 >= M-P. \n
 * @param[out] THETA
          THETA is REAL array, dimension (Q) \n
          The entries of the bidiagonal blocks B11, B21 are defined by
          THETA and PHI. See Further Details. \n
 * @param[out] PHI
          PHI is REAL array, dimension (Q-1) \n
          The entries of the bidiagonal blocks B11, B21 are defined by
          THETA and PHI. See Further Details. \n
 * @param[out] TAUP1
          TAUP1 is REAL array, dimension (P) \n
          The scalar factors of the elementary reflectors that define
          P1. \n
 * @param[out] TAUP2
          TAUP2 is REAL array, dimension (M-P) \n
          The scalar factors of the elementary reflectors that define
          P2. \n
 * @param[out] TAUQ1
          TAUQ1 is REAL array, dimension (Q) \n
          The scalar factors of the elementary reflectors that define
          Q1. \n
 * @param[out] PHANTOM
          PHANTOM is REAL array, dimension (M) \n
          The routine computes an M-by-1 column vector Y that is
          orthogonal to the columns of [ X11; X21 ]. PHANTOM(1:P) and
          PHANTOM(P+1:M) contain Householder vectors for Y(1:P) and
          Y(P+1:M), respectively. \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= M-Q.
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array,   returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void orbdb4(integer *m, integer *p, integer *q, T *x11, integer *ldx11, T *x21, integer *ldx21,
                T *theta, T *phi, T *taup1, T *phantom, T *taup2, T *tauq1, T *work, integer *lwork,
                integer *info)
    {
        orbdb4(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work,
               lwork, info);
    }
    template <typename T, typename Ta>
    void unbdb4(integer *m, integer *p, integer *q, T *x11, integer *ldx11, T *x21, integer *ldx21,
                Ta *theta, Ta *phi, T *taup1, T *phantom, T *taup2, T *tauq1, T *work,
                integer *lwork, integer *info)
    {
        unbdb4(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work,
               lwork, info);
    }
    /** @}*/ // end of

    /** @defgroup unbdb5 {un,or}bdb5
     * @ingroup CS
     * @{
     */
    /*! @brief ORBDB5 orthogonalizes the column vector

 * @details
 * \b Purpose:
    \verbatim
    ORBDB5 orthogonalizes the column vector
         X = [ X1 ]
             [ X2 ]
    with respect to the columns of
         Q = [ Q1 ] .
             [ Q2 ]
    The columns of Q must be orthonormal.

    If the projection is zero according to Kahan's "twice is enough"
    criterion, then some other vector from the orthogonal complement
    is   returned. This vector is chosen in an arbitrary but deterministic
    way.
    \endverbatim

 * @param[in] M1
          M1 is INTEGER \n
          The dimension of X1 and the number of rows in Q1. 0 <= M1. \n
 * @param[in] M2
          M2 is INTEGER \n
          The dimension of X2 and the number of rows in Q2. 0 <= M2. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns in Q1 and Q2. 0 <= N. \n
 * @param[in,out] X1
          X1 is REAL array, dimension (M1) \n
          On entry, the top part of the vector to be orthogonalized.
          On exit, the top part of the projected vector. \n
 * @param[in] INCX1
          INCX1 is INTEGER \n
          Increment for entries of X1. \n
 * @param[in,out] X2
          X2 is REAL array, dimension (M2) \n
          On entry, the bottom part of the vector to be
          orthogonalized. On exit, the bottom part of the projected
          vector. \n
 * @param[in] INCX2
          INCX2 is INTEGER \n
          Increment for entries of X2. \n
 * @param[in] Q1
          Q1 is REAL array, dimension (LDQ1, N) \n
          The top part of the orthonormal basis matrix. \n
 * @param[in] LDQ1
          LDQ1 is INTEGER \n
          The leading dimension of Q1. LDQ1 >= M1. \n
 * @param[in] Q2
          Q2 is REAL array, dimension (LDQ2, N) \n
          The bottom part of the orthonormal basis matrix. \n
 * @param[in] LDQ2
          LDQ2 is INTEGER \n
          The leading dimension of Q2. LDQ2 >= M2. \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= N. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void orbdb5(integer *m1, integer *m2, integer *n, T *x1, integer *incx1, T *x2, integer *incx2,
                T *q1, integer *ldq1, T *q2, integer *ldq2, T *work, integer *lwork, integer *info)
    {
        orbdb5(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
    }
    template <typename T>
    void unbdb5(integer *m1, integer *m2, integer *n, T *x1, integer *incx1, T *x2, integer *incx2,
                T *q1, integer *ldq1, T *q2, integer *ldq2, T *work, integer *lwork, integer *info)
    {
        unbdb5(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
    }
    /** @}*/ // end of

    /** @defgroup orbdb6 {un,or}bdb6
     * @ingroup CS
     * @{
     */
    /*! @brief ORBDB6 orthogonalizes the column vector

 * @details
 * \b Purpose:
    \verbatim
    ORBDB6 orthogonalizes the column vector
         X = [ X1 ]
             [ X2 ]
    with respect to the columns of
         Q = [ Q1 ] .
             [ Q2 ]
    The columns of Q must be orthonormal.

    If the projection is zero according to Kahan's "twice is enough"
    criterion, then the zero vector is   returned.
    \endverbatim

 * @param[in] M1
          M1 is INTEGER \n
          The dimension of X1 and the number of rows in Q1. 0 <= M1. \n
 * @param[in] M2
          M2 is INTEGER \n
          The dimension of X2 and the number of rows in Q2. 0 <= M2. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns in Q1 and Q2. 0 <= N. \n
 * @param[in,out] X1
          X1 is REAL array, dimension (M1) \n
          On entry, the top part of the vector to be orthogonalized.
          On exit, the top part of the projected vector. \n
 * @param[in] INCX1
          INCX1 is INTEGER \n
          Increment for entries of X1. \n
 * @param[in,out] X2
          X2 is REAL array, dimension (M2) \n
          On entry, the bottom part of the vector to be
          orthogonalized. On exit, the bottom part of the projected
          vector. \n
 * @param[in] INCX2
          INCX2 is INTEGER \n
          Increment for entries of X2. \n
 * @param[in] Q1
          Q1 is REAL array, dimension (LDQ1, N) \n
          The top part of the orthonormal basis matrix. \n
 * @param[in] LDQ1
          LDQ1 is INTEGER \n
          The leading dimension of Q1. LDQ1 >= M1. \n
 * @param[in] Q2
          Q2 is REAL array, dimension (LDQ2, N) \n
          The bottom part of the orthonormal basis matrix. \n
 * @param[in] LDQ2
          LDQ2 is INTEGER \n
          The leading dimension of Q2. LDQ2 >= M2. \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= N. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void orbdb6(integer *m1, integer *m2, integer *n, T *x1, integer *incx1, T *x2, integer *incx2,
                T *q1, integer *ldq1, T *q2, integer *ldq2, T *work, integer *lwork, integer *info)
    {
        orbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
    }
    template <typename T>
    void unbdb6(integer *m1, integer *m2, integer *n, T *x1, integer *incx1, T *x2, integer *incx2,
                T *q1, integer *ldq1, T *q2, integer *ldq2, T *work, integer *lwork, integer *info)
    {
        unbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
    }
    /** @}*/ // end of orbdb6

    /** @defgroup lapmr lapmr
     * @ingroup CS
     * @{
     */
    /*! @brief LAPMR rearranges rows of a matrix as specified by a permutation vector.

 * @details
 * \b Purpose:
    \verbatim
     LAPMR rearranges the rows of the M by N matrix X as specified
     by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.
     If FORWRD = .TRUE.,  forward permutation:
          X(K(I),*) is moved X(I,*) for I = 1,2,...,M.

     If FORWRD = .FALSE., backward permutation:
          X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.
    \endverbatim

 * @param[in] FORWRD
          FORWRD is LOGICAL \n
          = .TRUE., forward permutation \n
          = .FALSE., backward permutation \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix X. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix X. N >= 0. \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,N) \n
          On entry, the M by N matrix X. \n
          On exit, X contains the permuted matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X, LDX >= MAX(1,M). \n
 * @param[in,out] K
          K is INTEGER array, dimension (M) \n
          On entry, K contains the permutation vector. K is used as
          internal workspace, but reset to its original value on
          output. \n

 * */
    template <typename T>
    void lapmr(logical *forwrd, integer *m, integer *n, T *x, integer *ldx, integer *k)
    {
        lapmr(forwrd, m, n, x, ldx, k);
    }
    /** @}*/ // end of lapmr

    /** @defgroup lapmt lapmt
     * @ingroup CS
     * @{
     */
    /*! @brief LAPMT performs a forward or backward permutation of the columns of a matrix.

 * @details
 * \b Purpose:
    \verbatim
     LAPMT rearranges the columns of the M by N matrix X as specified
     by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
     If FORWRD = .TRUE.,  forward permutation:
          X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.

     If FORWRD = .FALSE., backward permutation:
          X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.
    \endverbatim

 * @param[in] FORWRD
          FORWRD is LOGICAL \n
          = .TRUE., forward permutation \n
          = .FALSE., backward permutation \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix X. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix X. N >= 0. \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,N) \n
          On entry, the M by N matrix X. \n
          On exit, X contains the permuted matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X, LDX >= MAX(1,M). \n
 * @param[in,out] K
          K is INTEGER array, dimension (M) \n
          On entry, K contains the permutation vector. K is used as
          internal workspace, but reset to its original value on
          output. \n

 * */
    template <typename T>
    void lapmt(logical *forwrd, integer *m, integer *n, T *x, integer *ldx, integer *k)
    {
        lapmt(forwrd, m, n, x, ldx, k);
    }
    /** @}*/ // end of lapmt
    /** @} */ // end of CS

    /** @defgroup HH Household reflector
     * @ingroup Orthogonal
     * @{
     */
    /** @defgroup larf larf
     * @ingroup HH
     * @{
     */
    /*! @brief LARF applies an elementary reflector to a general rectangular matrix

 * @details
 * \b Purpose:
    \verbatim
    LARF applies a real elementary reflector H to a real m by n matrix
    C, from either the left or the right. H is represented in the form

          H = I - tau * v * v**T

    where tau is a real scalar and v is a real vector.

    If tau = 0, then H is taken to be the unit matrix.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': form  H * C \n
          = 'R': form  C * H \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. \n
 * @param[in] V
          V is REAL array, dimension \n
                     (1 + (M-1)*abs(INCV)) if SIDE = 'L' \n
                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R' \n
          The vector v in the representation of H. V is not used if
          TAU = 0. \n
 * @param[in] INCV
          INCV is INTEGER \n
          The increment between elements of v. INCV <> 0. \n
 * @param[in] TAU
          TAU is REAL \n
          The value tau in the representation of H. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the m by n matrix C. \n
          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
          or C * H if SIDE = 'R'. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out] WORK
          WORK is REAL array, dimension \n
                         (N) if SIDE = 'L' \n
                      or (M) if SIDE = 'R'  \n

 *  * */
    template <typename T>
    void larf(char *side, integer *m, integer *n, T *v, integer *incv, T *tau, T *c, integer *ldc,
              T *work)
    {
        larf(side, m, n, v, incv, tau, c, ldc, work);
    }
    /** @}*/ // end of larf

    /** @defgroup larf1f larf1f
     * @ingroup HH
     * @{
     */
    /*! @brief LARF1F applies an elementary reflector to a general rectangular matrix, \n
     assuming v(1) = 1
 * @details
 * \b Purpose:
    \verbatim
    LARF1F applies a real elementary reflector H to a real m by n matrix
    C, from either the left or the right. H is represented in the form

          H = I - tau * v * v**T

    where tau is a real scalar and v is a real vector.

    If tau = 0, then H is taken to be the unit matrix.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': form  H * C \n
          = 'R': form  C * H \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. \n
 * @param[in] V
          V is REAL array, dimension \n
                     (1 + (M-1)*abs(INCV)) if SIDE = 'L' \n
                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R' \n
          The vector v in the representation of H. V is not used if
          TAU = 0. V(1) is not referenced or modified. \n
 * @param[in] INCV
          INCV is INTEGER \n
          The increment between elements of v. INCV <> 0. \n
 * @param[in] TAU
          TAU is REAL \n
          The value tau in the representation of H. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the m by n matrix C. \n
          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
          or C * H if SIDE = 'R'. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out] WORK
          WORK is REAL array, dimension \n
                         (N) if SIDE = 'L' \n
                      or (M) if SIDE = 'R'  \n

 *  * */
    template <typename T>
    void larf1f(char *side, integer *m, integer *n, T *v, integer *incv, T *tau, T *c, integer *ldc,
                T *work)
    {
        larf1f(side, m, n, v, incv, tau, c, ldc, work);
    }
    /** @}*/ // end of larf1f

    /** @defgroup larf1l larf1l
     * @ingroup HH
     * @{
     */
    /*! @brief LARF1L applies an elementary reflector to a general rectangular matrix, \n
     assuming v(lastv) = 1
 * @details
 * \b Purpose:
    \verbatim
    LARF1L applies a real elementary reflector H to a real m by n matrix
    C, from either the left or the right. H is represented in the form

          H = I - tau * v * v**T

    where tau is a real scalar and v is a real vector.

    If tau = 0, then H is taken to be the unit matrix.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': form  H * C \n
          = 'R': form  C * H \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. \n
 * @param[in] V
          V is REAL array, dimension \n
                     (1 + (M-1)*abs(INCV)) if SIDE = 'L' \n
                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R' \n
          The vector v in the representation of H. V is not used if
          TAU = 0. \n
 * @param[in] INCV
          INCV is INTEGER \n
          The increment between elements of v. INCV <> 0. \n
 * @param[in] TAU
          TAU is REAL \n
          The value tau in the representation of H. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the m by n matrix C. \n
          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
          or C * H if SIDE = 'R'. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out] WORK
          WORK is REAL array, dimension \n
                         (N) if SIDE = 'L' \n
                      or (M) if SIDE = 'R'  \n

 *  * */
    template <typename T>
    void larf1l(char *side, integer *m, integer *n, T *v, integer *incv, T *tau, T *c, integer *ldc,
                T *work)
    {
        larf1l(side, m, n, v, incv, tau, c, ldc, work);
    }
    /** @}*/ // end of larf1l

    /** @defgroup larfx larfx
     * @ingroup HH
     * @{
     */
    /*! @brief LARFX applies an elementary reflector to a general rectangular matrix, \n
     with loop unrolling when the reflector has order ≤ 10
 * @details
 * \b Purpose:
    \verbatim
     LARFX applies a real elementary reflector H to a real m by n
     matrix C, from either the left or the right. H is represented in the
     form
           H = I - tau * v * v**T

     where tau is a real scalar and v is a real vector.

     If tau = 0, then H is taken to be the unit matrix

     This version uses inline code if H has order < 11.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': form  H * C \n
          = 'R': form  C * H \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. \n
 * @param[in] V
          V is REAL array, dimension (M) if SIDE = 'L'
                                     or (N) if SIDE = 'R' \n
          The vector v in the representation of H. \n
 * @param[in] TAU
          TAU is REAL \n
          The value tau in the representation of H. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the m by n matrix C. \n
          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
          or C * H if SIDE = 'R'. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= (1,M). \n
 * @param[out]	WORK
          WORK is REAL array, dimension
                      (N) if SIDE = 'L'
                      or (M) if SIDE = 'R' \n
          WORK is not referenced if H has order < 11. \n

 * */
    template <typename T>
    void larfx(char *side, integer *m, integer *n, T *v, T *tau, T *c, integer *ldc, T *work)
    {
        larfx(side, m, n, v, tau, c, ldc, work);
    }
    /** @}*/ // end of larfx

    /** @defgroup larfy larfy
     * @ingroup HH
     * @{
     */
    /*! @brief LARFY applies an elementary reflector, or Householder matrix, H, \n
     to an n x n symmetric matrix C, from both the left and the right
 * @details
 * \b Purpose:
    \verbatim
    LARFY applies an elementary reflector, or Householder matrix, H,
    to an n x n symmetric matrix C, from both the left and the right.

    H is represented in the form

       H = I - tau * v * v'

    where  tau  is a scalar and  v  is a vector.

    If  tau  is  zero, then  H  is taken to be the unit matrix.
    \endverbatim

  * @param[in] UPLO
           UPLO is CHARACTER*1 \n
           Specifies whether the upper or lower triangular part of the
           symmetric matrix C is stored. \n
           = 'U':  Upper triangle \n
           = 'L':  Lower triangle \n
  * @param[in] N
           N is INTEGER \n
           The number of rows and columns of the matrix C.  N >= 0. \n
  * @param[in] V
           V is REAL array, dimension
                   (1 + (N-1)*abs(INCV)) \n
           The vector v as described above. \n
  * @param[in] INCV
           INCV is INTEGER \n
           The increment between successive elements of v.  INCV must
           not be zero. \n
  * @param[in] TAU
           TAU is REAL \n
           The value tau as described above. \n
  * @param[in,out] C
           C is REAL array, dimension (LDC, N) \n
           On entry, the matrix C.
           On exit, C is overwritten by H * C * H'. \n
  * @param[in] LDC
           LDC is INTEGER \n
           The leading dimension of the array C.  LDC >= fla_max(1, N). \n
  * @param[out] WORK
           WORK is REAL array, dimension (N)   \n

 *  * */
    template <typename T>
    void larfy(char *uplo, integer *n, T *v, integer *incv, T *tau, T *c, integer *ldc, T *work)
    {
        larfy(uplo, n, v, incv, tau, c, ldc, work);
    }
    /** @}*/ // end of larfy

    /** @defgroup larfb larfb
     * @ingroup HH
     * @{
     */
    /*! @brief LARFB applies a block reflector or its transpose to a general rectangular matrix

 * @details
 * \b Purpose:
    \verbatim
     LARFB applies a real block reflector H or its transpose H**T to a
     real m by n matrix C, from either the left or the right.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply H or H**T from the Left \n
          = 'R': apply H or H**T from the Right \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N': apply H (No transpose) \n
          = 'T': apply H**T (Transpose) \n
 * @param[in] DIRECT
          DIRECT is CHARACTER*1 \n
          Indicates how H is formed from a product of elementary
          reflectors \n
          = 'F': H = H(1) H(2) . . . H(k) (Forward) \n
          = 'B': H = H(k) . . . H(2) H(1) (Backward) \n
 * @param[in] STOREV
          STOREV is CHARACTER*1 \n
          Indicates how the vectors which define the elementary
          reflectors are stored: \n
          = 'C': Columnwise \n
          = 'R': Rowwise \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. \n
 * @param[in] K
          K is INTEGER \n
          The order of the matrix T (= the number of elementary
          reflectors whose product defines the block reflector). \n
          If SIDE = 'L', M >= K >= 0; \n
          if SIDE = 'R', N >= K >= 0. \n
 * @param[in] V
          V is REAL array, dimension \n
                                (LDV,K) if STOREV = 'C' \n
                                (LDV,M) if STOREV = 'R' and SIDE = 'L' \n
                                (LDV,N) if STOREV = 'R' and SIDE = 'R' \n
          The matrix V. See Further Details. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V. \n
          If STOREV = 'C' and SIDE = 'L', LDV >= fla_max(1,M); \n
          if STOREV = 'C' and SIDE = 'R', LDV >= fla_max(1,N); \n
          if STOREV = 'R', LDV >= K. \n
 * @param[in] T
          T is REAL array, dimension (LDT,K) \n
          The triangular k by k matrix T in the representation of the
          block reflector. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= K. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the m by n matrix C. \n
          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (LDWORK,K) \n
 * @param[in]	LDWORK
          LDWORK is INTEGER \n
          The leading dimension of the array WORK. \n
          If SIDE = 'L', LDWORK >= fla_max(1,N); \n
          if SIDE = 'R', LDWORK >= fla_max(1,M). \n

   * */
    template <typename T>
    void larfb(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
               integer *k, T *v, integer *ldv, T *t, integer *ldt, T *c, integer *ldc, T *work,
               integer *ldwork)
    {
        larfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    }
    /** @}*/ // end of larfb

    /** @defgroup larfg larfg
     * @ingroup HH
     * @{
     */
    /*! @brief Generation of an elementary reflector (Householder matrix)
    *
    * @details
    * \b Purpose:
    * \verbatim
        Generation of a real elementary reflector (Householder matrix) H of order n, such that

        H * ( alpha) = ( beta),   H**T * H = I.
        (   X  )   (   0 )

        where alpha and beta are scalars, and X is an (n-1)-element real vector.
        H is represented in the form

        H = I - tau * ( 1) * ( 1 v**T) ,
        ( v)

        where tau is a real scalar and v is a real (n-1)-element vector.

        If the elements of X are all zero, then tau = 0 and H is taken to be the unit matrix.

        Otherwise  1 <= tau <= 2.
    \endverbatim

    * @param[in] n
              n is integer* \n
              The order of the elementary reflector. \n
    * @param[in,out] alpha
              alpha is float/double/COMPLEX/COMPLEX*16* \n
              On entry, the value alpha. \n
              On exit, it is overwritten with the value beta. \n
    * @param[in,out] x
              x is float/double/COMPLEX/COMPLEX*16 array, dimension (1+(N-2)*abs(incx)) \n
              On entry, the vector x. \n
              On exit, it is overwritten with the vector v. \n
    * @param[in] incx
              incx is integer* \n
              The increment between elements of X. incx > 0. \n
    * @param[out] tau
              tau is float/double/COMPLEX/COMPLEX*16* \n
              The value tau. \n

    *     *  */
    template <typename T>
    void larfg(integer *n, T *alpha, T *x, integer *incx, T *tau)
    {
        larfg(n, alpha, x, incx, tau);
    }
    /** @}*/ // end of larfg

    /** @defgroup larfgp larfgp
     * @ingroup HH
     * @{
     */
    /*! @brief Generation of an elementary reflector (Householder matrix) with non-negative beta
    *
    * @details
    * \b Purpose:
    * \verbatim
        Generation of a real elementary reflector H (Householder matrix) with non-negative beta
        of order n, such that

        H * ( alpha) = ( beta),   H**T * H = I.
        (   x  )   (   0 )

        where alpha and beta are scalars, beta is non-negative, and x is an (n-1)-element real
    vector. H is represented in the form

        H = I - tau * ( 1) * ( 1 v**T) ,
        ( v)

        where tau is a real scalar and v is a real (n-1)-element vector.

        If the elements of x are all zero, then tau = 0 and H is taken to be the unit matrix.
    \endverbatim

    * @param[in] n
              n is integer* \n
              The order of the elementary reflector. \n
    * @param[in,out] alpha
              alpha is float/double/COMPLEX/COMPLEX*16* \n
              On entry, the value alpha. \n
              On exit, it is overwritten with the value beta. \n
    * @param[in,out] x
              x is float/double/COMPLEX/COMPLEX*16 array, dimension (1+(N-2)*abs(incx)) \n
              On entry, the vector x. \n
              On exit, it is overwritten with the vector v. \n
    * @param[in] incx
              incx is integer* \n
              The increment between elements of X. incx > 0. \n
    * @param[out] tau
              tau is float/double/COMPLEX/COMPLEX*16* \n
              The value tau. \n

    *     *  */
    template <typename T>
    void larfgp(integer *n, T *alpha, T *x, integer *incx, T *tau)
    {
        larfgp(n, alpha, x, incx, tau);
    }
    /** @}*/ // end of larfgp

    /** @defgroup larft larft
     * @ingroup HH
     * @{
     */
    /*! @brief Triangular factor of a block reflector
    *
    * @details
    * \b Purpose:
    * \verbatim
        Formation of the triangular factor T of a real block reflector H of order n, which is
        defined as a product of k elementary reflectors.

        If direct = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;

        If direct = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.

        If storev = 'C', the vector which defines the elementary reflector H(i) is stored in the
        i-th column of the array V, and

        H  =  I - V * T * V**T

        If storev = 'R', the vector which defines the elementary reflector H(i) is stored in the
        i-th row of the array V, and

        H  =  I - V**T * T * V
    \endverbatim

    * @param[in] direct
              direct is char* \n
              Specifies the order in which the elementary reflectors are
              multiplied to form the block reflector: \n
              = 'F': H = H(1) H(2) . . . H(k) (Forward) \n
              = 'B': H = H(k) . . . H(2) H(1) (Backward) \n
    * @param[in] storev
              storev is char* \n
              Specifies how the vectors which define the elementary
              reflectors are stored (see also Further Details): \n
              = 'C': columnwise \n
              = 'R': rowwise \n
    * @param[in] n
              n is integer* \n
              The order of the block reflector H. n >= 0. \n
    * @param[in] k
              k is integer* \n
              The order of the triangular factor T (= the number of
              elementary reflectors). K >= 1. \n
    * @param[in] v
              v is float/double/COMPLEX/COMPLEX*16 array, dimension \n
              (ldv,k) if storev = 'C' \n
              (ldv,k) if storev = 'R' \n
              The matrix v. See further details. \n
    * @param[in] ldv
              ldv is integer* \n
              The leading dimension of the array V. \n
              If storev = 'C', ldv >= fla_max(1,n); if storev = 'R', ldv >= k. \n
    * @param[in] tau
              tau is float/double/COMPLEX/COMPLEX*16 array, dimension (k) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i). \n
    * @param[out] t
              t is float/double/COMPLEX/COMPLEX*16 array, dimension (ldt,k) \n
              The k by k triangular factor T of the block reflector. \n
              If direct = 'F', T is upper triangular; if direct = 'B', T is
              lower triangular. The rest of the array is not used. \n
    * @param[in] ldt
              ldt is integer* \n
              The leading dimension of the array T. ldt >= k.
              *
              * \n
              * **Further Details**
              * \verbatim
                      The shape of the matrix V and the storage of the vectors which define the
    H(i) is best illustrated by the following example with n = 5 and k = 3. The elements equal
    to 1 are not stored.

                      direct = 'F' and storev = 'C':   direct = 'F' and storev = 'R':

                      V = (  1)  V = (  1 v1 v1 v1 v1)
                      ( v1  1)   (  1 v2 v2 v2)
                      ( v1 v2  1)   (  1 v3 v3)
                      ( v1 v2 v3)
                      ( v1 v2 v3)

                      direct = 'B' and storev = 'C':   direct = 'B' and storev = 'R':

                      V = ( v1 v2 v3)  V = ( v1 v1  1)
                      ( v1 v2 v3)   ( v2 v2 v2  1)
                      (  1 v2 v3)   ( v3 v3 v3 v3  1)
                      (  1 v3)
                      (  1)
              \endverbatim

    *     *  */
    template <typename T>
    void larft(char *direct, char *storev, integer *n, integer *k, T *v, integer *ldv, T *tau, T *t,
               integer *ldt)
    {
        larft(direct, storev, n, k, v, ldv, tau, t, ldt);
    }
    /** @}*/ // end of larft
    /** @}*/ // end of HH

    /** @defgroup Jacob Givens/Jacobi plane rotations
     * @ingroup Orthogonal
     * @{
     */
    /** @defgroup lartg_jacobi lartg
     * @ingroup Jacobi
     * @{
     */
    /*! @brief LARTG generates a plane rotation with real cosine and real sine

 * @details
 * \b Purpose:
    \verbatim
    LARTG generate a plane rotation so that

       [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
       [ -SN  CS  ]     [ G ]     [ 0 ]

    This is a slower, more accurate version of the BLAS1 routine SROTG,
    with the following other differences:
       F and G are unchanged on   return.
       If G=0, then CS=1 and SN=0.
       If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
          floating point operations (saves work in SBDSQR when
          there are zeros on the diagonal).

    If F exceeds G in magnitude, CS will be positive.
    \endverbatim

 * @param[in] F
          F is REAL \n
          The first component of vector to be rotated. \n
 * @param[in] G
          G is REAL \n
          The second component of vector to be rotated. \n
 * @param[out] CS
          CS is REAL \n
          The cosine of the rotation. \n
 * @param[out] SN
          SN is REAL \n
          The sine of the rotation. \n
 * @param[out] R
          R is REAL \n
          The nonzero component of the rotated vector. \n
 \n
  This version has a few statements commented out for thread safety
  (machine parameters are computed on each entry). 10 feb 03, SJH.

 *  * */
    template <typename T>
    void lartg(T *f, T *g, T *cs, T *sn, T *r__)
    {
        lartg(f, g, cs, sn, r__);
    }
    template <typename T, typename Ta>
    void lartg(T *f, T *g, Ta *cs, T *sn, T *r__)
    {
        lartg(f, g, cs, sn, r__);
    }
    /** @}*/ // end of lartg_jacobi

    /** @defgroup lartgp lartgp
     * @ingroup Jacobi
     * @{
     */
    /*! @brief LARTGP generates a plane rotation so that the diagonal is nonnegative

 * @details
 * \b Purpose:
    \verbatim
     LARTGP generates a plane rotation so that

        [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
        [ -SN  CS  ]     [ G ]     [ 0 ]

     This is a slower, more accurate version of the Level 1 BLAS routine SROTG,
     with the following other differences:
        F and G are unchanged on return.
        If G=0, then CS=(+/-)1 and SN=0.
        If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1.

     The sign is chosen so that R >= 0.
    \endverbatim

 * @param[in] F
          F is REAL \n
          The first component of vector to be rotated. \n
 * @param[in] G
          G is REAL \n
          The second component of vector to be rotated. \n
 * @param[out] CS
          CS is REAL \n
          The cosine of the rotation. \n
 * @param[out] SN
          SN is REAL \n
          The sine of the rotation. \n
 * @param[out] R
          R is REAL \n
          The nonzero component of the rotated vector. \n

  This version has a few statements commented out for thread safety
  (machine parameters are computed on each entry). 10 feb 03, SJH. \n

 * */
    template <typename T>
    void lartgp(T *f, T *g, T *cs, T *sn, T *r)
    {
        lartgp(f, g, cs, sn, r);
    }
    /** @}*/ // end of lartgp

    /** @defgroup lasr lasr
     * @ingroup Jacobi
     * @{
     */
    /*! @brief LASR applies a sequence of plane rotations to a general rectangular matrix

 * @details
 * \b Purpose:
    \verbatim
    LASR applies a sequence of plane rotations to a real matrix A,
    from either the left or the right.

    When SIDE = 'L', the transformation takes the form

       A := P*A

    and when SIDE = 'R', the transformation takes the form

       A := A*P**T

    where P is an orthogonal matrix consisting of a sequence of z plane
    rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
    and P**T is the transpose of P.

    When DIRECT = 'F' (Forward sequence), then

       P = P(z-1) * ... * P(2) * P(1)

    and when DIRECT = 'B' (Backward sequence), then

       P = P(1) * P(2) * ... * P(z-1)

    where P(k) is a plane rotation matrix defined by the 2-by-2 rotation

       R(k) = ( c(k)  s(k))
            = (-s(k)  c(k)).

    When PIVOT = 'V' (Variable pivot), the rotation is performed
    for the plane (k,k+1), i.e., P(k) has the form

       P(k) = ( 1                                           )
              (      ...                                    )
              (             1                               )
              (                  c(k)  s(k)                 )
              (                 -s(k)  c(k)                 )
              (                               1             )
              (                                    ...      )
              (                                           1 )

    where R(k) appears as a rank-2 modification to the identity matrix in
    rows and columns k and k+1.

    When PIVOT = 'T' (Top pivot), the rotation is performed for the
    plane (1,k+1), so P(k) has the form

       P(k) = ( c(k)                    s(k)                )
              (        1                                    )
              (             ...                             )
              (                    1                        )
              (-s(k)                    c(k)                )
              (                                1            )
              (                                     ...     )
              (                                            1)

    where R(k) appears in rows and columns 1 and k+1.

    Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
    performed for the plane (k,z), giving P(k) the form

       P(k) = (1                                            )
              (     ...                                     )
              (            1                                )
              (                 c(k)                    s(k))
              (                        1                    )
              (                             ...             )
              (                                    1        )
              (                -s(k)                    c(k))

    where R(k) appears in rows and columns k and z.  The rotations are
    performed without ever forming P(k) explicitly.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          Specifies whether the plane rotation matrix P is applied to
          A on the left or the right. \n
          = 'L':  Left, compute A := P*A \n
          = 'R':  Right, compute A:= A*P**T \n
 * @param[in] PIVOT
          PIVOT is CHARACTER*1 \n
          Specifies the plane for which P(k) is a plane rotation
          matrix. \n
          = 'V':  Variable pivot, the plane (k,k+1) \n
          = 'T':  Top pivot, the plane (1,k+1) \n
          = 'B':  Bottom pivot, the plane (k,z) \n
 * @param[in] DIRECT
          DIRECT is CHARACTER*1 \n
          Specifies whether P is a forward or backward sequence of
          plane rotations. \n
          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1) \n
          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1) \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  If m <= 1, an immediate
          return is effected. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  If n <= 1, an
          immediate   return is effected. \n
 * @param[in] C
          C is REAL array, dimension \n
                  (M-1) if SIDE = 'L' \n
                  (N-1) if SIDE = 'R' \n
          The cosines c(k) of the plane rotations. \n
 * @param[in] S
          S is REAL array, dimension \n
                  (M-1) if SIDE = 'L' \n
                  (N-1) if SIDE = 'R' \n
          The sines s(k) of the plane rotations.  The 2-by-2 plane
          rotation part of the matrix P(k), R(k), has the form \n
          R(k) = ( c(k)  s(k)) \n
                 (-s(k)  c(k)). \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          The M-by-N matrix A.  On exit, A is overwritten by P*A if
          SIDE = 'R' or by A*P**T if SIDE = 'L'. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M).  \n

 *  * */
    template <typename T>
    void lasr(char *side, char *pivot, char *direct, integer *m, integer *n, T *c, T *s, T *a,
              integer *lda)
    {
        lasr(side, pivot, direct, m, n, c, s, a, lda);
    }
    template <typename T, typename Ta>
    void lasr(char *side, char *pivot, char *direct, integer *m, integer *n, Ta *c, Ta *s, T *a,
              integer *lda)
    {
        lasr(side, pivot, direct, m, n, c, s, a, lda);
    }
    /** @}*/ // end of lasr

    /** @defgroup largv largv
     * @ingroup Jacobi
     * @{
     */
    /*! @brief LARGV generates a vector of plane rotations with real cosines and real sines

 * @details
 * \b Purpose:
    \verbatim
    LARGV generates a vector of real plane rotations, determined by
    elements of the real vectors x and y. For i = 1,2,...,n

       ( c(i)  s(i)) (x(i)) = (a(i))
       (-s(i)  c(i)) (y(i)) = (  0 )
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of plane rotations to be generated. \n
 * @param[in,out] X
          X is REAL array,
                         dimension (1+(N-1)*INCX) \n
          On entry, the vector x.
          On exit, x(i) is overwritten by a(i), for i = 1,...,n. \n
 * @param[in] INCX
          INCX is INTEGER \n
          The increment between elements of X. INCX > 0. \n
 * @param[in,out] Y
          Y is REAL array,
                         dimension (1+(N-1)*INCY) \n
          On entry, the vector y.
          On exit, the sines of the plane rotations. \n
 * @param[in] INCY
          INCY is INTEGER \n
          The increment between elements of Y. INCY > 0. \n
 * @param[out] C
          C is REAL array, dimension (1+(N-1)*INCC) \n
          The cosines of the plane rotations. \n
 * @param[in] INCC
          INCC is INTEGER \n
          The increment between elements of C. INCC > 0. \n

 *  * */
    template <typename T>
    void largv(integer *n, T *x, integer *incx, T *y, integer *incy, T *c, integer *incc)
    {
        largv(n, x, incx, y, incy, c, incc);
    }
    template <typename T, typename Ta>
    void largv(integer *n, T *x, integer *incx, T *y, integer *incy, Ta *c, integer *incc)
    {
        largv(n, x, incx, y, incy, c, incc);
    }
    /** @}*/ // end of largv

    /** @defgroup lartv lartv
     * @ingroup Jacobi
     * @{
     */
    /*! @brief LARTV applies a vector of plane rotations with real cosines and \n
     real sines to the elements of a pair of vectors
 * @details
 * \b Purpose:
    \verbatim
    LARTV applies a vector of real plane rotations to elements of the
    real vectors x and y. For i = 1,2,...,n

       (x(i)) := ( c(i)  s(i)) (x(i))
       (y(i))    (-s(i)  c(i)) (y(i))
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of plane rotations to be applied. \n
 * @param[in,out] X
          X is REAL array,
                         dimension (1+(N-1)*INCX) \n
          The vector x. \n
 * @param[in] INCX
          INCX is INTEGER \n
          The increment between elements of X. INCX > 0. \n
 * @param[in,out] Y
          Y is REAL array,
                         dimension (1+(N-1)*INCY) \n
          The vector y. \n
 * @param[in] INCY
          INCY is INTEGER \n
          The increment between elements of Y. INCY > 0. \n
 * @param[in] C
          C is REAL array, dimension (1+(N-1)*INCC) \n
          The cosines of the plane rotations. \n
 * @param[in] S
          S is REAL array, dimension (1+(N-1)*INCC) \n
          The sines of the plane rotations. \n
 * @param[in] INCC
          INCC is INTEGER \n
          The increment between elements of C and S. INCC > 0. \n

 *  * */
    template <typename T>
    void lartv(integer *n, T *x, integer *incx, T *y, integer *incy, T *c, T *s, integer *incc)
    {
        lartv(n, x, incx, y, incy, c, s, incc);
    }
    template <typename T, typename Ta>
    void lartv(integer *n, T *x, integer *incx, T *y, integer *incy, Ta *c, T *s, integer *incc)
    {
        lartv(n, x, incx, y, incy, c, s, incc);
    }
    /** @}*/ // end of lartv

    /** @defgroup lar2v lar2v
     * @ingroup Jacobi
     * @{
     */
    /*! @brief LAR2V applies a vector of plane rotations with real cosines and real \n
     sines from both sides to a sequence of 2-by-2 symmetric/Hermitian matrices
 * @details
 * \b Purpose:
    \verbatim
    LAR2V applies a vector of real plane rotations from both sides to
    a sequence of 2-by-2 real symmetric matrices, defined by the elements
    of the vectors x, y and z. For i = 1,2,...,n

    (x(i)  z(i)) := ( c(i)  s(i)) (x(i)  z(i)) (c(i) -s(i))
    (z(i)  y(i))    (-s(i)  c(i)) (z(i)  y(i)) (s(i)  c(i))
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of plane rotations to be applied. \n
 * @param[in,out] X
          X is REAL array,
                         dimension (1+(N-1)*INCX) \n
          The vector x. \n
 * @param[in,out] Y
          Y is REAL array,
                         dimension (1+(N-1)*INCX) \n
          The vector y. \n
 * @param[in,out] Z
          Z is REAL array,
                         dimension (1+(N-1)*INCX) \n
          The vector z. \n
 * @param[in] INCX
          INCX is INTEGER \n
          The increment between elements of X, Y and Z. INCX > 0. \n
 * @param[in] C
          C is REAL array, dimension (1+(N-1)*INCC) \n
          The cosines of the plane rotations. \n
 * @param[in] S
          S is REAL array, dimension (1+(N-1)*INCC) \n
          The sines of the plane rotations. \n
 * @param[in] INCC
          INCC is INTEGER \n
          The increment between elements of C and S. INCC > 0.  \n

 *  * */
    template <typename T>
    void lar2v(integer *n, T *x, T *y, T *z, integer *incx, T *c, T *s, integer *incc)
    {
        lar2v(n, x, y, z, incx, c, s, incc);
    }
    template <typename T, typename Ta>
    void lar2v(integer *n, T *x, T *y, T *z, integer *incx, Ta *c, T *s, integer *incc)
    {
        lar2v(n, x, y, z, incx, c, s, incc);
    }

    /** @}*/ // end of lar2v

    /** @defgroup  lacrt lacrt
     * @ingroup Jacobi
     * @{
     */
    /*! @brief LACRT performs a linear transformation of a pair of complex vectors

 * @details
 * \b Purpose:
    \verbatim
    LACRT performs the operation

       ( c  s)(x)  ==> (x)
       (-s  c)(y)      (y)

    where c and s are complex and the vectors x and y are complex.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of elements in the vectors CX and CY. \n
 * @param[in,out] CX
          CX is COMPLEX array, dimension (N) \n
          On input, the vector x.
          On output, CX is overwritten with c*x + s*y. \n
 * @param[in] INCX
          INCX is INTEGER \n
          The increment between successive values of CX.  INCX <> 0. \n
 * @param[in,out] CY
          CY is COMPLEX array, dimension (N) \n
          On input, the vector y.
          On output, CY is overwritten with -s*x + c*y. \n
 * @param[in] INCY
          INCY is INTEGER \n
          The increment between successive values of CY.  INCY <> 0. \n
 * @param[in] C
          C is COMPLEX \n
 * @param[in] S
          S is COMPLEX \n
          C and S define the matrix \n
             [  C   S  ]. \n
             [ -S   C  ] \n

 * */
    template <typename T>
    void lacrt(integer *n, T *cx, integer *incx, T *cy, integer *incy, T *c, T *s)
    {
        lacrt(n, cx, incx, cy, incy, c, s);
    }
    /** @}*/ // end of lacr
    /** @}*/ // end of Jacobi
    /** @} */ // end of Orthogonal
}
#endif