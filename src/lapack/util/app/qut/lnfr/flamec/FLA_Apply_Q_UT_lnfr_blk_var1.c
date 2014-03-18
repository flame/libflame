
#include "FLAME.h"

FLA_Error FLA_Apply_Q_UT_lnfr_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
/*
  Apply a unitary matrix Q to a matrix B from the left,

    B :=  Q B

  where Q is the forward product of Householder transformations:

    Q  =  H(0) H(1) ... H(k-1)

  where H(i) corresponds to the Householder vector stored above the diagonal
  in the ith row of A. Thus, the operation becomes:

    B :=  Q B
       =  H(0) H(1) ... H(k-1) B

  From this, we can see that we must move through A from bottom-right to top-
  left, since the Householder vector for H(k-1) was stored in the last row
  of A. We intend to apply blocks of reflectors at a time, where a block
  reflector H of b consecutive Householder transforms may be expressed as:

    H  =  ( H(i) H(i+1) ... H(i+b-1) )
       =  ( I - U inv(T) U' )

  where:
    - U^T is the strictly upper trapezoidal (with implicit unit diagonal) matrix
      of Householder vectors, stored above the diagonal of A in rows i through
      i+b-1, corresponding to H(i) through H(i+b-1).
    - T is the upper triangular block Householder matrix corresponding to
      Householder vectors i through i+b-1.

  Consider applying H to B as an intermediate step towards applying all of Q:

    B  :=  H B
        =  ( I - U inv(T) U' ) B
        =  B - U inv(T) U' B

  We must move from bottom-right to top-left. So, we partition:

    U^T -> ( U11 U12 )  B -> / B1 \  T -> ( T2 T1 )
                             \ B2 / 

  where:
    - U11 is stored in strictly upper triangle of A11 with implicit unit
      diagonal.
    - U12 is stored in A12.
    - T1 is an upper triangular block of row-panel matrix T.

  Substituting repartitioned U, B, and T, we have:

    / B1 \  :=   / B1 \ - ( U11 U12 )^T inv(T1) conj( U11 U12 ) / B1 \
    \ B2 /       \ B2 /                                         \ B2 /
             =   / B1 \ - / U11^T \ inv(T1) conj( U11 U12 ) / B1 \
                 \ B2 /   \ U12^T /                         \ B2 /
             =   / B1 \ - / U11^T \ inv(T1) ( conj(U11) B1 + conj(U12) B2 )
                 \ B2 /   \ U12^T /

  Thus, B1 is updated as:

      B1    :=     B1   -   U11^T inv(T1) ( conj(U11) B1 + conj(U12) B2 )

  And B2 is updated as:

      B2    :=     B2   -   U12^T inv(T1) ( conj(U11) B1 + conj(U12) B2 )

  Note that:

    inv(T1) ( conj(U11) B1 + conj(U12) B2 )

  is common to both updates, and thus may be computed and stored in
  workspace, and then re-used.

  -FGVZ
*/
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj T1T,
          T2B;

  FLA_Obj WTL,  WTR,
          WBL,  WBR;

  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  dim_t   b_alg, b;
  dim_t   m_BR, n_BR;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  // If m < n, then we have to initialize our partitionings carefully so
  // that we begin in the proper location in A and B (since we traverse
  // matrix A from BR to TL).
  if ( FLA_Obj_length( A ) < FLA_Obj_width( A ) )
  {
    m_BR = 0;
    n_BR = FLA_Obj_width( A ) - FLA_Obj_length( A );
  }
  else if ( FLA_Obj_length( A ) > FLA_Obj_width( A ) )
  {
    m_BR = FLA_Obj_length( A ) - FLA_Obj_width( A );
    n_BR = 0;
  }
  else
  {
    m_BR = 0;
    n_BR = 0;
  }

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     m_BR, n_BR, FLA_BR );

  // A and T are dependent; we determine T matrix w.r.t. A
  FLA_Part_1x2( T,    &TL,  &TR,      FLA_Obj_min_dim( A ), FLA_LEFT );

  // Be carefule that A contains reflector in row-wise; 
  // corresponding B should be partitioned with n_BR.
  FLA_Part_2x1( B,    &BT, 
                      &BB,            n_BR, FLA_BOTTOM );

  while ( FLA_Obj_min_dim( ATL ) > 0 ){

    b = min( b_alg, FLA_Obj_min_dim( ATL ) );

    // Since T was filled from left to right, and since we need to access them
    // in reverse order, we need to handle the case where the last block is
    // smaller than the other b x b blocks.
    if ( FLA_Obj_width( TR ) == 0 && FLA_Obj_width( T ) % b_alg > 0 )
      b = FLA_Obj_width( T ) % b_alg;

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, &A01, /**/ &A02,
                                                &A10, &A11, /**/ &A12,
                        /* ************* */   /* ******************** */
                           ABL, /**/ ABR,       &A20, &A21, /**/ &A22,
                           b, b, FLA_TL );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, &T1, /**/ &T2,
                           b, FLA_LEFT );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                                              &B1, 
                        /* ** */            /* ** */
                           BB,                &B2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,    &T1T, 
                         &T2B,     b, FLA_TOP );

    FLA_Part_2x2( W,     &WTL, &WTR,
                         &WBL, &WBR,     b, FLA_Obj_width( B1 ), FLA_TL );

    // WTL = B1;

    FLA_Copyt_internal( FLA_NO_TRANSPOSE, B1, WTL,
                        FLA_Cntl_sub_copyt( cntl ) );

    // U11 = triuu( A11 );
    // U12 = A12;
    // 
    // WTL = inv( triu(T1T) ) * ( conj(U11) * B1 + conj(U12) * B2 );

    FLA_Trmm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_CONJ_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm1( cntl ) );

    FLA_Gemm_internal( FLA_CONJ_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_ONE, A12, B2, FLA_ONE, WTL,
                       FLA_Cntl_sub_gemm1( cntl ) );

    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T1T, WTL,
                       FLA_Cntl_sub_trsm( cntl ) );

    // B2 = B2 - U12^T * WTL;
    // B1 = B1 - U11^T * WTL;

    FLA_Gemm_internal( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A12, WTL, FLA_ONE, B2,
                       FLA_Cntl_sub_gemm2( cntl ) );

    FLA_Trmm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_MINUS_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm2( cntl ) );

    FLA_Axpyt_internal( FLA_NO_TRANSPOSE, FLA_ONE, WTL, B1,
                        FLA_Cntl_sub_axpyt( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                            /* ************** */  /* ****************** */
                                                     A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                              FLA_BR );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, /**/ T1, T2,
                              FLA_RIGHT );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                            /* ** */           /* ** */
                                                  B1, 
                              &BB,                B2,     FLA_BOTTOM );
  }

  return FLA_SUCCESS;
}

