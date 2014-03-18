
#include "FLAME.h"

FLA_Error FLA_Apply_Q_UT_rnfc_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
/*
  Apply a unitary matrix Q to a matrix B from the right,

    B :=  B Q

  where Q is the forward product of Householder transformations:

    Q  =  H(0) H(1) ... H(k-1)

  where H(i) corresponds to the Householder vector stored below the diagonal
  in the ith column of A. Thus, the operation becomes:

    B :=  B Q
       =  B H(0) H(1) ... H(k-1)

  From this, we can see that we must move through A from top-left to bottom-
  right, since the Householder vector for H(0) was stored in the first column
  of A. We intend to apply blocks of reflectors at a time, where a block
  reflector H of b consecutive Householder transforms may be expressed as:

    H  =  ( H(i) H(i+1) ... H(i+b-1) )
       =  ( I - U inv(T) U' )

  where:
    - U is the strictly lower trapezoidal (with implicit unit diagonal) matrix
      of Householder vectors, stored below the diagonal of A in columns i through
      i+b-1, corresponding to H(i) through H(i+b-1).
    - T is the upper triangular block Householder matrix corresponding to
      Householder vectors i through i+b-1.

  Consider applying H to B as an intermediate step towards applying all of Q:

    B  :=  B H
        =  B ( I - U inv(T) U' )
        =  B - B U inv(T) U'

  We must move from top-left to bottom-right. So, we partition:

    U -> / U11 \  B -> ( B1 B2 )  T -> ( T1 T2 )
         \ U21 /
  where:
    - U11 is stored in strictly lower triangle of A11 with implicit unit
      diagonal.
    - U21 is stored in A21.
    - T1 is an upper triangular block of row-panel matrix T.

  Substituting repartitioned U, B, and T, we have:

    ( B1 B2 ) :=  ( B1 B2 ) - ( B1 B2 ) / U11 \ inv(T1) / U11 \'
                                        \ U21 /         \ U21 /
               =  ( B1 B2 ) - ( B1 B2 ) / U11 \ inv(T1) ( U11' U21' )
                                        \ U21 /
               =  ( B1 B2 ) - ( B1 U11 + B2 U21 ) inv(T1) ( U11' U21' )

  Thus, B1 is updated as:

      B1    :=      B1   -   ( B1 U11 + B2 U21 ) inv(T1) U11'

  And B2 is updated as:

      B2    :=      B2   -   ( B1 U11 + B2 U21 ) inv(T1) U21'

  Note that:

    ( B1 U11 + B2 U21 ) inv(T1)

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

  FLA_Obj BL,    BR,       B0,  B1,  B2;

  dim_t   b_alg, b;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    b = min( b_alg, FLA_Obj_min_dim( ABR ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &T2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,    &T1T, 
                         &T2B,     b, FLA_TOP );

    FLA_Part_2x2( W,     &WTL, &WTR,
                         &WBL, &WBR,     b, FLA_Obj_length( B1 ), FLA_TL );

    // WTL = B1^T;

    FLA_Copyt_internal( FLA_TRANSPOSE, B1, WTL,
                        FLA_Cntl_sub_copyt( cntl ) );

    // U11 = trilu( A11 );
    // U21 = A21;
    // Let WTL^T be conformal to B1.
    //
    // WTL^T = ( B1 * U11 + B2 * U21 ) * inv( triu(T1T) );
    // WTL   = inv( triu(T1T)^T ) * ( U11^T * B1^T + U21^T * B2^T );

    FLA_Trmm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                       FLA_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm1( cntl ) );

    FLA_Gemm_internal( FLA_TRANSPOSE, FLA_TRANSPOSE, 
                       FLA_ONE, A21, B2, FLA_ONE, WTL,
                       FLA_Cntl_sub_gemm1( cntl ) );

    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T1T, WTL,
                       FLA_Cntl_sub_trsm( cntl ) );

    // B2 = B2 - WTL^T * U21';
    // B1 = B1 - WTL^T * U11';
    //    = B1 - ( conj(U11) * WTL )^T;

    FLA_Gemm_internal( FLA_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_MINUS_ONE, WTL, A21, FLA_ONE, B2,
                       FLA_Cntl_sub_gemm2( cntl ) );

    FLA_Trmm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                       FLA_CONJ_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_MINUS_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm2( cntl ) );

    FLA_Axpyt_internal( FLA_TRANSPOSE, FLA_ONE, WTL, B1,
                        FLA_Cntl_sub_axpyt( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ T2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

