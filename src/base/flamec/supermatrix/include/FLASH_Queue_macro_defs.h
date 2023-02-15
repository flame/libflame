/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef FLASH_QUEUE_MACRO_DEFS_H
#define FLASH_QUEUE_MACRO_DEFS_H


#ifdef FLA_ENABLE_SUPERMATRIX


#define FLASH_OBJ_PTR_ID( A )  ( A ).base->id

// FLASH_Verbose
#define FLASH_QUEUE_VERBOSE_NONE                     0
#define FLASH_QUEUE_VERBOSE_READABLE                 1
#define FLASH_QUEUE_VERBOSE_GRAPHVIZ                 2

// FLASH_Data_aff
#define FLASH_QUEUE_AFFINITY_NONE                    0
#define FLASH_QUEUE_AFFINITY_2D_BLOCK_CYCLIC         1
#define FLASH_QUEUE_AFFINITY_1D_ROW_BLOCK_CYCLIC     2
#define FLASH_QUEUE_AFFINITY_1D_COLUMN_BLOCK_CYCLIC  3
#define FLASH_QUEUE_AFFINITY_ROUND_ROBIN             4

/*
Reminder to create a macro to enqueue when SuperMatrix is configured, and
also to create a macro for when it is not below to return an error code.
*/

// LAPACK-level

#define ENQUEUE_FLASH_LU_piv_macro( A, p, cntl ) \
        FLASH_Queue_push( (void *) FLA_LU_piv_macro_task, \
                          (void *) cntl, \
                          "LU   ", \
                          FALSE, \
                          FALSE, \
                          0, 0, 0, 2, \
                          A, p )

#define ENQUEUE_FLASH_Apply_pivots_macro( side, trans, p, A, cntl ) \
        FLASH_Queue_push( (void *) FLA_Apply_pivots_macro_task, \
                          (void *) cntl, \
                          "Pivot", \
                          FALSE, \
                          FALSE, \
                          2, 0, 1, 1, \
                          side, trans, \
                          p, A )

#define ENQUEUE_FLASH_LU_piv( A, p, cntl ) \
        FLASH_Queue_push( (void *) FLA_LU_piv_task, \
                          (void *) cntl, \
                          "LU   ", \
                          FALSE, \
                          TRUE, \
                          0, 1, 0, 1, \
                          p, A )

#define ENQUEUE_FLASH_LU_piv_copy( A, p, U, cntl ) \
        FLASH_Queue_push( (void *) FLA_LU_piv_copy_task, \
                          (void *) cntl, \
                          "LU   ", \
                          FALSE, \
                          TRUE, \
                          0, 1, 0, 2, \
                          p, A, U )

#define ENQUEUE_FLASH_Trsm_piv( A, C, p, cntl ) \
        FLASH_Queue_push( (void *) FLA_Trsm_piv_task, \
                          (void *) cntl, \
                          "Trsm ", \
                          FALSE, \
                          TRUE, \
                          0, 1, 1, 1, \
                          p, A, C )

#define ENQUEUE_FLASH_SA_LU( U, D, p, L, nb_alg, cntl ) \
        FLASH_Queue_push( (void *) FLA_SA_LU_task, \
                          (void *) cntl, \
                          "SA_LU", \
                          FALSE, \
                          FALSE, \
                          1, 2, 0, 2, \
                          nb_alg, \
                          p, L, D, U )

#define ENQUEUE_FLASH_SA_FS( L, D, p, C, E, nb_alg, cntl ) \
        FLASH_Queue_push( (void *) FLA_SA_FS_task, \
                          (void *) cntl, \
                          "SA_FS", \
                          FALSE, \
                          FALSE, \
                          1, 2, 1, 2, \
                          nb_alg, \
                          L, p, D, E, C )

#define ENQUEUE_FLASH_LU_nopiv( A, cntl ) \
        FLASH_Queue_push( (void *) FLA_LU_nopiv_task, \
                          (void *) cntl, \
                          "LU   ", \
                          FALSE, \
                          FALSE, \
                          0, 0, 0, 1, \
                          A )

#define ENQUEUE_FLASH_Trinv( uplo, diag, A, cntl ) \
        FLASH_Queue_push( (void *) FLA_Trinv_task, \
                          (void *) cntl, \
                          "Trinv", \
                          FALSE, \
                          TRUE, \
                          2, 0, 0, 1, \
                          uplo, diag, \
                          A )

#define ENQUEUE_FLASH_Ttmm( uplo, A, cntl ) \
        FLASH_Queue_push( (void *) FLA_Ttmm_task, \
                          (void *) cntl, \
                          "Ttmm ", \
                          FALSE, \
                          FALSE, \
                          1, 0, 0, 1, \
                          uplo, \
                          A )

#define ENQUEUE_FLASH_Chol( uplo, A, cntl ) \
        FLASH_Queue_push( (void *) FLA_Chol_task, \
                          (void *) cntl, \
                          "Chol ", \
                          FALSE, \
                          TRUE, \
                          1, 0, 0, 1, \
                          uplo, \
                          A )

#define ENQUEUE_FLASH_Sylv( transA, transB, isgn, A, B, C, scale, cntl ) \
        FLASH_Queue_push( (void *) FLA_Sylv_task, \
                          (void *) cntl, \
                          "Sylv ", \
                          FALSE, \
                          FALSE, \
                          2, 2, 2, 1, \
                          transA, transB, \
                          isgn, scale, \
                          A, B, C )

#define ENQUEUE_FLASH_Lyap( trans, isgn, A, C, scale, cntl ) \
        FLASH_Queue_push( (void *) FLA_Lyap_task, \
                          (void *) cntl, \
                          "Lyap ", \
                          FALSE, \
                          FALSE, \
                          1, 2, 1, 1, \
                          trans, \
                          isgn, scale, \
                          A, C )

#define ENQUEUE_FLASH_QR_UT_macro( A, T, cntl ) \
        FLASH_Queue_push( (void *) FLA_QR_UT_macro_task, \
                          (void *) cntl, \
                          "QR   ", \
                          FALSE, \
                          FALSE, \
                          0, 0, 0, 2, \
                          A, T )

#define ENQUEUE_FLASH_QR_UT( A, T, cntl ) \
        FLASH_Queue_push( (void *) FLA_QR_UT_task, \
                          (void *) cntl, \
                          "QR   ", \
                          FALSE, \
                          FALSE, \
                          0, 1, 0, 1, \
                          T, A )

#define ENQUEUE_FLASH_QR_UT_copy( A, T, U, cntl ) \
        FLASH_Queue_push( (void *) FLA_QR_UT_copy_task, \
                          (void *) cntl, \
                          "QR   ", \
                          FALSE, \
                          FALSE, \
                          0, 1, 0, 2, \
                          T, A, U )

#define ENQUEUE_FLASH_QR2_UT( B, D, T, cntl ) \
        FLASH_Queue_push( (void *) FLA_QR2_UT_task, \
                          (void *) cntl, \
                          "QR2  ", \
                          FALSE, \
                          FALSE, \
                          0, 1, 0, 2, \
                          T, D, B )

#define ENQUEUE_FLASH_LQ_UT_macro( A, T, cntl ) \
        FLASH_Queue_push( (void *) FLA_LQ_UT_macro_task, \
                          (void *) cntl, \
                          "LQ   ", \
                          FALSE, \
                          FALSE, \
                          0, 0, 0, 2, \
                          A, T )

#define ENQUEUE_FLASH_CAQR2_UT( B, D, T, cntl ) \
        FLASH_Queue_push( (void *) FLA_CAQR2_UT_task, \
                          (void *) cntl, \
                          "CAQR2", \
                          FALSE, \
                          FALSE, \
                          0, 1, 0, 2, \
                          T, D, B )

#define ENQUEUE_FLASH_Apply_Q_UT( side, trans, direct, storev, A, T, W, B, cntl ) \
        FLASH_Queue_push( (void *) FLA_Apply_Q_UT_task, \
                          (void *) cntl, \
                          "ApQ  ", \
                          FALSE, \
                          FALSE, \
                          4, 1, 1, 2, \
                          side, trans, direct, storev, \
                          T, A, B, W )

#define ENQUEUE_FLASH_Apply_Q2_UT( side, trans, direct, storev, D, T, W, C, E, cntl ) \
        FLASH_Queue_push( (void *) FLA_Apply_Q2_UT_task, \
                          (void *) cntl, \
                          "ApQ2 ", \
                          FALSE, \
                          FALSE, \
                          4, 1, 1, 3, \
                          side, trans, direct, storev, \
                          T, D, E, C, W )

#define ENQUEUE_FLASH_Apply_CAQ2_UT( side, trans, direct, storev, D, T, W, C, E, cntl ) \
        FLASH_Queue_push( (void *) FLA_Apply_CAQ2_UT_task, \
                          (void *) cntl, \
                          "ApCQ2", \
                          FALSE, \
                          FALSE, \
                          4, 1, 1, 3, \
                          side, trans, direct, storev, \
                          T, D, E, C, W )

#define ENQUEUE_FLASH_UDdate_UT( R, C, D, T, cntl ) \
        FLASH_Queue_push( (void *) FLA_UDdate_UT_task, \
                          (void *) cntl, \
                          "UD   ", \
                          FALSE, \
                          FALSE, \
                          0, 0, 0, 4, \
                          R, C, D, T )

#define ENQUEUE_FLASH_Apply_QUD_UT( side, trans, direct, storev, T, W, R, U, C, V, D, cntl ) \
        FLASH_Queue_push( (void *) FLA_Apply_QUD_UT_task, \
                          (void *) cntl, \
                          "ApQUD", \
                          FALSE, \
                          FALSE, \
                          4, 0, 3, 4, \
                          side, trans, direct, storev, \
                          T, U, V, W, R, C, D )

#define ENQUEUE_FLASH_Eig_gest( inv, uplo, A, Y, B, cntl ) \
        FLASH_Queue_push( (void *) FLA_Eig_gest_task, \
                          (void *) cntl, \
                          "Eig  ", \
                          FALSE, \
                          TRUE, \
                          2, 0, 1, 2, \
                          inv, uplo, \
                          B, Y, A )

// Level-3 BLAS

#define ENQUEUE_FLASH_Gemm( transA, transB, alpha, A, B, beta, C, cntl ) \
        FLASH_Queue_push( (void *) FLA_Gemm_task, \
                          (void *) cntl, \
                          "Gemm ", \
                          TRUE, \
                          TRUE, \
                          2, 2, 2, 1, \
                          transA, transB, \
                          alpha, beta, \
                          A, B, C )

#define ENQUEUE_FLASH_Hemm( side, uplo, alpha, A, B, beta, C, cntl ) \
        FLASH_Queue_push( (void *) FLA_Hemm_task, \
                          (void *) cntl, \
                          "Hemm ", \
                          TRUE, \
                          TRUE, \
                          2, 2, 2, 1, \
                          side, uplo, \
                          alpha, beta, \
                          A, B, C )

#define ENQUEUE_FLASH_Herk( uplo, transA, alpha, A, beta, C, cntl ) \
        FLASH_Queue_push( (void *) FLA_Herk_task, \
                          (void *) cntl, \
                          "Herk ", \
                          TRUE, \
                          TRUE, \
                          2, 2, 1, 1, \
                          uplo, transA, \
                          alpha, beta, \
                          A, C )

#define ENQUEUE_FLASH_Her2k( uplo, transA, alpha, A, B, beta, C, cntl ) \
        FLASH_Queue_push( (void *) FLA_Her2k_task, \
                          (void *) cntl, \
                          "Her2k", \
                          TRUE, \
                          TRUE, \
                          2, 2, 2, 1, \
                          uplo, transA, \
                          alpha, beta, \
                          A, B, C )

#define ENQUEUE_FLASH_Symm( side, uplo, alpha, A, B, beta, C, cntl ) \
        FLASH_Queue_push( (void *) FLA_Symm_task, \
                          (void *) cntl, \
                          "Symm ", \
                          TRUE, \
                          TRUE, \
                          2, 2, 2, 1, \
                          side, uplo, \
                          alpha, beta, \
                          A, B, C )

#define ENQUEUE_FLASH_Syrk( uplo, transA, alpha, A, beta, C, cntl ) \
        FLASH_Queue_push( (void *) FLA_Syrk_task, \
                          (void *) cntl, \
                          "Syrk ", \
                          TRUE, \
                          TRUE, \
                          2, 2, 1, 1, \
                          uplo, transA, \
                          alpha, beta, \
                          A, C )

#define ENQUEUE_FLASH_Syr2k( uplo, transA, alpha, A, B, beta, C, cntl ) \
        FLASH_Queue_push( (void *) FLA_Syr2k_task, \
                          (void *) cntl, \
                          "Syr2k", \
                          TRUE, \
                          TRUE, \
                          2, 2, 2, 1, \
                          uplo, transA, \
                          alpha, beta, \
                          A, B, C )

#define ENQUEUE_FLASH_Trmm( side, uplo, trans, diag, alpha, A, C, cntl ) \
        FLASH_Queue_push( (void *) FLA_Trmm_task, \
                          (void *) cntl, \
                          "Trmm ", \
                          TRUE, \
                          TRUE, \
                          4, 1, 1, 1, \
                          side, uplo, trans, diag, \
                          alpha, \
                          A, C )

#define ENQUEUE_FLASH_Trsm( side, uplo, trans, diag, alpha, A, C, cntl ) \
        FLASH_Queue_push( (void *) FLA_Trsm_task, \
                          (void *) cntl, \
                          "Trsm ", \
                          TRUE, \
                          TRUE, \
                          4, 1, 1, 1, \
                          side, uplo, trans, diag, \
                          alpha, \
                          A, C )

// Level-2 BLAS

#define ENQUEUE_FLASH_Gemv( trans, alpha, A, x, beta, y, cntl ) \
        FLASH_Queue_push( (void *) FLA_Gemv_task, \
                          (void *) cntl, \
                          "Gemv ", \
                          TRUE, \
                          TRUE, \
                          1, 2, 2, 1, \
                          trans, \
                          alpha, beta, \
                          A, x, y )

#define ENQUEUE_FLASH_Trsv( uplo, trans, diag, A, x, cntl ) \
        FLASH_Queue_push( (void *) FLA_Trsv_task, \
                          (void *) cntl, \
                          "Trsv ", \
                          TRUE, \
                          TRUE, \
                          3, 0, 1, 1, \
                          uplo, trans, diag, \
                          A, x )

// Level-1 BLAS

#define ENQUEUE_FLASH_Axpy( alpha, A, B, cntl ) \
        FLASH_Queue_push( (void *) FLA_Axpy_task, \
                          (void *) cntl, \
                          "Axpy ", \
                          TRUE, \
                          TRUE, \
                          0, 1, 1, 1, \
                          alpha, \
                          A, B )

#define ENQUEUE_FLASH_Axpyt( trans, alpha, A, B, cntl ) \
        FLASH_Queue_push( (void *) FLA_Axpyt_task, \
                          (void *) cntl, \
                          "Axpyt", \
                          FALSE, \
                          FALSE, \
                          1, 1, 1, 1, \
                          trans, \
                          alpha, \
                          A, B )

#define ENQUEUE_FLASH_Copy( A, B, cntl ) \
        FLASH_Queue_push( (void *) FLA_Copy_task, \
                          (void *) cntl, \
                          "Copy ", \
                          TRUE, \
                          TRUE, \
                          0, 0, 1, 1, \
                          A, B )

#define ENQUEUE_FLASH_Copyt( trans, A, B, cntl ) \
        FLASH_Queue_push( (void *) FLA_Copyt_task, \
                          (void *) cntl, \
                          "Copyt", \
                          FALSE, \
                          FALSE, \
                          1, 0, 1, 1, \
                          trans, \
                          A, B )

#define ENQUEUE_FLASH_Copyr( uplo, A, B, cntl ) \
        FLASH_Queue_push( (void *) FLA_Copyr_task, \
                          (void *) cntl, \
                          "Copyr", \
                          FALSE, \
                          TRUE, \
                          1, 0, 1, 1, \
                          uplo, \
                          A, B )

#define ENQUEUE_FLASH_Scal( alpha, A, cntl ) \
        FLASH_Queue_push( (void *) FLA_Scal_task, \
                          (void *) cntl, \
                          "Scal ", \
                          TRUE, \
                          TRUE, \
                          0, 1, 0, 1, \
                          alpha, \
                          A )

#define ENQUEUE_FLASH_Scalr( uplo, alpha, A, cntl ) \
        FLASH_Queue_push( (void *) FLA_Scalr_task, \
                          (void *) cntl, \
                          "Scalr", \
                          TRUE, \
                          TRUE, \
                          1, 1, 0, 1, \
                          uplo, \
                          alpha, \
                          A )

// Base

#define ENQUEUE_FLASH_Obj_create_buffer( rs, cs, A, cntl ) \
        FLASH_Queue_push( (void *) FLA_Obj_create_buffer_task, \
                          (void *) cntl, \
                          "Buff ", \
                          FALSE, \
                          FALSE, \
                          2, 0, 0, 1, \
                          rs, cs, \
                          A )

#define ENQUEUE_FLASH_Obj_free_buffer( A, cntl ) \
        FLASH_Queue_push( (void *) FLA_Obj_free_buffer_task, \
                          (void *) cntl, \
                          "Free ", \
                          FALSE, \
                          FALSE, \
                          0, 0, 0, 1, \
                          A )

#else

// LAPACK-level

#define ENQUEUE_FLASH_LU_piv_macro( A, p, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Apply_pivots_macro( side, trans, p, A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_LU_piv( A, p, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_LU_piv_copy( A, p, U, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Trsm_piv( A, C, p, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_SA_LU( U, D, p, L, nb_alg, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_SA_FS( L, D, p, C, E, nb_alg, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_LU_nopiv( A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Trinv( uplo, diag, A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Ttmm( uplo, A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Chol( uplo, A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Sylv( transA, transB, isgn, A, B, C, scale, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Lyap( trans, isgn, A, C, scale, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_QR_UT_macro( A, T, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_QR_UT( A, T, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_QR_UT_copy( A, T, U, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_QR2_UT( B, D, T, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_LQ_UT_macro( A, T, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_CAQR2_UT( B, D, T, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_UDdate_UT( R, C, D, T, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Apply_Q_UT( side, trans, direct, storev, A, T, W, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Apply_Q2_UT( side, trans, direct, storev, D, T, W, C, E, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Apply_CAQ2_UT( side, trans, direct, storev, D, T, W, C, E, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Apply_QUD_UT( side, trans, direct, storev, T, W, R, U, C, V, D, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Eig_gest( inv, uplo, A, Y, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

// Level-3 BLAS

#define ENQUEUE_FLASH_Gemm( transA, transB, alpha, A, B, beta, C, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Hemm( side, uplo, alpha, A, B, beta, C, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Herk( uplo, transA, alpha, A, beta, C, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Her2k( uplo, transA, alpha, A, B, beta, C, cntl  ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Symm( side, uplo, alpha, A, B, beta, C, cntl  ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Syrk( uplo, transA, alpha, A, beta, C, cntl  ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Syr2k( uplo, transA, alpha, A, B, beta, C, cntl  ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Trmm( side, uplo, trans, diag, alpha, A, C, cntl  ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Trsm( side, uplo, trans, diag, alpha, A, C, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

// Level-2 BLAS

#define ENQUEUE_FLASH_Gemv( transA, alpha, A, x, beta, y, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Trsv( uplo, trans, diag, A, x, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

// Level-1 BLAS

#define ENQUEUE_FLASH_Axpy( alpha, A, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Axpyt( trans, alpha, A, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Copy( A, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Copyt( trans, A, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Copyr( uplo, A, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Scal( alpha, A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Scalr( uplo, alpha, A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

// Base

#define ENQUEUE_FLASH_Obj_create_buffer( rs, cs, A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Obj_free_buffer( A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#endif // FLA_ENABLE_SUPERMATRIX


#endif // FLASH_QUEUE_MACRO_DEFS_H
