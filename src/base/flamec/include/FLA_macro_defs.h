			/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
    Modifications Copyright (C) 2021-2026, Advanced Micro Devices, Inc. All rights reserved.
*/



// --- Miscellaneous macro definitions -----------------------------------------
#ifndef FLA_MACRO_DEFS_H
#define FLA_MACRO_DEFS_H

#undef  NULL
#define NULL 0

#ifdef FLA_ENABLE_WINDOWS_BUILD
  #define restrict  __restrict
#endif

// --- Macro to enable/disable Thread Local Storage (TLS) for global variables -
#if defined(WINDOWS_FLA_SHARED_BUILD) && defined(WINDOWS_FLA_TEST)
 #define ENABLE_THREAD_LOCAL_STORAGE 0
 #define LIBFLAME_IMPORT __declspec(dllimport)
#else
 #if defined(FLA_ENABLE_AMD_OPT) && !defined FLA_ENABLE_SUPERMATRIX
   #define ENABLE_THREAD_LOCAL_STORAGE 1
 #else
   #define ENABLE_THREAD_LOCAL_STORAGE 0
 #endif
 #define LIBFLAME_IMPORT
#endif

#if ENABLE_THREAD_LOCAL_STORAGE
   #if _WIN32
      #if defined(FLA_ENABLE_WINDOWS_BUILD)
         #define TLS_CLASS_SPEC __declspec(thread)
      #else
         #define TLS_CLASS_SPEC
      #endif
   #else
     #define TLS_CLASS_SPEC __thread
   #endif
#else
   #define TLS_CLASS_SPEC
#endif

// --- Type-related macro definitions ------------------------------------------

// FLA_Bool
#undef  TRUE
#undef  FALSE
#define TRUE  1
#define FALSE 0

// FLA_Error (non-specific)
#define FLA_SUCCESS           (-1)
#define FLA_FAILURE           (-2)

// FLA_Quadrant
#define FLA_TL                 11
#define FLA_TR                 12
#define FLA_BL                 21
#define FLA_BR                 22

// FLA_Datatype
#define FLA_FLOAT             100
#define FLA_DOUBLE            101
#define FLA_COMPLEX           102
#define FLA_DOUBLE_COMPLEX    103
#define FLA_INT               104
#define FLA_CONSTANT          105

// FLA_Elemtype
#define FLA_MATRIX            150
#define FLA_SCALAR            151

// FLA_Side
#define FLA_TOP               200
#define FLA_BOTTOM            201
#define FLA_LEFT              210
#define FLA_RIGHT             211
#define FLA_SIDE_MASK         0x1

// FLA_Uplo
#define FLA_LOWER_TRIANGULAR  300
#define FLA_UPPER_TRIANGULAR  301
#define FLA_ZERO_MATRIX       310
#define FLA_FULL_MATRIX       311
#define FLA_UPLO_MASK         0x1

// FLA_Trans
#define FLA_NO_TRANSPOSE      400
#define FLA_TRANSPOSE         401
#define FLA_CONJ_TRANSPOSE    402
#define FLA_CONJ_NO_TRANSPOSE 403
#define FLA_TRANS_MASK        0x3

// FLA_Conj
#define FLA_NO_CONJUGATE      450
#define FLA_CONJUGATE         451

// FLA_Diag
#define FLA_UNIT_DIAG         500
#define FLA_NONUNIT_DIAG      501
#define FLA_ZERO_DIAG         502
#define FLA_DIAG_MASK         0x3

// FLA_Dimension
#define FLA_DIMENSION_M       600
#define FLA_DIMENSION_K       601
#define FLA_DIMENSION_N       602
#define FLA_DIMENSION_MIN     603

// FLA_Dimension_index
#define FLA_DIM_M_INDEX         0
#define FLA_DIM_K_INDEX         1
#define FLA_DIM_N_INDEX         2
#define FLA_DIM_MIN_INDEX       3
#define FLA_DIM_INDEX_MASK    0x3

// FLA_Pivot_type
#define FLA_NATIVE_PIVOTS     700
#define FLA_LAPACK_PIVOTS     701

// FLA_Direct
#define FLA_FORWARD           800
#define FLA_BACKWARD          801

// FLA_Store
#define FLA_COLUMNWISE        900
#define FLA_ROWWISE           901

// FLA_Matrix_type
#define FLA_FLAT             1000
#define FLA_HIER             1001

// FLA_Precision
#define FLA_SINGLE_PRECISION 1100
#define FLA_DOUBLE_PRECISION 1101

// FLA_Domain
#define FLA_REAL_DOMAIN      1200
#define FLA_COMPLEX_DOMAIN   1201

// FLA_Inv
#define FLA_NO_INVERSE       1300
#define FLA_INVERSE          1301

// FLA_Evd_type
#define FLA_EVD_WITHOUT_VECTORS         1400
#define FLA_EVD_WITH_VECTORS            1401
#define FLA_EVD_OF_TRIDIAG_WITH_VECTORS 1402

// FLA_Svd_type
#define FLA_SVD_VECTORS_ALL           1500
#define FLA_SVD_VECTORS_MIN_COPY      1501
#define FLA_SVD_VECTORS_MIN_OVERWRITE 1502
#define FLA_SVD_VECTORS_NONE          1503

// FLA_Machval
#define FLA_MACH_START                1600
#define FLA_MACH_EPS                  1600
#define FLA_MACH_SFMIN                1601
#define FLA_MACH_BASE                 1602
#define FLA_MACH_PREC                 1603
#define FLA_MACH_NDIGMANT             1604
#define FLA_MACH_RND                  1605
#define FLA_MACH_EMIN                 1606
#define FLA_MACH_RMIN                 1607
#define FLA_MACH_EMAX                 1608
#define FLA_MACH_RMAX                 1609
#define FLA_MACH_EPS2                 1610
#define FLA_MACH_N_VALS                 11

// FLA_Diag_off
#define FLA_SUPER_DIAGONAL     ( 1)
#define FLA_MAIN_DIAGONAL        0
#define FLA_SUB_DIAGONAL       (-1)

// FLAME threading model
#define FLA_OPENMP              1
#define FLA_PTHREADS            2

// FLAME vector intrinsics types
#define FLA_NO_INTRINSICS       0
#define FLA_SSE_INTRINSICS      3

// FLAME internal error checking level
#define FLA_FULL_ERROR_CHECKING 2
#define FLA_MIN_ERROR_CHECKING  1
#define FLA_NO_ERROR_CHECKING   0

// FLA_Datatype_index
#define FLA_S_INDEX             0
#define FLA_D_INDEX             1
#define FLA_C_INDEX             2
#define FLA_Z_INDEX             3
#define FLA_DTYPE_INDEX_MASK  0x3

#define FLA_SMALL_LU_BLOCKSIZE 64
// Default blocksize if none are available.
#ifndef FLA_DEFAULT_M_BLOCKSIZE
  #define FLA_DEFAULT_M_BLOCKSIZE  128
#endif
#ifndef FLA_DEFAULT_K_BLOCKSIZE
  #define FLA_DEFAULT_K_BLOCKSIZE  128
#endif
#ifndef FLA_DEFAULT_N_BLOCKSIZE
  #define FLA_DEFAULT_N_BLOCKSIZE  128
#endif

// Blocksize for hierarchical storage of matrices for FLASH based
// external interfaces
#define FLA_EXT_HIER_BLOCKSIZE 128

// QR and LQ factorizations typically has an inner blocksize that corresponds
// to the length of the S (or T) block Householder matrix. For consistency, we
// define the ratio of the inner blocksize to the outer blocksize here, as it
// is used in several places. Note that other operations have analagous inner
// blocksizes, which we also define in terms of the outer storage blocksize,
// or in some cases such as Hessenberg, tridiagonal, and bidiagonal reductions,
// in terms of the system-wide default blocksize.
#define FLA_QR_INNER_TO_OUTER_B_RATIO      (0.25)
#define FLA_LQ_INNER_TO_OUTER_B_RATIO      (0.25)
#define FLA_LU_INNER_TO_OUTER_B_RATIO      (0.25)
#define FLA_UDDATE_INNER_TO_OUTER_B_RATIO  (0.25)
#define FLA_HESS_INNER_TO_OUTER_B_RATIO    (0.25)
#define FLA_TRIDIAG_INNER_TO_OUTER_B_RATIO (0.25)
#define FLA_BIDIAG_INNER_TO_OUTER_B_RATIO  (0.25)
#define FLA_CAQR_INNER_TO_OUTER_B_RATIO    (0.25)

// Matrix size thresholds for choosing the unblocked small-QR kernel
#define FLA_GEQRF_STHRESH (32)
#define FLA_CGEQRF_STHRESH (48)
#define FLA_ZGEQRF_STHRESH (22)
#define FLA_GEQRF_BLOCK_SIZE (32)

// Cache line size in bytes, used for memory alignment
#define FLA_CACHE_LINE_SIZE_BYTES (64)

// Matrix size thresholds for choosing multithreaded QR for large matrices
#define FLA_GEQRF_MT_LARGE_M_THRESH (128)
#define FLA_GEQRF_MT_LARGE_N_THRESH (128)
#define FLA_GEQRF_MT_LARGE_PANEL_SIZE (16)

#define FLA_GEQRF_MT_THRESHOLD_8_THREADS_REAL (15848840)
#define FLA_CGEQRF_MT_THRESHOLD_8_THREADS (1500000)
#define FLA_ZGEQRF_MT_THRESHOLD_8_THREADS (4000000)

#define FLA_SGEQRF_MT_MIN_ELEMS (50000)
#define FLA_GEQRF_MT_HUGE_MIN_THREADS (32)
#define FLA_GEQRF_MT_WIDE_ASPECT (2)
#define FLA_GEQRF_MT_WIDE_MIN_ELEMS (1000000)
#define FLA_CGEQRF_MT_HUGE_THRESH (9000000)
#define FLA_ZGEQRF_MT_HUGE_THRESH (25000000)

// Matrix size thresholds for choosing transposed QR for LQ
#define FLA_DELQF_TRAN_THRESH (1280)

//POTRF, threshold numbers to chose paths for performance
#define FLA_POTRF_FLOAT_SMALL         (70)
#define FLA_POTRF_DOUBLE_SMALL        (75)
#define FLA_DPOTRF_LOWER_SMALL        (18)
#define FLA_DPOTRF_LOWER_MEDIUM       (75)
#define FLA_DPOTRF_UPPER_MEDIUM       (50)
#define FLA_POTRF_BLOCK_SIZE          (256)

// POTRI threshold numbers to chose paths for performance
#define FLA_POTRI_DOUBLE_SMALL        (65)

//GTSV, threshold numbers to chose paths for performance.
// The optimized kernels vectorize across the NRHS columns (8-wide for AVX512,
// 4-wide for AVX2), so the SIMD path is only chosen when there are enough
// right-hand sides to fill the vector registers.
#define FLA_DGTSV_AVX512_MIN_NRHS     (8)
#define FLA_DGTSV_AVX2_MIN_NRHS       (4)

//GETRF , threshold numbers to chose paths for performance
#define TALL_RATIO_THRESHOLD 1.5
#define WIDE_RATIO_THRESHOLD 0.67
#define FLA_GETRF_SMALL               (85)
#define FLA_GETRF_FLOAT               (1440)
#define FLA_GETRF_COMPLEX             (670)
#define FLA_GETRF_DOUBLE_COMPLEX      (1630)

#define FLA_SGETRF_SMALL_THRESH1      (117)
#define FLA_SGETRF_MEDIUM_THRESH0     (480)

#define FLA_SGETRF_SMALL_THRESH0      (7)
#define FLA_SGETRF_SMALL_AVX2_THRESH0      (31)
#define FLA_SGETRF_SMALL_AVX512_THRESH0    (156)

#define FLA_DGETRF_SMALL_THRESH0      (8)
#define FLA_DGETRF_SMALL_AVX2_THRESH0      (81)
#define FLA_DGETRF_SMALL_AVX512_THRESH0    (160)
#define FLA_DGETRF_SMALL_AVX512_THRESH1    (128)

#define FLA_ZGETRF_SMALL_THRESH      (100)

#define FLA_CGETRF_SMALL_THRESH0 (10)
#define FLA_CGETRF_SMALL_THRESH (100)

#define FLA_SSYTD2_SMALL_THRESH      (256)
// GETRFNPI , these thresholds are used to chose between 3 algorithms to get best
// results in terms of perfmormance
#define FLA_MN_SIZE            (373321)
#define FLA_NFACT_PERCENT      (0.51)
#define FLA_FULL_DGER_CONSTANT (200)

// ORGQR , threshold numbers to chose paths for performance
#define FLA_DGEMV_DGER_SIMD_SMALL_THRESH   (128)
#define FLA_DGEMV_DGER_SIMD_SMALL_THRESH_M (6)
#define FLA_DGEMV_DGER_SIMD_AVX512_THRESH_M (25)
#define FLA_ZGEMV_ZGER_SIMD_AXV2_THRESH_M (8)
#define FLA_ZGEMV_ZGER_SIMD_AXV512_THRESH_M (49)
#define FLA_SSCAL_INLINE_SMALL          (512)
#define FLA_DSCAL_INLINE_SMALL          (512)
#define FLA_ZGERC_INLINE_SMALL_THRESH0  (128)
#define FLA_ZGERC_INLINE_SMALL_THRESH1  (32)
#define FLA_ZSCAL_INLINE_SMALL          (128)
#define FLA_ZUNGQR_BLOCK_SMALL_THRESH   (32)

/* DLARF , threshold numbers to choose paths for performance */
#define FLA_DLARF_L_ST_BLOCKED_THRESH_M (1000)
#define FLA_DLARF_L_ST_BLOCKED_THRESH_N (100)
/* DLARF1F/DLARF1L applied from the left , the blocked path is already fused per
   column, so against it the SIMD kernel only saves the DDOT/DAXPY call overhead.
   That is worth it while C stays inside L3, past which both are memory bound and
   the tuned BLAS kernels stream better */
#define FLA_DLARF1_L_FUSED_SIMD_THRESH (1048576)
/* DLARF1F/DLARF1L applied from the right , the fused kernel makes two passes over
   each row block of C and touches one cache line per column, so it needs C to be
   narrow enough for a row block to fit the L1 way count and small enough overall
   to keep the reuse between the two passes. Both have to hold: C that is narrow
   but tall loses the reuse just as C that is small but wide does */
#define FLA_DLARF1_R_SIMD_COLS_THRESH (96)
#define FLA_DLARF1_R_SIMD_THRESH (4096)
/* DLARF1F/DLARF1L applied from the right , the fused kernel takes a block of
   rows of C and makes two passes over all *n of its columns, relying on the
   block staying L1 resident in between. Successive columns are *ldc doubles
   apart, and x86 L1 data caches have 64 sets of 64 byte lines, so addresses
   4096 bytes apart collide. When *ldc is a multiple of 64 doubles a row block
   reaches at most 8 of the 64 sets, the reuse the fusion depends on is lost
   and the kernel runs up to 20x slower. GEMV + GER stream over C in two
   independent passes and need no such reuse, so hand these cases to them */
#define FLA_DLARF1_R_SIMD_LDC_ALIAS_MASK (63)
/* DLARF1F/DLARF1L applied from the right , minimum number of rows for which the
   AVX512 kernel is faster than the AVX2 one. Below this the AVX512 kernel only
   executes its 256/128 bit tails, which are identical to the AVX2 kernel */
#define FLA_DLARF1_R_SIMD_AVX512_THRESH_M (8)
#define FLA_DLARF_L_THRESH_UNBLOCKED (213760)
#define FLA_DLARF_L_THRESH_THREAD_8 (1048576)
#define FLA_DLARF_L_THRESH_THREAD_64 (45198729)

#define FLA_DLARF_R_THRESH_UNBLOCKED (15000)
#define FLA_DLARF_R_THRESH_THREAD_8 (6553600)
#define FLA_DLARF_R_THRESH_THREAD_64 (45198729)
#define FLA_DLARF_R_BLOCK_SIZE (16)

/* DLARFT , threshold numbers to chose paths for performance */
#define FLA_DLARFT_BLOCK_NB (16)

#define FLA_DLASR_L_SIMD_AVX512_THRESH_N (8)

/* ZLARF1F threshold numbers */
#ifndef FLA_ZLARF1_SMALL_DIM_THRESH
#define FLA_ZLARF1_SMALL_DIM_THRESH (34)
#endif
#ifndef FLA_ZLARF1_R_SIMD_AVX512_THRESH_M
#define FLA_ZLARF1_R_SIMD_AVX512_THRESH_M (8)
#endif

#define FLA_ZLANGEM_SIMD_AVX512_THRESH_M (15)
#define FLA_CLANGEM_SIMD_AVX512_THRESH_M (15)

// DNRM2 threshold numbers to chose AVX kernel
#define FLA_DNRM2_SMALL_THRESH0  (15)

// GEQP3 , threshold numbers to chose paths for performance
#define FLA_IDAMAX_INLINE_SMALL_THRESH   (128)
#define FLA_DGEQP3_BLOCK_SMALL_THRESH    (32)
#define FLA_ZGEQP3_BLOCK_SMALL_THRESH    (32)

// GESVD, DGESDD thresholds for small size optimization
#define FLA_GESVD_SMALL_SIZE_THRESH0    (16)
#define FLA_GESVD_SMALL_SIZE_THRESH1    (128)
#define FLA_GESVD_SMALL_SIZE_THRESH2    (386)
#define FLA_DGESDD_SMALL_SIZE_THRESH     (40)

// TRTRI, threshold numbers to chose paths for performance
#define FLA_TRTRI_SMALL_THRESH0      (60)

// DGEBRD thresholds for lapack_API based code path for large sizes
#define FLA_DGEBRD_LARGE_SIZE_THRESH    (10000)
#define FLA_DGEBRD_SMALL_SIZE_THRESH    (160)

// DGEBD2 thresholds for lapack_API based code path for small sizes
#define FLA_DGEBD2_SMALL_SIZE_THRESH   (140)

// DORMQR/DORMLQ threshold for small vs large sizes
#define FLA_DORMQR_SMALL_SIZE_THRESH (100)
#define FLA_DORMLQ_SMALL_SIZE_THRESH (50)

// DORMLQ, work (m*n*k) thresholds selecting 8 / 16 / 32 / 64 threads.
// Below each one the work per thread is too small to cover its share of the
// threading overhead.
#define FLA_DORMLQ_THREAD_THRESH0 (1200000000LL)
#define FLA_DORMLQ_THREAD_THRESH1 (6000000000LL)
#define FLA_DORMLQ_THREAD_THRESH2 (50000000000LL)

// DORMLQ, smallest share of the split dimension (n for side='L', m for side='R')
// worth giving a thread, for the same reason.
#define FLA_DORMLQ_MIN_SPLIT_PER_THREAD (12)

// --- Error-related macro definitions -----------------------------------------

// Useful when determining the relative index base of the error codes.
#define FLA_ERROR_CODE_MIN                    (-10)

// FLA_Error values.
#define FLA_INVALID_SIDE                      (-10)
#define FLA_INVALID_UPLO                      (-11)
#define FLA_INVALID_TRANS                     (-12)
#define FLA_INVALID_TRANS_GIVEN_DATATYPE      (-13)
#define FLA_INVALID_CONJ                      (-14)
#define FLA_INVALID_DIRECT                    (-15)
#define FLA_INVALID_STOREV                    (-16)
#define FLA_INVALID_DATATYPE                  (-17)
#define FLA_INVALID_INTEGER_DATATYPE          (-18)
#define FLA_INVALID_REAL_DATATYPE             (-19)
#define FLA_INVALID_COMPLEX_DATATYPE          (-20)
#define FLA_OBJECT_NOT_INTEGER                (-21)
#define FLA_OBJECT_NOT_REAL                   (-22)
#define FLA_OBJECT_NOT_COMPLEX                (-23)
#define FLA_OBJECT_NOT_SQUARE                 (-24)
#define FLA_OBJECT_NOT_SCALAR                 (-25)
#define FLA_OBJECT_NOT_VECTOR                 (-26)
#define FLA_INCONSISTENT_DATATYPES            (-27)
#define FLA_NONCONFORMAL_DIMENSIONS           (-28)
#define FLA_UNEQUAL_VECTOR_DIMS               (-29)
#define FLA_INVALID_HESSENBERG_INDICES        (-30)
#define FLA_NULL_POINTER                      (-32)
#define FLA_SPECIFIED_OBJ_DIM_MISMATCH        (-33)
#define FLA_INVALID_PIVOT_TYPE                (-35)
#define FLA_MALLOC_RETURNED_NULL_POINTER      (-37)
#define FLA_OBJECT_BASE_BUFFER_MISMATCH       (-38)
#define FLA_OBJECTS_NOT_VERTICALLY_ADJ        (-39)
#define FLA_OBJECTS_NOT_HORIZONTALLY_ADJ      (-40)
#define FLA_ADJACENT_OBJECT_DIM_MISMATCH      (-41)
#define FLA_OBJECTS_NOT_VERTICALLY_ALIGNED    (-42)
#define FLA_OBJECTS_NOT_HORIZONTALLY_ALIGNED  (-43)
#define FLA_INVALID_FLOATING_DATATYPE         (-44)
#define FLA_OBJECT_NOT_FLOATING_POINT         (-45)
#define FLA_INVALID_BLOCKSIZE_VALUE           (-46)
#define FLA_OPEN_RETURNED_ERROR               (-47)
#define FLA_LSEEK_RETURNED_ERROR              (-48)
#define FLA_CLOSE_RETURNED_ERROR              (-49)
#define FLA_UNLINK_RETURNED_ERROR             (-50)
#define FLA_READ_RETURNED_ERROR               (-51)
#define FLA_WRITE_RETURNED_ERROR              (-52)
#define FLA_INVALID_QUADRANT                  (-53)
#define FLA_NOT_YET_IMPLEMENTED               (-54)
#define FLA_EXPECTED_NONNEGATIVE_VALUE        (-55)
#define FLA_SUPERMATRIX_NOT_ENABLED           (-56)
#define FLA_UNDEFINED_ERROR_CODE              (-57)
#define FLA_INVALID_DIAG                      (-58)
#define FLA_INCONSISTENT_OBJECT_PRECISION     (-59)
#define FLA_INVALID_BLOCKSIZE_OBJ             (-60)
#define FLA_VECTOR_DIM_BELOW_MIN              (-61)
#define FLA_PTHREAD_CREATE_RETURNED_ERROR     (-63)
#define FLA_PTHREAD_JOIN_RETURNED_ERROR       (-64)
#define FLA_INVALID_ISGN_VALUE                (-65)
#define FLA_CHOL_FAILED_MATRIX_NOT_SPD        (-67)
#define FLA_INVALID_ELEMTYPE                  (-68)
#define FLA_POSIX_MEMALIGN_FAILED             (-69)
#define FLA_INVALID_SUBMATRIX_DIMS            (-70)
#define FLA_INVALID_SUBMATRIX_OFFSET          (-71)
#define FLA_OBJECT_NOT_SCALAR_ELEMTYPE        (-72)
#define FLA_OBJECT_NOT_MATRIX_ELEMTYPE        (-73)
#define FLA_ENCOUNTERED_NON_POSITIVE_NTHREADS (-74)
#define FLA_INVALID_CONJ_GIVEN_DATATYPE       (-75)
#define FLA_INVALID_COMPLEX_TRANS             (-76)
#define FLA_INVALID_REAL_TRANS                (-77)
#define FLA_INVALID_BLAS_TRANS                (-78)
#define FLA_INVALID_NONCONSTANT_DATATYPE      (-79)
#define FLA_OBJECT_NOT_NONCONSTANT            (-80)
#define FLA_OBJECT_DATATYPES_NOT_EQUAL        (-82)
#define FLA_DIVIDE_BY_ZERO                    (-83)
#define FLA_OBJECT_ELEMTYPES_NOT_EQUAL        (-84)
#define FLA_INVALID_PIVOT_INDEX_RANGE         (-85)
#define FLA_HOUSEH_PANEL_MATRIX_TOO_SMALL     (-86)
#define FLA_INVALID_OBJECT_LENGTH             (-87)
#define FLA_INVALID_OBJECT_WIDTH              (-88)
#define FLA_INVALID_ERROR_CHECKING_LEVEL      (-89)
#define FLA_ATTEMPTED_OVER_REPART_2X2         (-90)
#define FLA_ATTEMPTED_OVER_REPART_2X1         (-91)
#define FLA_ATTEMPTED_OVER_REPART_1X2         (-92)
#define FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED   (-93)
#define FLA_INVALID_ROW_STRIDE                (-94)
#define FLA_INVALID_COL_STRIDE                (-95)
#define FLA_INVALID_STRIDE_COMBINATION        (-96)
#define FLA_INVALID_VECTOR_DIM                (-97)
#define FLA_EXPECTED_ROW_VECTOR               (-98)
#define FLA_EXPECTED_COL_VECTOR               (-99)
#define FLA_INVALID_INVERSE                   (-100)
#define FLA_MALLOC_GPU_RETURNED_NULL_POINTER  (-101)
#define FLA_INVALID_EVD_TYPE                  (-102)
#define FLA_INVALID_SVD_TYPE                  (-103)
#define FLA_INVALID_MACHVAL                   (-104)
#define FLA_INVALID_DIAG_OFFSET               (-105)
#define FLA_EXPECTED_COL_STORAGE              (-106)
#define FLA_EXPECTED_ROW_STORAGE              (-107)
#define FLA_LAPAC2FLAME_INVALID_RETURN        (-108)
#define FLA_INVALID_SVD_TYPE_COMBINATION      (-109)
#define FLA_INVALID_SVD_TYPE_AND_TRANS_COMBINATION (-110)
#define FLA_OBJECT_NOT_COMPARABLE             (-111)
#define FLA_MEMORY_ALLOCATION_FAILURE         (-112)

// Necessary when computing whether an error code is defined.
#define FLA_ERROR_CODE_MAX                    (-112)

// Internal string matrix limits.
#define FLA_MAX_NUM_ERROR_MSGS                 150
#define FLA_MAX_ERROR_MSG_LENGTH               200

// Error code translation and output macro definition.
#define FLA_Check_error_code( code ) \
        FLA_Check_error_code_helper( code, __FILE__, __LINE__ )



// --- Common functions implemented as macros ----------------------------------

#define fla_min( x, y )    \
({                         \
   __typeof__(x) _x = (x); \
   __typeof__(y) _y = (y); \
   _x < _y ? _x : _y;      \
})

#define fla_max( x, y )    \
({                         \
   __typeof__(x) _x = (x); \
   __typeof__(y) _y = (y); \
   _x > _y ? _x : _y;      \
})

#undef signof
#define signof( a, b ) ( (b) >= 0 ? (a) : -(a) )

#ifndef __cplusplus

 #undef exchange
 #define exchange( a, b, temp ) { temp = a; a = b; b = temp; }

#endif

// --- Other macro definitions -------------------------------------------------

#define FLA_NEGATE( a ) \
        ( a.base == FLA_ONE.base ? FLA_MINUS_ONE : FLA_ONE )

#endif

