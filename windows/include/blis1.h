

#ifndef BLIS1_H
#define BLIS1_H

// Allow C++ users to include this header file in their source code. However,
// we make the extern "C" conditional on whether we're using a C++ compiler,
// since regular C compilers don't understand the extern "C" construct.
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h> // skipped
#include <stdlib.h> // skipped
#include <math.h> // skipped

// Determine whether or not we are using BLIS from libflame.
//#define BLIS1_FROM_LIBFLAME

#ifdef BLIS1_FROM_LIBFLAME

  // If using libflame, pull in its header files so that
  // vector intrinsics-related macro constants are set properly.
  //#include "FLAME.h"
// begin FLA_config.h







#define F77_FUNC(name,NAME) name ## _


#define F77_FUNC_(name,NAME) name ## _

















#define FLA_ENABLE_BLIS1_USE_OF_FLA_MALLOC 1











#define FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES 1








#define FLA_ENABLE_INTERNAL_ERROR_CHECKING 1


#define FLA_ENABLE_LAPACK2FLAME 1














#define FLA_ENABLE_NON_CRITICAL_CODE 1


#define FLA_ENABLE_PORTABLE_TIMER 1

















#define FLA_INTERNAL_ERROR_CHECKING_LEVEL 2





#define FLA_MULTITHREADING_MODEL 0


#define FLA_PORTABLE_TIMER_IS_CLOCK_GETTIME 1








#define FLA_VECTOR_INTRINSIC_TYPE 0


#define HAVE_ASSERT_H 1


#define HAVE_FCNTL_H 1





#define HAVE_INTTYPES_H 1


#define HAVE_LIBM 1


#define HAVE_MATH_H 1


#define HAVE_MEMORY_H 1


#define HAVE_SIGNAL_H 1


#define HAVE_STDINT_H 1


#define HAVE_STDLIB_H 1


#define HAVE_STRINGS_H 1


#define HAVE_STRING_H 1


#define HAVE_SYS_STAT_H 1


#define HAVE_SYS_TIME_H 1


#define HAVE_SYS_TYPES_H 1


#define HAVE_UNISTD_H 1





#define PACKAGE_BUGREPORT ""


#define PACKAGE_NAME ""


#define PACKAGE_STRING ""


#define PACKAGE_TARNAME ""


#define PACKAGE_URL ""


#define PACKAGE_VERSION ""


#define PROTOTYPES 1


#define STDC_HEADERS 1


#define TIME_WITH_SYS_TIME 1


#define _GNU_SOURCE 1


#define __PROTOTYPES 1





#ifndef __cplusplus

#endif



// end FLA_config.h
// begin FLA_macro_defs.h




// --- Miscellaneous macro definitions -----------------------------------------

#undef  NULL
#define NULL 0

#ifdef FLA_ENABLE_WINDOWS_BUILD
  #define restrict  __restrict
#endif

// --- Macro to enable/disable Thread Local Storage (TLS) for global variables -
#define ENABLE_THREAD_LOCAL_STORAGE 1

#if ENABLE_THREAD_LOCAL_STORAGE
#define TLS_CLASS_SPEC __thread
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

// Necessary when computing whether an error code is defined.
#define FLA_ERROR_CODE_MAX                    (-111)

// Internal string matrix limits.
#define FLA_MAX_NUM_ERROR_MSGS                 150
#define FLA_MAX_ERROR_MSG_LENGTH               200

// Error code translation and output macro definition.
#define FLA_Check_error_code( code ) \
        FLA_Check_error_code_helper( code, __FILE__, __LINE__ )



// --- Common functions implemented as macros ----------------------------------

#undef min
#define min( x, y ) ( (x) < (y) ? (x) : (y) )

#undef max
#define max( x, y ) ( (x) > (y) ? (x) : (y) )

#undef signof
#define signof( a, b ) ( (b) >= 0 ? (a) : -(a) )

#undef exchange
#define exchange( a, b, temp ) { temp = a; a = b; b = temp; }

// --- Other macro definitions -------------------------------------------------

#define FLA_NEGATE( a ) \
        ( a.base == FLA_ONE.base ? FLA_MINUS_ONE : FLA_ONE )


// end FLA_macro_defs.h
// begin FLA_type_defs.h


#ifndef FLA_TYPE_DEFS_H
#define FLA_TYPE_DEFS_H

#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
#ifdef FLA_ENABLE_TIDSP
#include <ti/omp/omp.h> // skipped
#else
#include <omp.h> // skipped
#endif
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
#include <pthread.h> // skipped
#endif


// --- Complex type definitions -----------------------------------------------

#ifndef _DEFINED_SCOMPLEX
#define _DEFINED_SCOMPLEX
typedef struct scomplex
{
  float real, imag;
} scomplex;
#endif

#ifndef _DEFINED_DCOMPLEX
#define _DEFINED_DCOMPLEX
typedef struct dcomplex
{
  double real, imag;
} dcomplex;
#endif


// --- Parameter and return type definitions ----------------------------------

typedef int FLA_Bool;
typedef int FLA_Error;
typedef int FLA_Quadrant;
typedef int FLA_Datatype;
typedef int FLA_Elemtype;
typedef int FLA_Side;
typedef int FLA_Uplo;
typedef int FLA_Trans;
typedef int FLA_Conj;
typedef int FLA_Diag;
typedef int FLA_Dimension;
typedef int FLA_Pivot_type;
typedef int FLA_Direct;
typedef int FLA_Store;
typedef int FLA_Matrix_type;
typedef int FLA_Precision;
typedef int FLA_Domain;
typedef int FLA_Inv;
typedef int FLA_Evd_type;
typedef int FLA_Svd_type;
typedef int FLA_Machval;
typedef int FLA_Diag_off;

#ifndef _DEFINED_DIM_T
#define _DEFINED_DIM_T
typedef unsigned long dim_t;
#endif

// --- Intrinsic/assembly definitions ----------------------------------------

#if FLA_VECTOR_INTRINSIC_TYPE == FLA_SSE_INTRINSICS

#include "pmmintrin.h" // skipped

//typedef double v2df __attribute__ ((vector_size (16)));

typedef union
{
    __m128  v; 
    float   f[4];
} v4sf_t;

typedef union
{
    __m128d v; 
    double  d[2];
} v2df_t;

#endif

// --- FLAME object definitions -----------------------------------------------

typedef struct FLA_Lock_s     FLA_Lock;

//#ifdef FLA_ENABLE_MULTITHREADING
struct FLA_Lock_s
{
  // Implementation-specific lock object
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_lock_t       lock;
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_mutex_t  lock;
#endif
};
//#endif

#ifdef FLA_ENABLE_SUPERMATRIX
typedef int                   FLASH_Verbose;
typedef int                   FLASH_Data_aff;

typedef struct FLASH_Queue_s  FLASH_Queue;
typedef struct FLASH_Task_s   FLASH_Task;
typedef struct FLASH_Dep_s    FLASH_Dep;
#endif
typedef struct FLASH_Thread_s FLASH_Thread;

typedef struct FLA_Obj_struct
{
  // Basic object description fields
  FLA_Datatype  datatype;
  FLA_Elemtype  elemtype;
  dim_t         m;
  dim_t         n;
  dim_t         rs;
  dim_t         cs;
  dim_t         m_inner;
  dim_t         n_inner;
  unsigned long id;
  dim_t         m_index;
  dim_t         n_index;

  dim_t         n_elem_alloc;
  void*         buffer;
  int           buffer_info;

  FLA_Uplo      uplo;

#ifdef FLA_ENABLE_SUPERMATRIX
  // Fields for supermatrix
  int           n_read_blocks;
  int           n_write_blocks;

  // All the tasks that previously read this block, anti-dependency
  int           n_read_tasks;
  FLASH_Dep*    read_task_head;
  FLASH_Dep*    read_task_tail;

  // Task that last overwrote this block, flow dependency
  FLASH_Task*   write_task;
#endif
} FLA_Base_obj;

typedef struct FLA_Obj_view
{
  // Basic object view description fields
  dim_t         offm;
  dim_t         offn;
  dim_t         m;
  dim_t         n;
  dim_t         m_inner;
  dim_t         n_inner;

  FLA_Base_obj* base;

} FLA_Obj;

#ifdef FLA_ENABLE_SUPERMATRIX
struct FLASH_Queue_s
{
  // Number of tasks currently in queue
  unsigned int  n_tasks;

  // Pointers to head (front) and tail (back) of queue
  FLASH_Task*   head;
  FLASH_Task*   tail;
};

struct FLASH_Task_s
{
  // Execution information
  int           n_ready;

  // Labels
  int           order;
  int           queue;
  int           height;
  int           thread;
  int           cache;
  FLA_Bool      hit;
      
  // Function pointer
  void*         func;

  // Control tree pointer
  void*         cntl;

  // Name of task
  char*         name;

  // GPU enabled task
  FLA_Bool      enabled_gpu;

  // Integer arguments
  int           n_int_args;
  int*          int_arg;

  // Constant FLA_Obj arguments
  int           n_fla_args;
  FLA_Obj*      fla_arg;

  // Input FLA_Obj arguments
  int           n_input_args;
  FLA_Obj*      input_arg;

  // Output FLA_Obj argument
  int           n_output_args;
  FLA_Obj*      output_arg;

  // Number of blocks within all macroblocks
  int           n_macro_args;

  // Number of write after read dependencies
  int           n_war_args;

  // Dependence information
  int           n_dep_args;
  FLASH_Dep*    dep_arg_head;
  FLASH_Dep*    dep_arg_tail;
  
  // Support for a doubly linked list of tasks
  FLASH_Task*   prev_task;
  FLASH_Task*   next_task;

  // Support for a doubly linked list for wait queue
  FLASH_Task*   prev_wait;
  FLASH_Task*   next_wait;
};

struct FLASH_Dep_s
{
  // Task yielding dependency
  FLASH_Task*   task;

  // Support for linked list of FLASH_Deps
  FLASH_Dep*    next_dep;
};
#endif // FLA_ENABLE_SUPERMATRIX

struct FLASH_Thread_s
{
  // The thread's unique identifier
  int       id;

  // Pointer to variables needed to execute SuperMatrix mechanism
  void*     args;

#if FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  // The thread object. Only needed for the POSIX threads implementation.
  pthread_t pthread_obj;
#endif
};

#endif // FLA_TYPE_DEFS_H
// end FLA_type_defs.h

  // --- Pass-through macros for BLIS ---
  #ifdef FLA_ENABLE_CBLAS_INTERFACES
    #define BLIS1_ENABLE_CBLAS_INTERFACES
  #endif
  #ifdef FLA_ENABLE_WINDOWS_BUILD
    #define BLIS1_ENABLE_WINDOWS_BUILD
  #endif
  #ifdef FLA_ENABLE_UPPERCASE_F77
    #define BLIS1_ENABLE_UPPERCASE_F77
  #endif
  #ifdef FLA_ENABLE_VECTOR_INTRINSICS
    #define BLIS1_ENABLE_VECTOR_INTRINSICS
  #endif

  #define BLIS1_VECTOR_INTRINSIC_TYPE FLA_VECTOR_INTRINSIC_TYPE

#else

  // --- BLIS configuration options ---

  // #define BLIS1_ENABLE_USE_OF_FLA_MALLOC
  // #define BLIS1_ENABLE_CBLAS_INTERFACES
  // #define BLIS1_ENABLE_WINDOWS_BUILD
  // #define BLIS1_ENABLE_UPPERCASE_F77
  // #define BLIS1_ENABLE_VECTOR_INTRINSICS
  //   #define BLIS1_VECTOR_INTRINSIC_TYPE BLIS1_NO_INTRINSICS
  //   #define BLIS1_VECTOR_INTRINSIC_TYPE BLIS1_SSE_INTRINSICS

#endif

// begin blis_macro_defs.h


#ifndef BLIS1_MACRO_DEFS_H
#define BLIS1_MACRO_DEFS_H

// --- Constants ---------------------------------------------------------------

#define BLIS1_NO_INTRINSICS  0
#define BLIS1_SSE_INTRINSICS 3

// --- boolean ---

#undef FALSE
#define FALSE 0

#undef TRUE
#define TRUE 1



// --- Functional macros -------------------------------------------------------

// --- Type-agnostic ---

// min, max, abs

#define bl1_min( a, b )  ( (a) < (b) ? (a) : (b) )
#define bl1_max( a, b )  ( (a) > (b) ? (a) : (b) )
#define bl1_abs( a )     ( (a) <= 0 ? -(a) : (a) )

// fmin, fmax, fabs

#define bl1_fmin( a, b ) bl1_min( a, b )
#define bl1_fmax( a, b ) bl1_max( a, b )
#define bl1_fabs( a )    ( (a) <= 0.0 ? -(a) : (a) )

// fminabs, fmaxabs
#define bl1_fminabs( a, b ) \
\
    bl1_fmin( bl1_fabs( a ), \
              bl1_fabs( b ) )

#define bl1_fmaxabs( a, b ) \
\
    bl1_fmax( bl1_fabs( a ), \
              bl1_fabs( b ) )

// --- Type-dependent ---

// --- neg1 ---

// void bl1_sneg1( float* x );
#define bl1_sneg1( x ) \
*(x)     *= -1.0F;

// void bl1_dneg1( double* x );
#define bl1_dneg1( x ) \
*(x)     *= -1.0;

// void bl1_cneg1( scomplex* x );
#define bl1_cneg1( x ) \
(x)->real *= -1.0F; \
(x)->imag *= -1.0F;

// void bl1_zneg1( dcomplex* x );
#define bl1_zneg1( x ) \
(x)->real *= -1.0; \
(x)->imag *= -1.0;

// --- neg2 ---

// void bl1_sneg2( float* x, float* y );
#define bl1_sneg2( x, y ) \
*(y)      = -1.0F * *(x);

// void bl1_dneg2( double* x, double* y );
#define bl1_dneg2( x, y ) \
*(y)      = -1.0  * *(x);

// void bl1_cneg2( scomplex* x, scomplex* y );
#define bl1_cneg2( x, y ) \
(y)->real = -1.0F * (x)->real; \
(y)->imag = -1.0F * (x)->imag;

// void bl1_zneg2( dcomplex* x, dcomplex* y );
#define bl1_zneg2( x, y ) \
(y)->real = -1.0  * (x)->real; \
(y)->imag = -1.0  * (x)->imag;

// --- sqrte ---

// void bl1_ssqrte( float* alpha, int* error );
#define bl1_ssqrte( alpha, error ) \
if ( *(alpha)      <= 0.0F || isnan( *(alpha) ) ) {  *(error) = FLA_FAILURE; } \
else { *(alpha)      =  ( float ) sqrt( *(alpha) );  *(error) = FLA_SUCCESS; }

// void bl1_dsqrte( double* alpha, int* error );
#define bl1_dsqrte( alpha, error ) \
if ( *(alpha)      <= 0.0 || isnan( *(alpha) ) ) {   *(error) = FLA_FAILURE; } \
else { *(alpha)      = ( double ) sqrt( *(alpha) );  *(error) = FLA_SUCCESS; }

// void bl1_csqrte( scomplex* alpha, int* error );
#define bl1_csqrte( alpha, error ) \
if ( (alpha)->real <= 0.0F || isnan( (alpha)->real) ) \
{                     *(error) = FLA_FAILURE; } \
else { \
(alpha)->real =  ( float ) sqrt( (alpha)->real ); \
(alpha)->imag = 0.0F; *(error) = FLA_SUCCESS; }

// void bl1_zsqrte( dcomplex* alpha, int* error );
#define bl1_zsqrte( alpha, error ) \
if ( (alpha)->real <= 0.0 || isnan( (alpha)->real) )  \
{                     *(error) = FLA_FAILURE; } \
else { \
(alpha)->real = ( double ) sqrt( (alpha)->real ); \
(alpha)->imag = 0.0;  *(error) = FLA_SUCCESS; }

// --- absval2 ---

// void bl1_sabsval2( float* alpha, float* absval );
#define bl1_sabsval2( alpha, absval ) \
*(absval) = ( float ) fabs( ( double ) *(alpha) );

// void bl1_dabsval2( double* alpha, double* absval );
#define bl1_dabsval2( alpha, absval ) \
*(absval) = fabs( *(alpha) );

// void bl1_cabsval2( scomplex* x, scomplex* a );
#define bl1_cabsval2( x, a ) \
{ \
	float  s   = bl1_fmaxabs( (x)->real, (x)->imag ); \
	float  mag = sqrtf( s ) * \
	             sqrtf( ( (x)->real / s ) * (x)->real + \
	                    ( (x)->imag / s ) * (x)->imag ); \
	(a)->real   = mag; \
	(a)->imag   = 0.0F; \
}

// void bl1_csabsval2( scomplex* x, float* a );
#define bl1_csabsval2( x, a ) \
{ \
	float  s   = bl1_fmaxabs( (x)->real, (x)->imag ); \
	float  mag = sqrtf( s ) * \
	             sqrtf( ( (x)->real / s ) * (x)->real + \
	                    ( (x)->imag / s ) * (x)->imag ); \
	*(a)       = mag; \
}

// void bl1_zabsval2( dcomplex* x, dcomplex* a );
#define bl1_zabsval2( x, a ) \
{ \
	double s   = bl1_fmaxabs( (x)->real, (x)->imag ); \
	double mag = sqrt( s ) * \
	             sqrt( ( (x)->real / s ) * (x)->real + \
	                   ( (x)->imag / s ) * (x)->imag ); \
	(a)->real   = mag; \
	(a)->imag   = 0.0; \
}

// void bl1_zdabsval2( dcomplex* x, double* a );
#define bl1_zdabsval2( x, a ) \
{ \
	double s   = bl1_fmaxabs( (x)->real, (x)->imag ); \
	double mag = sqrt( s ) * \
	             sqrt( ( (x)->real / s ) * (x)->real + \
	                   ( (x)->imag / s ) * (x)->imag ); \
	*(a)       = mag; \
}


// --- absqr ---

// void bl1_sabsqr( float* alpha );
#define bl1_sabsqr( alpha ) \
*(alpha) = *(alpha) * *(alpha);

// void bl1_dabsqr( double* alpha );
#define bl1_dabsqr( alpha ) \
*(alpha) = *(alpha) * *(alpha);

// void bl1_cabsqr( scomplex* alpha );
#define bl1_cabsqr( alpha ) \
(alpha)->real = (alpha)->real * (alpha)->real + (alpha)->imag * (alpha)->imag; \
(alpha)->imag = 0.0F;

// void bl1_zabsqr( dcomplex* alpha );
#define bl1_zabsqr( alpha ) \
(alpha)->real = (alpha)->real * (alpha)->real + (alpha)->imag * (alpha)->imag; \
(alpha)->imag = 0.0;

// --- invscals ---

// void bl1_sinvscals( float* a, float* y );
#define bl1_sinvscals( a, y ) \
*(y) = *(y) / *(a);

// void bl1_dinvscals( double* a, double* y );
#define bl1_dinvscals( a, y ) \
*(y) = *(y) / *(a);

// void bl1_csinvscals( float* a, scomplex* y );
#define bl1_csinvscals( a, y ) \
{ \
(y)->real = (y)->real / *(a); \
(y)->imag = (y)->imag / *(a); \
}

// void bl1_cinvscals( scomplex* a, scomplex* y );
#define bl1_cinvscals( a, y ) \
{ \
	float  s     = bl1_fmaxabs( (a)->real, (a)->imag ); \
	float  ar_s  = (a)->real / s; \
	float  ai_s  = (a)->imag / s; \
	float  yrt   = (y)->real; \
	float  temp  = ( ar_s * (a)->real + ai_s * (a)->imag ); \
	(y)->real    = ( (yrt)     * ar_s + (y)->imag * ai_s ) / temp; \
	(y)->imag    = ( (y)->imag * ar_s - (yrt)     * ai_s ) / temp; \
}

// void bl1_zdinvscals( double* a, dcomplex* y );
#define bl1_zdinvscals( a, y ) \
{ \
(y)->real = (y)->real / *(a); \
(y)->imag = (y)->imag / *(a); \
}

// void bl1_zinvscals( dcomplex* a, dcomplex* y );
#define bl1_zinvscals( a, y ) \
{ \
	double s     = bl1_fmaxabs( (a)->real, (a)->imag ); \
	double ar_s  = (a)->real / s; \
	double ai_s  = (a)->imag / s; \
	double yrt   = (y)->real; \
	double temp  = ( ar_s * (a)->real + ai_s * (a)->imag ); \
	(y)->real    = ( (yrt)     * ar_s + (y)->imag * ai_s ) / temp; \
	(y)->imag    = ( (y)->imag * ar_s - (yrt)     * ai_s ) / temp; \
}

// --- div3 ---

// void bl1_sdiv3( float* x, float* y, float* a );
#define bl1_sdiv3( x, y, a ) \
*(a) = *(x) / *(y);

// void bl1_ddiv3( double* x, double* y, double* a );
#define bl1_ddiv3( x, y, a ) \
*(a) = *(x) / *(y);

// void bl1_cdiv3( scomplex* x, scomplex* y, scomplex* a );
// a = x / y;
#define bl1_cdiv3( x, y, a ) \
{ \
	*a = *x; \
	bl1_cinvscals( y, a ); \
}

// void bl1_zdiv3( dcomplex* x, dcomplex* y, dcomplex* a );
#define bl1_zdiv3( x, y, a ) \
{ \
	*a = *x; \
	bl1_zinvscals( y, a ); \
}

// --- add3 ---

// void bl1_sadd3( float* x, float* y, float* a );
#define bl1_sadd3( x, y, a ) \
*(a) = *(x) + *(y);

// void bl1_dadd3( double* x, double* y, double* a );
#define bl1_dadd3( x, y, a ) \
*(a) = *(x) + *(y);

// void bl1_cadd3( scomplex* x, scomplex* y, scomplex* a );
#define bl1_cadd3( x, y, a ) \
{ \
(a)->real = (x)->real + (y)->real; \
(a)->imag = (x)->imag + (y)->imag; \
}

// void bl1_zadd3( dcomplex* x, dcomplex* y, dcomplex* a );
#define bl1_zadd3( x, y, a ) \
{ \
(a)->real = (x)->real + (y)->real; \
(a)->imag = (x)->imag + (y)->imag; \
}

// --- copys ---

// void bl1_scopys( conj1_t conj, float* x, float* y );
#define bl1_scopys( conj, x, y ) \
*(y) = *(x);

// void bl1_dcopys( conj1_t conj, double* x, double* y );
#define bl1_dcopys( conj, x, y ) \
*(y) = *(x);

// void bl1_ccopys( conj1_t conj, scomplex* x, scomplex* y );
#define bl1_ccopys( conj, x, y ) \
*(y) = *(x); \
if ( bl1_is_conj( conj ) ) (y)->imag *= -1.0F;

// void bl1_zcopys( conj1_t conj, dcomplex* x, dcomplex* y );
#define bl1_zcopys( conj, x, y ) \
*(y) = *(x); \
if ( bl1_is_conj( conj ) ) (y)->imag *= -1.0;

// --- scals ---

// void bl1_sscals( float* a, float* y );
#define bl1_sscals( a, y ) \
*(y) = *(a) * *(y);

// void bl1_dscals( double* a, double* y );
#define bl1_dscals( a, y ) \
*(y) = *(a) * *(y);

// void bl1_csscals( float* a, scomplex* y );
#define bl1_csscals( a, y ) \
{ \
(y)->real = *(a) * (y)->real; \
(y)->imag = *(a) * (y)->imag; \
}

// void bl1_cscals( scomplex* a, scomplex* y );
#define bl1_cscals( a, y ) \
{ \
float tempr = (a)->real * (y)->real - (a)->imag * (y)->imag; \
float tempi = (a)->imag * (y)->real + (a)->real * (y)->imag; \
(y)->real = tempr; \
(y)->imag = tempi; \
}

// void bl1_zdscals( double* a, dcomplex* y );
#define bl1_zdscals( a, y ) \
{ \
(y)->real = *(a) * (y)->real; \
(y)->imag = *(a) * (y)->imag; \
}

// void bl1_zscals( dcomplex* a, dcomplex* y );
#define bl1_zscals( a, y ) \
{ \
double tempr = (a)->real * (y)->real - (a)->imag * (y)->imag; \
double tempi = (a)->imag * (y)->real + (a)->real * (y)->imag; \
(y)->real = tempr; \
(y)->imag = tempi; \
}

// --- mult3 ---

// void bl1_smult3( float* x, float* y, float* a );
#define bl1_smult3( x, y, a ) \
*(a) = *(x) * *(y);

// void bl1_dmult3( double* x, double* y, double* a );
#define bl1_dmult3( x, y, a ) \
*(a) = *(x) * *(y);

// void bl1_cmult3( scomplex* x, scomplex* y, scomplex* a );
#define bl1_cmult3( x, y, a ) \
{ \
float tempr = (x)->real * (y)->real - (x)->imag * (y)->imag; \
float tempi = (x)->imag * (y)->real + (x)->real * (y)->imag; \
(a)->real = tempr; \
(a)->imag = tempi; \
}

// void bl1_zmult3( dcomplex* x, dcomplex* y, dcomplex* a );
#define bl1_zmult3( x, y, a ) \
{ \
double tempr = (x)->real * (y)->real - (x)->imag * (y)->imag; \
double tempi = (x)->imag * (y)->real + (x)->real * (y)->imag; \
(a)->real = tempr; \
(a)->imag = tempi; \
}

// --- mult4 ---

// void bl1_smult4( float* alpha, float* x, float* y1, float* y2 );
#define bl1_smult4( alpha, x, y1, y2 ) \
*(y2) = *(y1) + *(alpha) * *(x);

// void bl1_dmult4( double* alpha, double* x, double* y1, double* y2 );
#define bl1_dmult4( alpha, x, y1, y2 ) \
*(y2) = *(y1) + *(alpha) * *(x);

// void bl1_cmult4( scomplex* alpha, scomplex* x, scomplex* y1, scomplex* y2 );
#define bl1_cmult4( alpha, x, y1, y2 ) \
{ \
(y2)->real = (y1)->real + (alpha)->real * (x)->real - (alpha)->imag * (x)->imag; \
(y2)->imag = (y1)->imag + (alpha)->imag * (x)->real + (alpha)->real * (x)->imag; \
}

// void bl1_zmult4( dcomplex* alpha, dcomplex* x, dcomplex* y1, dcomplex* y2 );
#define bl1_zmult4( alpha, x, y1, y2 ) \
{ \
(y2)->real = (y1)->real + (alpha)->real * (x)->real - (alpha)->imag * (x)->imag; \
(y2)->imag = (y1)->imag + (alpha)->imag * (x)->real + (alpha)->real * (x)->imag; \
}

// --- conjs ---

// void bl1_sconjs( float* a );
#define bl1_sconjs( a ) \
;

// void bl1_dconjs( double* a );
#define bl1_dconjs( a ) \
;

// void bl1_cconjs( scomplex* a );
#define bl1_cconjs( a ) \
(a)->imag *= -1.0F;

// void bl1_zconjs( dcomplex* a );
#define bl1_zconjs( a ) \
(a)->imag *= -1.0;

// --- copyconj ---

// void bl1_scopyconj( float* x, float* y );
#define bl1_scopyconj( x, y ) \
*(y) = *(x);

// void bl1_dcopyconj( double* x, double* y );
#define bl1_dcopyconj( x, y ) \
*(y) = *(x);

// void bl1_ccopyconj( scomplex* x, scomplex* y );
#define bl1_ccopyconj( x, y ) \
(y)->real =         (x)->real; \
(y)->imag = -1.0F * (x)->imag;

// void bl1_zcopyconj( dcomplex* x, dcomplex* y );
#define bl1_zcopyconj( x, y ) \
(y)->real =         (x)->real; \
(y)->imag = -1.0  * (x)->imag;

// --- eq1 ---

// void bl1_seq1( float* alpha );
#define bl1_seq1( alpha ) \
  ( *alpha == 1.0F )

// void bl1_deq1( double* alpha );
#define bl1_deq1( alpha ) \
  ( *alpha == 1.0 )

// void bl1_ceq1( scomplex* alpha );
#define bl1_ceq1( alpha ) \
  ( (alpha)->real == 1.0F && (alpha)->imag == 0.0F )

// void bl1_zeq1( dcomplex* alpha );
#define bl1_zeq1( alpha ) \
  ( (alpha)->real == 1.0 && (alpha)->imag == 0.0 )

// --- Swapping/toggle macros --------------------------------------------------

// --- swap_pointers ---

#define bl1_sswap_pointers( a, b ) \
{ \
float* temp = (a); \
(a) = (b); \
(b) = temp; \
}

#define bl1_dswap_pointers( a, b ) \
{ \
double* temp = (a); \
(a) = (b); \
(b) = temp; \
}

#define bl1_cswap_pointers( a, b ) \
{ \
void* temp = (a); \
(a) = (b); \
(b) = temp; \
}

#define bl1_zswap_pointers( a, b ) \
{ \
void* temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- swap_ints ---

#define bl1_swap_ints( a, b ) \
{ \
int temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- swap_trans ---

#define bl1_swap_trans( a, b ) \
{ \
trans1_t temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- swap_conj ---

#define bl1_swap_conj( a, b ) \
{ \
conj1_t temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- toggle_side ---

#define bl1_toggle_side( side ) \
{ \
if ( bl1_is_left( side ) ) side = BLIS1_RIGHT; \
else                       side = BLIS1_LEFT; \
}

// --- toggle_uplo ---

#define bl1_toggle_uplo( uplo ) \
{ \
if ( bl1_is_lower( uplo ) ) uplo = BLIS1_UPPER_TRIANGULAR; \
else                        uplo = BLIS1_LOWER_TRIANGULAR; \
}

// --- toggle_trans ---
#define bl1_toggle_trans( trans ) \
{ \
if      ( bl1_is_notrans( trans ) )     trans = BLIS1_TRANSPOSE; \
else if ( bl1_is_trans( trans ) )       trans = BLIS1_NO_TRANSPOSE; \
else if ( bl1_is_conjnotrans( trans ) ) trans = BLIS1_CONJ_TRANSPOSE; \
else                                    trans = BLIS1_CONJ_NO_TRANSPOSE; \
}

// --- toggle_conjtrans ---
#define bl1_toggle_conjtrans( trans ) \
{ \
if      ( bl1_is_notrans( trans ) )     trans = BLIS1_CONJ_TRANSPOSE; \
else                                    trans = BLIS1_NO_TRANSPOSE; \
}

// --- toggle_conj ---

#define bl1_toggle_conj( conj ) \
{ \
if ( bl1_is_conj( conj ) ) conj = BLIS1_NO_CONJUGATE; \
else                       conj = BLIS1_CONJUGATE; \
}

#endif // #ifndef BLIS1_MACRO_DEFS_H
// end blis_macro_defs.h
// begin blis_type_defs.h


#ifndef BLIS1_TYPE_DEFS_H
#define BLIS1_TYPE_DEFS_H

// --- Basic type definitions -------------------------------------------------



#define BLIS1_TRANS_BEGIN 100
#define BLIS1_UPLO_BEGIN  200
#define BLIS1_SIDE_BEGIN  300
#define BLIS1_DIAG_BEGIN  400
#define BLIS1_CONJ_BEGIN  500

typedef enum
{
	BLIS1_NO_TRANSPOSE = BLIS1_TRANS_BEGIN,
	BLIS1_TRANSPOSE,
	BLIS1_CONJ_NO_TRANSPOSE,
	BLIS1_CONJ_TRANSPOSE
} trans1_t;

typedef enum
{
	BLIS1_LOWER_TRIANGULAR = BLIS1_UPLO_BEGIN,
	BLIS1_UPPER_TRIANGULAR
} uplo1_t;

typedef enum
{
	BLIS1_LEFT = BLIS1_SIDE_BEGIN,
	BLIS1_RIGHT
} side1_t;

typedef enum
{
	BLIS1_NONUNIT_DIAG = BLIS1_DIAG_BEGIN,
	BLIS1_UNIT_DIAG,
	BLIS1_ZERO_DIAG
} diag1_t;

typedef enum
{
	BLIS1_NO_CONJUGATE = BLIS1_CONJ_BEGIN,
	BLIS1_CONJUGATE
} conj1_t;





// --- Intrinsic/assembly definitions ----------------------------------------





// Only define vector intrinsics types if they are not already provided by
// libflame.
#ifndef BLIS1_FROM_LIBFLAME

#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS

#include "pmmintrin.h" // skipped
typedef union
{
    __m128d v; 
    double  d[2];
} v2df_t;
#endif

#endif


// --- Complex type definitions -----------------------------------------------

// Only define complex types if they are not already provided by libflame.
//#ifndef BLIS1_ENABLE_USE_OF_LIBFLAME_TYPES
#ifndef BLIS1_FROM_LIBFLAME

typedef struct scomplex
{
  float real, imag;
} scomplex;

typedef struct dcomplex
{
  double real, imag;
} dcomplex;

#endif


#endif // BLIS1_TYPE_DEFS_H
// end blis_type_defs.h

// begin blis_prototypes_util.h


// --- Utility-level BLAS-like prototypes --------------------------------------

// --- constant-generating functions ---

float    bl1_s2( void );
double   bl1_d2( void );
scomplex bl1_c2( void );
dcomplex bl1_z2( void );
float    bl1_s1( void );
double   bl1_d1( void );
scomplex bl1_c1( void );
dcomplex bl1_z1( void );
float    bl1_s1h( void );
double   bl1_d1h( void );
scomplex bl1_c1h( void );
dcomplex bl1_z1h( void );
float    bl1_s0( void );
double   bl1_d0( void );
scomplex bl1_c0( void );
dcomplex bl1_z0( void );
float    bl1_sm1h( void );
double   bl1_dm1h( void );
scomplex bl1_cm1h( void );
dcomplex bl1_zm1h( void );
float    bl1_sm1( void );
double   bl1_dm1( void );
scomplex bl1_cm1( void );
dcomplex bl1_zm1( void );
float    bl1_sm2( void );
double   bl1_dm2( void );
scomplex bl1_cm2( void );
dcomplex bl1_zm2( void );

// --- allocv ---

void*     bl1_vallocv( unsigned int n_elem, unsigned int elem_size );
int*      bl1_iallocv( unsigned int n_elem );
float*    bl1_sallocv( unsigned int n_elem );
double*   bl1_dallocv( unsigned int n_elem );
scomplex* bl1_callocv( unsigned int n_elem );
dcomplex* bl1_zallocv( unsigned int n_elem );

// --- allocm ---

void*     bl1_vallocm( unsigned int m, unsigned int n, unsigned int elem_size );
int*      bl1_iallocm( unsigned int m, unsigned int n );
float*    bl1_sallocm( unsigned int m, unsigned int n );
double*   bl1_dallocm( unsigned int m, unsigned int n );
scomplex* bl1_callocm( unsigned int m, unsigned int n );
dcomplex* bl1_zallocm( unsigned int m, unsigned int n );

// --- apdiagmv ---

void bl1_sapdiagmv( side1_t side, conj1_t conj, int m, int n, float*    x, int incx, float*    a, int a_rs, int a_cs );
void bl1_dapdiagmv( side1_t side, conj1_t conj, int m, int n, double*   x, int incx, double*   a, int a_rs, int a_cs );
void bl1_csapdiagmv( side1_t side, conj1_t conj, int m, int n, float*    x, int incx, scomplex* a, int a_rs, int a_cs );
void bl1_capdiagmv( side1_t side, conj1_t conj, int m, int n, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs );
void bl1_zdapdiagmv( side1_t side, conj1_t conj, int m, int n, double*   x, int incx, dcomplex* a, int a_rs, int a_cs );
void bl1_zapdiagmv( side1_t side, conj1_t conj, int m, int n, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs );

// --- create_contigm ---

void bl1_screate_contigm( int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dcreate_contigm( int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_ccreate_contigm( int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zcreate_contigm( int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- create_contigmt ---

void bl1_screate_contigmt( trans1_t trans_dims, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dcreate_contigmt( trans1_t trans_dims, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_ccreate_contigmt( trans1_t trans_dims, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zcreate_contigmt( trans1_t trans_dims, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- create_contigmr ---

void bl1_screate_contigmr( uplo1_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dcreate_contigmr( uplo1_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_ccreate_contigmr( uplo1_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zcreate_contigmr( uplo1_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- create_contigmsr ---

void bl1_screate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dcreate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_ccreate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zcreate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_contigm ---

void bl1_sfree_contigm( float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dfree_contigm( double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_cfree_contigm( scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zfree_contigm( dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_saved_contigm ---

void bl1_sfree_saved_contigm( int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dfree_saved_contigm( int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_cfree_saved_contigm( int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zfree_saved_contigm( int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_saved_contigmr ---

void bl1_sfree_saved_contigmr( uplo1_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dfree_saved_contigmr( uplo1_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_cfree_saved_contigmr( uplo1_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zfree_saved_contigmr( uplo1_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_saved_contigmsr ---

void bl1_sfree_saved_contigmsr( side1_t side, uplo1_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dfree_saved_contigmsr( side1_t side, uplo1_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_cfree_saved_contigmsr( side1_t side, uplo1_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zfree_saved_contigmsr( side1_t side, uplo1_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- ewinvscalv ---

void bl1_sewinvscalv( conj1_t conj, int n, float*    x, int incx, float*    y, int incy );
void bl1_dewinvscalv( conj1_t conj, int n, double*   x, int incx, double*   y, int incy );
void bl1_csewinvscalv( conj1_t conj, int n, float*    x, int incx, scomplex* y, int incy );
void bl1_cewinvscalv( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zdewinvscalv( conj1_t conj, int n, double*   x, int incx, dcomplex* y, int incy );
void bl1_zewinvscalv( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- ewscalmt ---

void bl1_sewinvscalmt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dewinvscalmt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_csewinvscalmt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cewinvscalmt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zdewinvscalmt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zewinvscalmt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- ewscalv ---

void bl1_sewscalv( conj1_t conj, int n, float*    x, int incx, float*    y, int incy );
void bl1_dewscalv( conj1_t conj, int n, double*   x, int incx, double*   y, int incy );
void bl1_csewscalv( conj1_t conj, int n, float*    x, int incx, scomplex* y, int incy );
void bl1_cewscalv( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zdewscalv( conj1_t conj, int n, double*   x, int incx, dcomplex* y, int incy );
void bl1_zewscalv( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- ewscalmt ---

void bl1_sewscalmt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dewscalmt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_csewscalmt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cewscalmt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zdewscalmt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zewscalmt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- free ---

void bl1_vfree( void*     p );
void bl1_ifree( int*      p );
void bl1_sfree( float*    p );
void bl1_dfree( double*   p );
void bl1_cfree( scomplex* p );
void bl1_zfree( dcomplex* p );

// --- inverts ---

void bl1_sinverts( conj1_t conj, float*    alpha );
void bl1_dinverts( conj1_t conj, double*   alpha );
void bl1_cinverts( conj1_t conj, scomplex* alpha );
void bl1_zinverts( conj1_t conj, dcomplex* alpha );

// --- invert2s ---

void bl1_sinvert2s( conj1_t conj, float*    alpha, float*    beta );
void bl1_dinvert2s( conj1_t conj, double*   alpha, double*   beta );
void bl1_cinvert2s( conj1_t conj, scomplex* alpha, scomplex* beta );
void bl1_zinvert2s( conj1_t conj, dcomplex* alpha, dcomplex* beta );

// --- invertv ---

void bl1_sinvertv( conj1_t conj, int n, float*    x, int incx );
void bl1_dinvertv( conj1_t conj, int n, double*   x, int incx );
void bl1_cinvertv( conj1_t conj, int n, scomplex* x, int incx );
void bl1_zinvertv( conj1_t conj, int n, dcomplex* x, int incx );

// --- ident ---

void bl1_sident( int m, float*    a, int a_rs, int a_cs );
void bl1_dident( int m, double*   a, int a_rs, int a_cs );
void bl1_cident( int m, scomplex* a, int a_rs, int a_cs );
void bl1_zident( int m, dcomplex* a, int a_rs, int a_cs );

// --- maxabsv ---

void bl1_smaxabsv( int n, float*    x, int incx, float*  maxabs );
void bl1_dmaxabsv( int n, double*   x, int incx, double* maxabs );
void bl1_cmaxabsv( int n, scomplex* x, int incx, float*  maxabs );
void bl1_zmaxabsv( int n, dcomplex* x, int incx, double* maxabs );

// --- maxabsm ---

void bl1_smaxabsm( int m, int n, float*    a, int a_rs, int a_cs, float*  maxabs );
void bl1_dmaxabsm( int m, int n, double*   a, int a_rs, int a_cs, double* maxabs );
void bl1_cmaxabsm( int m, int n, scomplex* a, int a_rs, int a_cs, float*  maxabs );
void bl1_zmaxabsm( int m, int n, dcomplex* a, int a_rs, int a_cs, double* maxabs );

// --- maxabsmr ---

void bl1_smaxabsmr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, float*  maxabs );
void bl1_dmaxabsmr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, double* maxabs );
void bl1_cmaxabsmr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, float*  maxabs );
void bl1_zmaxabsmr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, double* maxabs );

// --- rands ---

void bl1_srands( float*    alpha );
void bl1_drands( double*   alpha );
void bl1_crands( scomplex* alpha );
void bl1_zrands( dcomplex* alpha );

// --- randv ---

void bl1_srandv( int n, float*    x, int incx );
void bl1_drandv( int n, double*   x, int incx );
void bl1_crandv( int n, scomplex* x, int incx );
void bl1_zrandv( int n, dcomplex* x, int incx );

// --- randm ---

void bl1_srandm( int m, int n, float*    a, int a_rs, int a_cs );
void bl1_drandm( int m, int n, double*   a, int a_rs, int a_cs );
void bl1_crandm( int m, int n, scomplex* a, int a_rs, int a_cs );
void bl1_zrandm( int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- randmr ---
void bl1_srandmr( uplo1_t uplo, diag1_t diag, int m, int n, float*    a, int a_rs, int a_cs );
void bl1_drandmr( uplo1_t uplo, diag1_t diag, int m, int n, double*   a, int a_rs, int a_cs );
void bl1_crandmr( uplo1_t uplo, diag1_t diag, int m, int n, scomplex* a, int a_rs, int a_cs );
void bl1_zrandmr( uplo1_t uplo, diag1_t diag, int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- set_contig_strides ---

void bl1_set_contig_strides( int m, int n, int* rs, int* cs );

// --- set_dims_with_side ---

void bl1_set_dim_with_side( side1_t side, int m, int n, int* dim_new );

// --- set_dims_with_trans ---

void bl1_set_dims_with_trans( trans1_t trans, int m, int n, int* m_new, int* n_new );

// --- setv ---

void bl1_isetv( int m, int*      sigma, int*      x, int incx );
void bl1_ssetv( int m, float*    sigma, float*    x, int incx );
void bl1_dsetv( int m, double*   sigma, double*   x, int incx );
void bl1_csetv( int m, scomplex* sigma, scomplex* x, int incx );
void bl1_zsetv( int m, dcomplex* sigma, dcomplex* x, int incx );

// --- setm ---

void bl1_isetm( int m, int n, int*      sigma, int*      a, int a_rs, int a_cs );
void bl1_ssetm( int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bl1_dsetm( int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bl1_csetm( int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zsetm( int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );

// --- setmr ---

void bl1_ssetmr( uplo1_t uplo, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bl1_dsetmr( uplo1_t uplo, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bl1_csetmr( uplo1_t uplo, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zsetmr( uplo1_t uplo, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );

// --- setdiag ---

void bl1_isetdiag( int offset, int m, int n, int*      sigma, int*      a, int a_rs, int a_cs );
void bl1_ssetdiag( int offset, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bl1_dsetdiag( int offset, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bl1_csetdiag( int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zsetdiag( int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );

// --- scalediag ---

void bl1_sscalediag( conj1_t conj, int offset, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bl1_dscalediag( conj1_t conj, int offset, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bl1_cscalediag( conj1_t conj, int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zscalediag( conj1_t conj, int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );
void bl1_csscalediag( conj1_t conj, int offset, int m, int n, float*    sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zdscalediag( conj1_t conj, int offset, int m, int n, double*   sigma, dcomplex* a, int a_rs, int a_cs );

// --- shiftdiag ---

void bl1_sshiftdiag( conj1_t conj, int offset, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bl1_dshiftdiag( conj1_t conj, int offset, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bl1_cshiftdiag( conj1_t conj, int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zshiftdiag( conj1_t conj, int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );
void bl1_csshiftdiag( conj1_t conj, int offset, int m, int n, float*    sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zdshiftdiag( conj1_t conj, int offset, int m, int n, double*   sigma, dcomplex* a, int a_rs, int a_cs );

// --- symmize ---

void bl1_ssymmize( conj1_t conj, uplo1_t uplo, int m, float*    a, int a_rs, int a_cs );
void bl1_dsymmize( conj1_t conj, uplo1_t uplo, int m, double*   a, int a_rs, int a_cs );
void bl1_csymmize( conj1_t conj, uplo1_t uplo, int m, scomplex* a, int a_rs, int a_cs );
void bl1_zsymmize( conj1_t conj, uplo1_t uplo, int m, dcomplex* a, int a_rs, int a_cs );

// end blis_prototypes_util.h
// begin blis_prototypes_query.h


// --- Query routine prototypes ------------------------------------------------

// --- trans ---

int bl1_does_trans( trans1_t trans );
int bl1_does_notrans( trans1_t trans );
int bl1_does_conj( trans1_t trans );

int bl1_is_notrans( trans1_t trans );
int bl1_is_trans( trans1_t trans );
int bl1_is_conjnotrans( trans1_t trans );
int bl1_is_conjtrans( trans1_t trans );

// --- conj ---

int bl1_is_noconj( conj1_t conj );
int bl1_is_conj( conj1_t conj );

// --- uplo ---

int bl1_is_lower( uplo1_t uplo );
int bl1_is_upper( uplo1_t uplo );

// --- side ---

int bl1_is_left( side1_t side );
int bl1_is_right( side1_t side );

// --- diag ---

int bl1_is_nonunit_diag( diag1_t diag );
int bl1_is_unit_diag( diag1_t diag );
int bl1_is_zero_diag( diag1_t diag );

// --- mapping-related ---

conj1_t bl1_proj_trans1_to_conj( trans1_t trans );

// --- storage-related ---

void bl1_check_storage_3m( int a_rs, int a_cs, int b_rs, int b_cs, int c_rs, int c_cs );
void bl1_check_storage_2m( int a_rs, int a_cs, int b_rs, int b_cs );
int bl1_is_row_or_col_storage( int rs, int cs );
int bl1_is_row_storage( int rs, int cs );
int bl1_is_col_storage( int rs, int cs );
int bl1_is_gen_storage( int rs, int cs );
int bl1_is_vector( int m, int n );

// --- vector-related ---

int bl1_vector_dim( int m, int n );
int bl1_vector_inc( trans1_t trans, int m, int n, int rs, int cs );

// --- dimension-related ---

int bl1_zero_dim1( int m );
int bl1_zero_dim2( int m, int n );
int bl1_zero_dim3( int m, int k, int n );

// end blis_prototypes_query.h
// begin blis_prototypes_misc.h


// --- Abort prototypes --------------------------------------------------------

void bl1_abort( void );
void bl1_abort_msg( char* message );

// --- Parameter-mapping prototypes --------------------------------------------

void bl1_param_map_to_netlib_trans( trans1_t blis_trans, void* blas_trans );
void bl1_param_map_to_netlib_uplo(  uplo1_t  blis_uplo,  void* blas_uplo );
void bl1_param_map_to_netlib_side(  side1_t  blis_side,  void* blas_side );
void bl1_param_map_to_netlib_diag(  diag1_t  blis_diag,  void* blas_diag );

// end blis_prototypes_misc.h

// begin blis_prototypes_level1.h


// --- Level-1 BLAS-like prototypes --------------------------------------------

// --- amax ---

void bl1_samax( int n, float*    x, int incx, int* index );
void bl1_damax( int n, double*   x, int incx, int* index );
void bl1_camax( int n, scomplex* x, int incx, int* index );
void bl1_zamax( int n, dcomplex* x, int incx, int* index );

// --- asum ---

void bl1_sasum( int n, float*    x, int incx, float*  norm );
void bl1_dasum( int n, double*   x, int incx, double* norm );
void bl1_casum( int n, scomplex* x, int incx, float*  norm );
void bl1_zasum( int n, dcomplex* x, int incx, double* norm );

// --- axpy ---

void bl1_saxpy( int n, float*    alpha, float*    x, int incx, float*    y, int incy );
void bl1_daxpy( int n, double*   alpha, double*   x, int incx, double*   y, int incy );
void bl1_caxpy( int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zaxpy( int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy );

// --- axpyv ---

void bl1_saxpyv( conj1_t conj, int n, float*    alpha, float*    x, int incx, float*    y, int incy );
void bl1_daxpyv( conj1_t conj, int n, double*   alpha, double*   x, int incx, double*   y, int incy );
void bl1_caxpyv( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zaxpyv( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy );

// --- axpymt ---

void bl1_saxpymt( trans1_t trans, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_daxpymt( trans1_t trans, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_caxpymt( trans1_t trans, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zaxpymt( trans1_t trans, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- axpymrt ---

void bl1_saxpymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_daxpymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_caxpymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zaxpymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- axpysv ---

void bl1_saxpysv( int n, float*    alpha0, float*    alpha1, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_daxpysv( int n, double*   alpha0, double*   alpha1, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_caxpysv( int n, scomplex* alpha0, scomplex* alpha1, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zaxpysv( int n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- axpysmt ---

void bl1_saxpysmt( trans1_t trans, int m, int n, float*    alpha0, float*    alpha1, float*    a, int a_rs, int a_cs, float*    beta, float*    b, int b_rs, int b_cs );
void bl1_daxpysmt( trans1_t trans, int m, int n, double*   alpha0, double*   alpha1, double*   a, int a_rs, int a_cs, double*   beta, double*   b, int b_rs, int b_cs );
void bl1_caxpysmt( trans1_t trans, int m, int n, scomplex* alpha0, scomplex* alpha1, scomplex* a, int a_rs, int a_cs, scomplex* beta, scomplex* b, int b_rs, int b_cs );
void bl1_zaxpysmt( trans1_t trans, int m, int n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* a, int a_rs, int a_cs, dcomplex* beta, dcomplex* b, int b_rs, int b_cs );

// --- conjv ---

void bl1_sconjv( int m, float* x, int incx );
void bl1_dconjv( int m, double* x, int incx );
void bl1_cconjv( int m, scomplex* x, int incx );
void bl1_zconjv( int m, dcomplex* x, int incx );

// --- conjm ---

void bl1_sconjm( int m, int n, float*    a, int a_rs, int a_cs );
void bl1_dconjm( int m, int n, double*   a, int a_rs, int a_cs );
void bl1_cconjm( int m, int n, scomplex* a, int a_rs, int a_cs );
void bl1_zconjm( int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- conjmr ---

void bl1_sconjmr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs );
void bl1_dconjmr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs );
void bl1_cconjmr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs );
void bl1_zconjmr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- copy ---

void bl1_scopy( int m, float*    x, int incx, float*    y, int incy );
void bl1_dcopy( int m, double*   x, int incx, double*   y, int incy );
void bl1_ccopy( int m, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zcopy( int m, dcomplex* x, int incx, dcomplex* y, int incy );

// --- copyv ---

void bl1_icopyv( conj1_t conj, int m, int*      x, int incx, int*      y, int incy );
void bl1_scopyv( conj1_t conj, int m, float*    x, int incx, float*    y, int incy );
void bl1_dcopyv( conj1_t conj, int m, double*   x, int incx, double*   y, int incy );
void bl1_ccopyv( conj1_t conj, int m, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zcopyv( conj1_t conj, int m, dcomplex* x, int incx, dcomplex* y, int incy );

void bl1_sdcopyv( conj1_t conj, int m, float*    x, int incx, double*   y, int incy );
void bl1_dscopyv( conj1_t conj, int m, double*   x, int incx, float*    y, int incy );
void bl1_sccopyv( conj1_t conj, int m, float*    x, int incx, scomplex* y, int incy );
void bl1_cscopyv( conj1_t conj, int m, scomplex* x, int incx, float*    y, int incy );
void bl1_szcopyv( conj1_t conj, int m, float*    x, int incx, dcomplex* y, int incy );
void bl1_zscopyv( conj1_t conj, int m, dcomplex* x, int incx, float*    y, int incy );
void bl1_dccopyv( conj1_t conj, int m, double*   x, int incx, scomplex* y, int incy );
void bl1_cdcopyv( conj1_t conj, int m, scomplex* x, int incx, double*   y, int incy );
void bl1_dzcopyv( conj1_t conj, int m, double*   x, int incx, dcomplex* y, int incy );
void bl1_zdcopyv( conj1_t conj, int m, dcomplex* x, int incx, double*   y, int incy );
void bl1_czcopyv( conj1_t conj, int m, scomplex* x, int incx, dcomplex* y, int incy );
void bl1_zccopyv( conj1_t conj, int m, dcomplex* x, int incx, scomplex* y, int incy );

// --- copymr ---

void bl1_scopymr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dcopymr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_ccopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zcopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bl1_sscopymr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_sdcopymr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dscopymr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_sccopymr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cscopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_szcopymr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zscopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_ddcopymr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dccopymr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cdcopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dzcopymr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zdcopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_cccopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_czcopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zccopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zzcopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- copymrt ---

void bl1_scopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_ccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bl1_sscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_sdcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_sccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_szcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_dscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_ddcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_dzcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_cscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_cdcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_cccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_czcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_zdcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_zccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zzcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- copymt ---

void bl1_icopymt( trans1_t trans, int m, int n, int*      a, int a_rs, int a_cs, int*      b, int b_rs, int b_cs );
void bl1_scopymt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dcopymt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_ccopymt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zcopymt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bl1_sscopymt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_sdcopymt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dscopymt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_sccopymt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cscopymt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_szcopymt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zscopymt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_ddcopymt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dccopymt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cdcopymt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dzcopymt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zdcopymt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_cccopymt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_czcopymt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zccopymt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zzcopymt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- dot ---

void bl1_cdot_in( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy, scomplex* rho );
void bl1_zdot_in( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* rho );

void bl1_sdot( conj1_t conj, int n, float*    x, int incx, float*    y, int incy, float*    rho );
void bl1_ddot( conj1_t conj, int n, double*   x, int incx, double*   y, int incy, double*   rho );
void bl1_cdot( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy, scomplex* rho );
void bl1_zdot( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* rho );

// --- dots ---

void bl1_sdots( conj1_t conj, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    beta, float*    rho );
void bl1_ddots( conj1_t conj, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   beta, double*   rho );
void bl1_cdots( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* beta, scomplex* rho );
void bl1_zdots( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* beta, dcomplex* rho );

// --- dot2s ---

void bl1_sdot2s( conj1_t conj, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    beta, float*    rho );
void bl1_ddot2s( conj1_t conj, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   beta, double*   rho );
void bl1_cdot2s( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* beta, scomplex* rho );
void bl1_zdot2s( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* beta, dcomplex* rho );

// --- fnorm ---

void bl1_sfnorm( int m, int n, float*    a, int a_rs, int a_cs, float*  norm );
void bl1_dfnorm( int m, int n, double*   a, int a_rs, int a_cs, double* norm );
void bl1_cfnorm( int m, int n, scomplex* a, int a_rs, int a_cs, float*  norm );
void bl1_zfnorm( int m, int n, dcomplex* a, int a_rs, int a_cs, double* norm );

// --- invscalv ---

void bl1_sinvscalv(  conj1_t conj, int n, float*    alpha, float*    x, int incx );
void bl1_dinvscalv(  conj1_t conj, int n, double*   alpha, double*   x, int incx );
void bl1_csinvscalv( conj1_t conj, int n, float*    alpha, scomplex* x, int incx );
void bl1_cinvscalv(  conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx );
void bl1_zdinvscalv( conj1_t conj, int n, double*   alpha, dcomplex* x, int incx );
void bl1_zinvscalv(  conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx );

// --- invscalm ---

void bl1_sinvscalm(  conj1_t conj, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs );
void bl1_dinvscalm(  conj1_t conj, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs );
void bl1_csinvscalm( conj1_t conj, int m, int n, float*    alpha, scomplex* a, int a_rs, int a_cs );
void bl1_cinvscalm(  conj1_t conj, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs );
void bl1_zdinvscalm( conj1_t conj, int m, int n, double*   alpha, dcomplex* a, int a_rs, int a_cs );
void bl1_zinvscalm(  conj1_t conj, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs );

// --- nrm2 ---

void bl1_snrm2( int n, float*    x, int incx, float*  norm );
void bl1_dnrm2( int n, double*   x, int incx, double* norm );
void bl1_cnrm2( int n, scomplex* x, int incx, float*  norm );
void bl1_znrm2( int n, dcomplex* x, int incx, double* norm );

// --- scal ---

void bl1_sscal(  int n, float*    alpha, float*    x, int incx );
void bl1_dscal(  int n, double*   alpha, double*   x, int incx );
void bl1_csscal( int n, float*    alpha, scomplex* x, int incx );
void bl1_cscal(  int n, scomplex* alpha, scomplex* x, int incx );
void bl1_zdscal( int n, double*   alpha, dcomplex* x, int incx );
void bl1_zscal(  int n, dcomplex* alpha, dcomplex* x, int incx );

// --- scalv ---

void bl1_sscalv(  conj1_t conj, int n, float*    alpha, float*    x, int incx );
void bl1_dscalv(  conj1_t conj, int n, double*   alpha, double*   x, int incx );
void bl1_csscalv( conj1_t conj, int n, float*    alpha, scomplex* x, int incx );
void bl1_cscalv(  conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx );
void bl1_zdscalv( conj1_t conj, int n, double*   alpha, dcomplex* x, int incx );
void bl1_zscalv(  conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx );

// --- scalm ---

void bl1_sscalm(  conj1_t conj, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs );
void bl1_dscalm(  conj1_t conj, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs );
void bl1_csscalm( conj1_t conj, int m, int n, float*    alpha, scomplex* a, int a_rs, int a_cs );
void bl1_cscalm(  conj1_t conj, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs );
void bl1_zdscalm( conj1_t conj, int m, int n, double*   alpha, dcomplex* a, int a_rs, int a_cs );
void bl1_zscalm(  conj1_t conj, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs );

// --- scalmr ---

void bl1_sscalmr(  uplo1_t uplo, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs );
void bl1_dscalmr(  uplo1_t uplo, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs );
void bl1_csscalmr( uplo1_t uplo, int m, int n, float*    alpha, scomplex* a, int a_rs, int a_cs );
void bl1_cscalmr(  uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs );
void bl1_zdscalmr( uplo1_t uplo, int m, int n, double*   alpha, dcomplex* a, int a_rs, int a_cs );
void bl1_zscalmr(  uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs );

// --- swap ---

void bl1_sswap( int n, float*    x, int incx, float*    y, int incy );
void bl1_dswap( int n, double*   x, int incx, double*   y, int incy );
void bl1_cswap( int n, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zswap( int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- swapv ---

void bl1_sswapv( int n, float*    x, int incx, float*    y, int incy );
void bl1_dswapv( int n, double*   x, int incx, double*   y, int incy );
void bl1_cswapv( int n, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zswapv( int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- swapmt ---

void bl1_sswapmt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dswapmt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_cswapmt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zswapmt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// end blis_prototypes_level1.h
// begin blis_prototypes_level2.h


// --- Level-2 BLAS-like prototypes --------------------------------------------

// --- gemv ---

void bl1_sgemv( trans1_t transa, conj1_t conjx, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_dgemv( trans1_t transa, conj1_t conjx, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_cgemv( trans1_t transa, conj1_t conjx, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zgemv( trans1_t transa, conj1_t conjx, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

void bl1_sgemv_blas( trans1_t transa, int m, int n, float*    alpha, float*    a, int lda, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_dgemv_blas( trans1_t transa, int m, int n, double*   alpha, double*   a, int lda, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_cgemv_blas( trans1_t transa, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zgemv_blas( trans1_t transa, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- ger ---

void bl1_sger( conj1_t conjx, conj1_t conjy, int m, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int a_rs, int a_cs );
void bl1_dger( conj1_t conjx, conj1_t conjy, int m, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int a_rs, int a_cs );
void bl1_cger( conj1_t conjx, conj1_t conjy, int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs );
void bl1_zger( conj1_t conjx, conj1_t conjy, int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs );

void bl1_sger_blas(  int m, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int lda );
void bl1_dger_blas(  int m, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int lda );
void bl1_cgerc_blas( int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bl1_cgeru_blas( int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bl1_zgerc_blas( int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );
void bl1_zgeru_blas( int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );

// --- hemv ---

void bl1_shemv( uplo1_t uplo, conj1_t conj, int m, float*    alpha, float*    a, int a_rs, int a_cs, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_dhemv( uplo1_t uplo, conj1_t conj, int m, double*   alpha, double*   a, int a_rs, int a_cs, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_chemv( uplo1_t uplo, conj1_t conj, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zhemv( uplo1_t uplo, conj1_t conj, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

void bl1_chemv_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zhemv_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- her ---

void bl1_sher( uplo1_t uplo, conj1_t conj, int m, float*  alpha, float*    x, int incx, float*    a, int a_rs, int a_cs );
void bl1_dher( uplo1_t uplo, conj1_t conj, int m, double* alpha, double*   x, int incx, double*   a, int a_rs, int a_cs );
void bl1_cher( uplo1_t uplo, conj1_t conj, int m, float*  alpha, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs );
void bl1_zher( uplo1_t uplo, conj1_t conj, int m, double* alpha, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs );

void bl1_cher_blas( uplo1_t uplo, int m, float*  alpha, scomplex* x, int incx, scomplex* a, int lda );
void bl1_zher_blas( uplo1_t uplo, int m, double* alpha, dcomplex* x, int incx, dcomplex* a, int lda );

// --- her2 ---

void bl1_sher2( uplo1_t uplo, conj1_t conj, int m, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int a_rs, int a_cs );
void bl1_dher2( uplo1_t uplo, conj1_t conj, int m, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int a_rs, int a_cs );
void bl1_cher2( uplo1_t uplo, conj1_t conj, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs );
void bl1_zher2( uplo1_t uplo, conj1_t conj, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs );

void bl1_cher2_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bl1_zher2_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );

// --- symv ---

void bl1_ssymv( uplo1_t uplo, int m, float*    alpha, float*    a, int a_rs, int a_cs, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_dsymv( uplo1_t uplo, int m, double*   alpha, double*   a, int a_rs, int a_cs, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_csymv( uplo1_t uplo, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zsymv( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

void bl1_ssymv_blas( uplo1_t uplo, int m, float*    alpha, float*    a, int lda, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_dsymv_blas( uplo1_t uplo, int m, double*   alpha, double*   a, int lda, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_csymv_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zsymv_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- syr ---

void bl1_ssyr( uplo1_t uplo, int m, float*    alpha, float*    x, int incx, float*    a, int a_rs, int a_cs );
void bl1_dsyr( uplo1_t uplo, int m, double*   alpha, double*   x, int incx, double*   a, int a_rs, int a_cs );
void bl1_csyr( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs );
void bl1_zsyr( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs );

void bl1_ssyr_blas( uplo1_t uplo, int m, float*    alpha, float*    x, int incx, float*    a, int lda );
void bl1_dsyr_blas( uplo1_t uplo, int m, double*   alpha, double*   x, int incx, double*   a, int lda );
void bl1_csyr_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* a, int lda );
void bl1_zsyr_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* a, int lda );

// --- syr2 ---

void bl1_ssyr2( uplo1_t uplo, int m, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int a_rs, int a_cs );
void bl1_dsyr2( uplo1_t uplo, int m, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int a_rs, int a_cs );
void bl1_csyr2( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs );
void bl1_zsyr2( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs );

void bl1_ssyr2_blas( uplo1_t uplo, int m, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int lda );
void bl1_dsyr2_blas( uplo1_t uplo, int m, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int lda );
void bl1_csyr2_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bl1_zsyr2_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );

// --- trmv ---

void bl1_strmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float*    a, int a_rs, int a_cs, float*    x, int incx );
void bl1_dtrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double*   a, int a_rs, int a_cs, double*   x, int incx );
void bl1_ctrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx );
void bl1_ztrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx );

void bl1_strmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float*    a, int lda, float*    x, int incx );
void bl1_dtrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double*   a, int lda, double*   x, int incx );
void bl1_ctrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* a, int lda, scomplex* x, int incx );
void bl1_ztrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* a, int lda, dcomplex* x, int incx );

// --- trsv ---

void bl1_strsv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float*    a, int a_rs, int a_cs, float*    x, int incx );
void bl1_dtrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double*   a, int a_rs, int a_cs, double*   x, int incx );
void bl1_ctrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx );
void bl1_ztrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx );

void bl1_strsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float*    a, int lda, float*    x, int incx );
void bl1_dtrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double*   a, int lda, double*   x, int incx );
void bl1_ctrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* a, int lda, scomplex* x, int incx );
void bl1_ztrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* a, int lda, dcomplex* x, int incx );

// --- trmvsx ---

void bl1_strmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy );
void bl1_dtrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy );
void bl1_ctrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_ztrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- trsvsx ---

void bl1_strsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy );
void bl1_dtrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy );
void bl1_ctrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_ztrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// end blis_prototypes_level2.h
// begin blis_prototypes_level3.h


// --- Level-3 BLAS-like prototypes --------------------------------------------

// --- gemm ---

void bl1_sgemm( trans1_t transa, trans1_t transb, int m, int k, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dgemm( trans1_t transa, trans1_t transb, int m, int k, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_cgemm( trans1_t transa, trans1_t transb, int m, int k, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_zgemm( trans1_t transa, trans1_t transb, int m, int k, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_sgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, float*    alpha, float*    a, int lda, float*    b, int ldb, float*    beta, float*    c, int ldc );
void bl1_dgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, double*   alpha, double*   a, int lda, double*   b, int ldb, double*   beta, double*   c, int ldc );
void bl1_cgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bl1_zgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- hemm ---

void bl1_shemm( side1_t side, uplo1_t uplo, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dhemm( side1_t side, uplo1_t uplo, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_chemm( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_zhemm( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_chemm_blas( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bl1_zhemm_blas( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- herk ---

void bl1_sherk( uplo1_t uplo, trans1_t trans, int m, int k, float*  alpha, float*    a, int a_rs, int a_cs, float*  beta, float*    c, int c_rs, int c_cs );
void bl1_dherk( uplo1_t uplo, trans1_t trans, int m, int k, double* alpha, double*   a, int a_rs, int a_cs, double* beta, double*   c, int c_rs, int c_cs );
void bl1_cherk( uplo1_t uplo, trans1_t trans, int m, int k, float*  alpha, scomplex* a, int a_rs, int a_cs, float*  beta, scomplex* c, int c_rs, int c_cs );
void bl1_zherk( uplo1_t uplo, trans1_t trans, int m, int k, double* alpha, dcomplex* a, int a_rs, int a_cs, double* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_cherk_blas( uplo1_t uplo, trans1_t trans, int m, int k, float*  alpha, scomplex* a, int lda, float*  beta, scomplex* c, int ldc );
void bl1_zherk_blas( uplo1_t uplo, trans1_t trans, int m, int k, double* alpha, dcomplex* a, int lda, double* beta, dcomplex* c, int ldc );

// --- her2k ---

void bl1_sher2k( uplo1_t uplo, trans1_t trans, int m, int k, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*  beta, float*    c, int c_rs, int c_cs );
void bl1_dher2k( uplo1_t uplo, trans1_t trans, int m, int k, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double* beta, double*   c, int c_rs, int c_cs );
void bl1_cher2k( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, float*  beta, scomplex* c, int c_rs, int c_cs );
void bl1_zher2k( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, double* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_cher2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, float*  beta, scomplex* c, int ldc );
void bl1_zher2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, double* beta, dcomplex* c, int ldc );

// --- symm ---

void bl1_ssymm( side1_t side, uplo1_t uplo, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dsymm( side1_t side, uplo1_t uplo, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_csymm( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_zsymm( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_ssymm_blas( side1_t side, uplo1_t uplo, int m, int n, float*    alpha, float*    a, int lda, float*    b, int ldb, float*    beta, float*    c, int ldc );
void bl1_dsymm_blas( side1_t side, uplo1_t uplo, int m, int n, double*   alpha, double*   a, int lda, double*   b, int ldb, double*   beta, double*   c, int ldc );
void bl1_csymm_blas( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bl1_zsymm_blas( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- syrk ---

void bl1_ssyrk( uplo1_t uplo, trans1_t trans, int m, int k, float*    alpha, float*    a, int a_rs, int a_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dsyrk( uplo1_t uplo, trans1_t trans, int m, int k, double*   alpha, double*   a, int a_rs, int a_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_csyrk( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_zsyrk( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_ssyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, float*    alpha, float*    a, int lda, float*    beta, float*    c, int ldc );
void bl1_dsyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, double*   alpha, double*   a, int lda, double*   beta, double*   c, int ldc );
void bl1_csyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int lda, scomplex* beta, scomplex* c, int ldc );
void bl1_zsyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* beta, dcomplex* c, int ldc );

// --- syr2k ---

void bl1_ssyr2k( uplo1_t uplo, trans1_t trans, int m, int k, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dsyr2k( uplo1_t uplo, trans1_t trans, int m, int k, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_csyr2k( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_zsyr2k( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_ssyr2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, float*    alpha, float*    a, int lda, float*    b, int ldb, float*    beta, float*    c, int ldc );
void bl1_dsyr2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, double*   alpha, double*   a, int lda, double*   b, int ldb, double*   beta, double*   c, int ldc );
void bl1_csyr2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bl1_zsyr2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- trmm ---

void bl1_strmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dtrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_ctrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_ztrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bl1_strmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int lda, float*    b, int ldb );
void bl1_dtrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int lda, double*   b, int ldb );
void bl1_ctrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb );
void bl1_ztrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb );

// --- trsm ---

void bl1_strsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dtrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_ctrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_ztrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bl1_strsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int lda, float*    b, int ldb );
void bl1_dtrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int lda, double*   b, int ldb );
void bl1_ctrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb );
void bl1_ztrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb );

// --- trmmsx ---

void bl1_strmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dtrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_ctrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_ztrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

// --- trsmsx ---

void bl1_strsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dtrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_ctrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_ztrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

// end blis_prototypes_level3.h

// begin blis_prototypes_fused1.h


// --- Fused Level-1 BLAS-like prototypes --------------------------------------

// --- axmyv2 ---

void bl1_saxmyv2( conj1_t conjx, int n, float*    alpha, float*    beta, float*    x, int inc_x, float*    y, int inc_y, float*    z, int inc_z );
void bl1_daxmyv2( conj1_t conjx, int n, double*   alpha, double*   beta, double*   x, int inc_x, double*   y, int inc_y, double*   z, int inc_z );
void bl1_caxmyv2( conj1_t conjx, int n, scomplex* alpha, scomplex* beta, scomplex* x, int inc_x, scomplex* y, int inc_y, scomplex* z, int inc_z );
void bl1_zaxmyv2( conj1_t conjx, int n, dcomplex* alpha, dcomplex* beta, dcomplex* x, int inc_x, dcomplex* y, int inc_y, dcomplex* z, int inc_z );

// --- axpyv2b ---

void bl1_saxpyv2b( int n, float*    beta1, float*    beta2, float*    a1, int inc_a1, float*    a2, int inc_a2, float*    w, int inc_w );
void bl1_daxpyv2b( int n, double*   beta1, double*   beta2, double*   a1, int inc_a1, double*   a2, int inc_a2, double*   w, int inc_w );
void bl1_caxpyv2b( int n, scomplex* beta1, scomplex* beta2, scomplex* a1, int inc_a1, scomplex* a2, int inc_a2, scomplex* w, int inc_w );
void bl1_zaxpyv2b( int n, dcomplex* beta1, dcomplex* beta2, dcomplex* a1, int inc_a1, dcomplex* a2, int inc_a2, dcomplex* w, int inc_w );

// --- axpyv3b ---

void bl1_saxpyv3b( int n, float*    beta1, float*    beta2, float*    beta3, float*    a1, int inc_a1, float*    a2, int inc_a2, float*    a3, int inc_a3, float*    w, int inc_w );
void bl1_daxpyv3b( int n, double*   beta1, double*   beta2, double*   beta3, double*   a1, int inc_a1, double*   a2, int inc_a2, double*   a3, int inc_a3, double*   w, int inc_w );
void bl1_caxpyv3b( int n, scomplex* beta1, scomplex* beta2, scomplex* beta3, scomplex* a1, int inc_a1, scomplex* a2, int inc_a2, scomplex* a3, int inc_a3, scomplex* w, int inc_w );
void bl1_zaxpyv3b( int n, dcomplex* beta1, dcomplex* beta2, dcomplex* beta3, dcomplex* a1, int inc_a1, dcomplex* a2, int inc_a2, dcomplex* a3, int inc_a3, dcomplex* w, int inc_w );

// --- axpyv2bdotaxpy ---

void bl1_saxpyv2bdotaxpy( int n, float*    beta, float*    u, int inc_u, float*    gamma, float*    z, int inc_z, float*    a, int inc_a, float*    x, int inc_x, float*    kappa, float*    rho, float*    w, int inc_w );
void bl1_daxpyv2bdotaxpy( int n, double*   beta, double*   u, int inc_u, double*   gamma, double*   z, int inc_z, double*   a, int inc_a, double*   x, int inc_x, double*   kappa, double*   rho, double*   w, int inc_w );
void bl1_caxpyv2bdotaxpy( int n, scomplex* beta, scomplex* u, int inc_u, scomplex* gamma, scomplex* z, int inc_z, scomplex* a, int inc_a, scomplex* x, int inc_x, scomplex* kappa, scomplex* rho, scomplex* w, int inc_w );
void bl1_zaxpyv2bdotaxpy( int n, dcomplex* beta, dcomplex* u, int inc_u, dcomplex* gamma, dcomplex* z, int inc_z, dcomplex* a, int inc_a, dcomplex* x, int inc_x, dcomplex* kappa, dcomplex* rho, dcomplex* w, int inc_w );

// --- dotsv2 ---

void bl1_sdotsv2( conj1_t conjxy, int n, float*    x, int inc_x, float*    y, int inc_y, float*    z, int inc_z, float*    beta, float*    rho_xz, float*    rho_yz );
void bl1_ddotsv2( conj1_t conjxy, int n, double*   x, int inc_x, double*   y, int inc_y, double*   z, int inc_z, double*   beta, double*   rho_xz, double*   rho_yz );
void bl1_cdotsv2( conj1_t conjxy, int n, scomplex* x, int inc_x, scomplex* y, int inc_y, scomplex* z, int inc_z, scomplex* beta, scomplex* rho_xz, scomplex* rho_yz );
void bl1_zdotsv2( conj1_t conjxy, int n, dcomplex* x, int inc_x, dcomplex* y, int inc_y, dcomplex* z, int inc_z, dcomplex* beta, dcomplex* rho_xz, dcomplex* rho_yz );

// --- dotsv3 ---

void bl1_sdotsv3( conj1_t conjxyw, int n, float*    x, int inc_x, float*    y, int inc_y, float*    w, int inc_w, float*    z, int inc_z, float*    beta, float*    rho_xz, float*    rho_yz, float*    rho_wz );
void bl1_ddotsv3( conj1_t conjxyw, int n, double*   x, int inc_x, double*   y, int inc_y, double*   w, int inc_w, double*   z, int inc_z, double*   beta, double*   rho_xz, double*   rho_yz, double*   rho_wz );
void bl1_cdotsv3( conj1_t conjxyw, int n, scomplex* x, int inc_x, scomplex* y, int inc_y, scomplex* w, int inc_w, scomplex* z, int inc_z, scomplex* beta, scomplex* rho_xz, scomplex* rho_yz, scomplex* rho_wz );
void bl1_zdotsv3( conj1_t conjxyw, int n, dcomplex* x, int inc_x, dcomplex* y, int inc_y, dcomplex* w, int inc_w, dcomplex* z, int inc_z, dcomplex* beta, dcomplex* rho_xz, dcomplex* rho_yz, dcomplex* rho_wz );

// --- dotaxpy ---

void bl1_sdotaxpy( int n, float*    a, int inc_a, float*    x, int inc_x, float*    kappa, float*    rho, float*    w, int inc_w );
void bl1_ddotaxpy( int n, double*   a, int inc_a, double*   x, int inc_x, double*   kappa, double*   rho, double*   w, int inc_w );
void bl1_cdotaxpy( int n, scomplex* a, int inc_a, scomplex* x, int inc_x, scomplex* kappa, scomplex* rho, scomplex* w, int inc_w );
void bl1_zdotaxpy( int n, dcomplex* a, int inc_a, dcomplex* x, int inc_x, dcomplex* kappa, dcomplex* rho, dcomplex* w, int inc_w );

// --- dotaxmyv2 ---

void bl1_sdotaxmyv2( int n, float*    alpha, float*    beta, float*    x, int inc_x, float*    u, int inc_u, float*    rho, float*    y, int inc_y, float*    z, int inc_z );
void bl1_ddotaxmyv2( int n, double*   alpha, double*   beta, double*   x, int inc_x, double*   u, int inc_u, double*   rho, double*   y, int inc_y, double*   z, int inc_z );
void bl1_cdotaxmyv2( int n, scomplex* alpha, scomplex* beta, scomplex* x, int inc_x, scomplex* u, int inc_u, scomplex* rho, scomplex* y, int inc_y, scomplex* z, int inc_z );
void bl1_zdotaxmyv2( int n, dcomplex* alpha, dcomplex* beta, dcomplex* x, int inc_x, dcomplex* u, int inc_u, dcomplex* rho, dcomplex* y, int inc_y, dcomplex* z, int inc_z );

// --- dotv2axpyv2b ---

void bl1_sdotv2axpyv2b( int n, float*    a1, int inc_a1, float*    a2, int inc_a2, float*    x,  int inc_x, float*    kappa1, float*    kappa2, float*    rho1, float*    rho2, float*    w, int inc_w );
void bl1_ddotv2axpyv2b( int n, double*   a1, int inc_a1, double*   a2, int inc_a2, double*   x,  int inc_x, double*   kappa1, double*   kappa2, double*   rho1, double*   rho2, double*   w, int inc_w );
void bl1_cdotv2axpyv2b( int n, scomplex* a1, int inc_a1, scomplex* a2, int inc_a2, scomplex* x,  int inc_x, scomplex* kappa1, scomplex* kappa2, scomplex* rho1, scomplex* rho2, scomplex* w, int inc_w );
void bl1_zdotv2axpyv2b( int n, dcomplex* a1, int inc_a1, dcomplex* a2, int inc_a2, dcomplex* x,  int inc_x, dcomplex* kappa1, dcomplex* kappa2, dcomplex* rho1, dcomplex* rho2, dcomplex* w, int inc_w );

// --- axpyv2bdots ---

void bl1_zaxpyv2bdots( int       n,
                       dcomplex* alpha1,
                       dcomplex* alpha2,
                       dcomplex* x1, int inc_x1,
                       dcomplex* x2, int inc_x2,
                       dcomplex* y,  int inc_y,
                       dcomplex* u,  int inc_u,
                       dcomplex* beta,
                       dcomplex* rho );
// end blis_prototypes_fused1.h

// begin blis_f77_name_mangling.h


// --- Define Fortran name-mangling macro --------------------------------------

// If the F77_FUNC name-mangling macro is undefined, then we we need to define
// it ourselves.
#ifndef F77_FUNC

  // Case 1: F77_FUNC is undefined because we're building for Windows.
  #ifdef BLIS1_ENABLE_WINDOWS_BUILD

    // Check whether we need to use uppercase Fortran routine names; otherwise
    // default to lowercase.
    #ifdef BLIS1_ENABLE_UPPERCASE_F77

      // Use uppercase routine names (no underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_upper
    #else

      // Use lowercase routine names (no underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_lower
    #endif

  // Case 2: F77_FUNC is undefined because we're in a Linux-like environment
  // that did not define it for us.
  #else

    // Check whether we need to use uppercase Fortran routine names; otherwise
    // default to lowercase.
    #ifdef BLIS1_ENABLE_UPPERCASE_F77

      // Use uppercase routine names (single underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_upper ## _
    #else

      // Use lowercase routine names (single underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_lower ## _
    #endif

  #endif // #ifdef BLIS1_ENABLE_WINDOWS_BUILD

#endif // #ifndef F77_FUNC

// end blis_f77_name_mangling.h

#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
// begin blis_prototypes_cblas.h


#include <stddef.h> // skipped


#define CBLAS_INDEX size_t  
enum CBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO      {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG      {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE      {CblasLeft=141, CblasRight=142};


float  cblas_sdsdot(const int N, const float alpha, const float *X,
                    const int incX, const float *Y, const int incY);
double cblas_dsdot(const int N, const float *X, const int incX, const float *Y,
                   const int incY);
float  cblas_sdot(const int N, const float  *X, const int incX,
                  const float  *Y, const int incY);
double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY);


void   cblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu);
void   cblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc);

void   cblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu);
void   cblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc);



float  cblas_snrm2(const int N, const float *X, const int incX);
float  cblas_sasum(const int N, const float *X, const int incX);

double cblas_dnrm2(const int N, const double *X, const int incX);
double cblas_dasum(const int N, const double *X, const int incX);

float  cblas_scnrm2(const int N, const void *X, const int incX);
float  cblas_scasum(const int N, const void *X, const int incX);

double cblas_dznrm2(const int N, const void *X, const int incX);
double cblas_dzasum(const int N, const void *X, const int incX);



CBLAS_INDEX cblas_isamax(const int N, const float  *X, const int incX);
CBLAS_INDEX cblas_idamax(const int N, const double *X, const int incX);
CBLAS_INDEX cblas_icamax(const int N, const void   *X, const int incX);
CBLAS_INDEX cblas_izamax(const int N, const void   *X, const int incX);




void cblas_sswap(const int N, float *X, const int incX, 
                 float *Y, const int incY);
void cblas_scopy(const int N, const float *X, const int incX, 
                 float *Y, const int incY);
void cblas_saxpy(const int N, const float alpha, const float *X,
                 const int incX, float *Y, const int incY);

void cblas_dswap(const int N, double *X, const int incX, 
                 double *Y, const int incY);
void cblas_dcopy(const int N, const double *X, const int incX, 
                 double *Y, const int incY);
void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);

void cblas_cswap(const int N, void *X, const int incX, 
                 void *Y, const int incY);
void cblas_ccopy(const int N, const void *X, const int incX, 
                 void *Y, const int incY);
void cblas_caxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);

void cblas_zswap(const int N, void *X, const int incX, 
                 void *Y, const int incY);
void cblas_zcopy(const int N, const void *X, const int incX, 
                 void *Y, const int incY);
void cblas_zaxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);



void cblas_srotg(float *a, float *b, float *c, float *s);
void cblas_srotmg(float *d1, float *d2, float *b1, const float b2, float *P);
void cblas_srot(const int N, float *X, const int incX,
                float *Y, const int incY, const float c, const float s);
void cblas_srotm(const int N, float *X, const int incX,
                float *Y, const int incY, const float *P);

void cblas_drotg(double *a, double *b, double *c, double *s);
void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double *P);
void cblas_drot(const int N, double *X, const int incX,
                double *Y, const int incY, const double c, const double s);
void cblas_drotm(const int N, double *X, const int incX,
                double *Y, const int incY, const double *P);



void cblas_sscal(const int N, const float alpha, float *X, const int incX);
void cblas_dscal(const int N, const double alpha, double *X, const int incX);
void cblas_cscal(const int N, const void *alpha, void *X, const int incX);
void cblas_zscal(const int N, const void *alpha, void *X, const int incX);
void cblas_csscal(const int N, const float alpha, void *X, const int incX);
void cblas_zdscal(const int N, const double alpha, void *X, const int incX);




void cblas_sgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY);
void cblas_sgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const float alpha,
                 const float *A, const int lda, const float *X,
                 const int incX, const float beta, float *Y, const int incY);
void cblas_strmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const float *A, const int lda, 
                 float *X, const int incX);
void cblas_stbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const float *A, const int lda, 
                 float *X, const int incX);
void cblas_stpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const float *Ap, float *X, const int incX);
void cblas_strsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const float *A, const int lda, float *X,
                 const int incX);
void cblas_stbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const float *A, const int lda,
                 float *X, const int incX);
void cblas_stpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const float *Ap, float *X, const int incX);

void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);
void cblas_dgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const double alpha,
                 const double *A, const int lda, const double *X,
                 const int incX, const double beta, double *Y, const int incY);
void cblas_dtrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const double *A, const int lda, 
                 double *X, const int incX);
void cblas_dtbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const double *A, const int lda, 
                 double *X, const int incX);
void cblas_dtpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const double *Ap, double *X, const int incX);
void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const double *A, const int lda, double *X,
                 const int incX);
void cblas_dtbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const double *A, const int lda,
                 double *X, const int incX);
void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const double *Ap, double *X, const int incX);

void cblas_cgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY);
void cblas_cgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const void *alpha,
                 const void *A, const int lda, const void *X,
                 const int incX, const void *beta, void *Y, const int incY);
void cblas_ctrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda, 
                 void *X, const int incX);
void cblas_ctbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda, 
                 void *X, const int incX);
void cblas_ctpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);
void cblas_ctrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda, void *X,
                 const int incX);
void cblas_ctbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX);
void cblas_ctpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);

void cblas_zgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY);
void cblas_zgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const void *alpha,
                 const void *A, const int lda, const void *X,
                 const int incX, const void *beta, void *Y, const int incY);
void cblas_ztrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda, 
                 void *X, const int incX);
void cblas_ztbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda, 
                 void *X, const int incX);
void cblas_ztpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);
void cblas_ztrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda, void *X,
                 const int incX);
void cblas_ztbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX);
void cblas_ztpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);



void cblas_ssymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const float alpha, const float *A,
                 const int lda, const float *X, const int incX,
                 const float beta, float *Y, const int incY);
void cblas_ssbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const int K, const float alpha, const float *A,
                 const int lda, const float *X, const int incX,
                 const float beta, float *Y, const int incY);
void cblas_sspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const float alpha, const float *Ap,
                 const float *X, const int incX,
                 const float beta, float *Y, const int incY);
void cblas_sger(const enum CBLAS_ORDER order, const int M, const int N,
                const float alpha, const float *X, const int incX,
                const float *Y, const int incY, float *A, const int lda);
void cblas_ssyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, float *A, const int lda);
void cblas_sspr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, float *Ap);
void cblas_ssyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, const float *Y, const int incY, float *A,
                const int lda);
void cblas_sspr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, const float *Y, const int incY, float *A);

void cblas_dsymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);
void cblas_dsbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const int K, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);
void cblas_dspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const double alpha, const double *Ap,
                 const double *X, const int incX,
                 const double beta, double *Y, const int incY);
void cblas_dger(const enum CBLAS_ORDER order, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda);
void cblas_dsyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *A, const int lda);
void cblas_dspr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *Ap);
void cblas_dsyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, const double *Y, const int incY, double *A,
                const int lda);
void cblas_dspr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, const double *Y, const int incY, double *A);



void cblas_chemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_chbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const int K, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_chpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *Ap,
                 const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_cgeru(const enum CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void cblas_cgerc(const enum CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void cblas_cher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha, const void *X, const int incX,
                void *A, const int lda);
void cblas_chpr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const float *alpha, const void *X,
                const int incX, void *A);
void cblas_cher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *A, const int lda);
void cblas_chpr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *Ap);

void cblas_zhemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_zhbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const int K, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_zhpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *Ap,
                 const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_zgeru(const enum CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void cblas_zgerc(const enum CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void cblas_zher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const void *X, const int incX,
                void *A, const int lda);
void cblas_zhpr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const double *alpha, const void *X,
                const int incX, void *A);
void cblas_zher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *A, const int lda);
void cblas_zhpr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *Ap);




void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc);
void cblas_ssymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *B, const int ldb, const float beta,
                 float *C, const int ldc);
void cblas_ssyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const float *A, const int lda,
                 const float beta, float *C, const int ldc);
void cblas_ssyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const float alpha, const float *A, const int lda,
                  const float *B, const int ldb, const float beta,
                  float *C, const int ldc);
void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb);
void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb);

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
void cblas_dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *B, const int ldb, const double beta,
                 double *C, const int ldc);
void cblas_dsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc);
void cblas_dsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const double alpha, const double *A, const int lda,
                  const double *B, const int ldb, const double beta,
                  double *C, const int ldc);
void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);
void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);

void cblas_cgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc);
void cblas_csymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_csyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const void *alpha, const void *A, const int lda,
                 const void *beta, void *C, const int ldc);
void cblas_csyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const void *beta,
                  void *C, const int ldc);
void cblas_ctrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);
void cblas_ctrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);

void cblas_zgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc);
void cblas_zsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_zsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const void *alpha, const void *A, const int lda,
                 const void *beta, void *C, const int ldc);
void cblas_zsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const void *beta,
                  void *C, const int ldc);
void cblas_ztrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);
void cblas_ztrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);



void cblas_chemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_cherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const void *A, const int lda,
                 const float beta, void *C, const int ldc);
void cblas_cher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const float beta,
                  void *C, const int ldc);

void cblas_zhemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_zherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const void *A, const int lda,
                 const double beta, void *C, const int ldc);
void cblas_zher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const double beta,
                  void *C, const int ldc);
// end blis_prototypes_cblas.h
#else
// begin blis_prototypes_blas.h


// --- Name-mangling macro definitions -----------------------------------------

// --- Name-mangle level-1 BLAS routines ---------------------------

#define F77_isamax F77_FUNC( isamax , ISAMAX )
#define F77_idamax F77_FUNC( idamax , IDAMAX )
#define F77_icamax F77_FUNC( icamax , ICAMAX )
#define F77_izamax F77_FUNC( izamax , IZAMAX )
#define F77_sasum  F77_FUNC( sasum  , SASUM  )
#define F77_dasum  F77_FUNC( dasum  , DASUM  )
#define F77_scasum F77_FUNC( scasum , SCASUM )
#define F77_dzasum F77_FUNC( dzasum , DZASUM )
#define F77_saxpy  F77_FUNC( saxpy  , SAXPY  )
#define F77_daxpy  F77_FUNC( daxpy  , DAXPY  )
#define F77_caxpy  F77_FUNC( caxpy  , CAXPY  )
#define F77_zaxpy  F77_FUNC( zaxpy  , ZAXPY  )
#define F77_scopy  F77_FUNC( scopy  , SCOPY  )
#define F77_dcopy  F77_FUNC( dcopy  , DCOPY  )
#define F77_ccopy  F77_FUNC( ccopy  , CCOPY  )
#define F77_zcopy  F77_FUNC( zcopy  , ZCOPY  )
#define F77_sdot   F77_FUNC( sdot   , SDOT   )
#define F77_ddot   F77_FUNC( ddot   , DDOT   )
#define F77_cdotu  F77_FUNC( cdotu  , CDOTU  )
#define F77_cdotc  F77_FUNC( cdotc  , CDOTC  )
#define F77_zdotu  F77_FUNC( zdotu  , ZDOTU  )
#define F77_zdotc  F77_FUNC( zdotc  , ZDOTC  )
#define F77_snrm2  F77_FUNC( snrm2  , SNRM2  )
#define F77_dnrm2  F77_FUNC( dnrm2  , DNRM2  )
#define F77_scnrm2 F77_FUNC( scnrm2 , SCNRM2 )
#define F77_dznrm2 F77_FUNC( dznrm2 , DZNRM2 )
#define F77_sscal  F77_FUNC( sscal  , SSCAL  )
#define F77_dscal  F77_FUNC( dscal  , DSCAL  )
#define F77_cscal  F77_FUNC( cscal  , CSCAL  )
#define F77_csscal F77_FUNC( csscal , CSSCAL )
#define F77_zscal  F77_FUNC( zscal  , ZSCAL  )
#define F77_zdscal F77_FUNC( zdscal , ZDSCAL )
#define F77_sswap  F77_FUNC( sswap  , SSWAP  )
#define F77_dswap  F77_FUNC( dswap  , DSWAP  )
#define F77_cswap  F77_FUNC( cswap  , CSWAP  )
#define F77_zswap  F77_FUNC( zswap  , ZSWAP  )

// --- Name-mangle level-2 BLAS routines ---------------------------

#define F77_sgemv  F77_FUNC( sgemv  , SGEMV  )
#define F77_dgemv  F77_FUNC( dgemv  , DGEMV  )
#define F77_cgemv  F77_FUNC( cgemv  , CGEMV  )
#define F77_zgemv  F77_FUNC( zgemv  , ZGEMV  )
#define F77_sger   F77_FUNC( sger   , SGER   )
#define F77_dger   F77_FUNC( dger   , DGER   )
#define F77_cgerc  F77_FUNC( cgerc  , CGERC  )
#define F77_cgeru  F77_FUNC( cgeru  , CGERU  )
#define F77_zgerc  F77_FUNC( zgerc  , ZGERC  )
#define F77_zgeru  F77_FUNC( zgeru  , ZGERU  )
#define F77_chemv  F77_FUNC( chemv  , CHEMV  )
#define F77_zhemv  F77_FUNC( zhemv  , ZHEMV  )
#define F77_cher   F77_FUNC( cher   , CHER   )
#define F77_zher   F77_FUNC( zher   , ZHER   )
#define F77_cher2  F77_FUNC( cher2  , CHER2  )
#define F77_zher2  F77_FUNC( zher2  , ZHER2  )
#define F77_ssymv  F77_FUNC( ssymv  , SSYMV  )
#define F77_dsymv  F77_FUNC( dsymv  , DSYMV  )
#define F77_ssyr   F77_FUNC( ssyr   , SSYR   )
#define F77_dsyr   F77_FUNC( dsyr   , DSYR   )
#define F77_ssyr2  F77_FUNC( ssyr2  , SSYR2  )
#define F77_dsyr2  F77_FUNC( dsyr2  , DSYR2  )
#define F77_strmv  F77_FUNC( strmv  , STRMV  )
#define F77_dtrmv  F77_FUNC( dtrmv  , DTRMV  )
#define F77_ctrmv  F77_FUNC( ctrmv  , CTRMV  )
#define F77_ztrmv  F77_FUNC( ztrmv  , ZTRMV  )
#define F77_strsv  F77_FUNC( strsv  , STRSV  )
#define F77_dtrsv  F77_FUNC( dtrsv  , DTRSV  )
#define F77_ctrsv  F77_FUNC( ctrsv  , CTRSV  )
#define F77_ztrsv  F77_FUNC( ztrsv  , ZTRSV  )

// --- Name-mangle level-3 BLAS routines ---------------------------

#define F77_sgemm  F77_FUNC( sgemm  , SGEMM  )
#define F77_dgemm  F77_FUNC( dgemm  , DGEMM  )
#define F77_cgemm  F77_FUNC( cgemm  , CGEMM  )
#define F77_zgemm  F77_FUNC( zgemm  , ZGEMM  )
#define F77_chemm  F77_FUNC( chemm  , CHEMM  )
#define F77_zhemm  F77_FUNC( zhemm  , ZHEMM  )
#define F77_cherk  F77_FUNC( cherk  , CHERK  )
#define F77_zherk  F77_FUNC( zherk  , ZHERK  )
#define F77_cher2k F77_FUNC( cher2k , CHER2K )
#define F77_zher2k F77_FUNC( zher2k , ZHER2K )
#define F77_ssymm  F77_FUNC( ssymm  , SSYMM  )
#define F77_dsymm  F77_FUNC( dsymm  , DSYMM  )
#define F77_csymm  F77_FUNC( csymm  , CSYMM  )
#define F77_zsymm  F77_FUNC( zsymm  , ZSYMM  )
#define F77_ssyrk  F77_FUNC( ssyrk  , SSYRK  )
#define F77_dsyrk  F77_FUNC( dsyrk  , DSYRK  )
#define F77_csyrk  F77_FUNC( csyrk  , CSYRK  )
#define F77_zsyrk  F77_FUNC( zsyrk  , ZSYRK  )
#define F77_ssyr2k F77_FUNC( ssyr2k , SSYR2K )
#define F77_dsyr2k F77_FUNC( dsyr2k , DSYR2K )
#define F77_csyr2k F77_FUNC( csyr2k , CSYR2K )
#define F77_zsyr2k F77_FUNC( zsyr2k , ZSYR2K )
#define F77_strmm  F77_FUNC( strmm  , STRMM  )
#define F77_dtrmm  F77_FUNC( dtrmm  , DTRMM  )
#define F77_ctrmm  F77_FUNC( ctrmm  , CTRMM  )
#define F77_ztrmm  F77_FUNC( ztrmm  , ZTRMM  )
#define F77_strsm  F77_FUNC( strsm  , STRSM  )
#define F77_dtrsm  F77_FUNC( dtrsm  , DTRSM  )
#define F77_ctrsm  F77_FUNC( ctrsm  , CTRSM  )
#define F77_ztrsm  F77_FUNC( ztrsm  , ZTRSM  )


// --- Prototypes --------------------------------------------------------------

// --- Level-1 BLAS prototypes -------------------

// --- amax ---
int      F77_isamax ( int* n, float*    x, int* incx );
int      F77_idamax ( int* n, double*   x, int* incx );
int      F77_icamax ( int* n, scomplex* x, int* incx );
int      F77_izamax ( int* n, dcomplex* x, int* incx );
// --- asum ---
float    F77_sasum  ( int* n, float*    x, int* incx );
double   F77_dasum  ( int* n, double*   x, int* incx );
float    F77_scasum ( int* n, scomplex* x, int* incx );
double   F77_dzasum ( int* n, dcomplex* x, int* incx );
// --- axpy ---
void     F77_saxpy  ( int* n, float*    alpha, float*    x, int* incx,  float*    y, int* incy );
void     F77_daxpy  ( int* n, double*   alpha, double*   x, int* incx,  double*   y, int* incy );
void     F77_caxpy  ( int* n, scomplex* alpha, scomplex* x, int* incx,  scomplex* y, int* incy );
void     F77_zaxpy  ( int* n, dcomplex* alpha, dcomplex* x, int* incx,  dcomplex* y, int* incy );
// --- copy ---
void     F77_scopy  ( int* n, float*    x, int* incx, float*    y, int* incy );
void     F77_dcopy  ( int* n, double*   x, int* incx, double*   y, int* incy );
void     F77_ccopy  ( int* n, scomplex* x, int* incx, scomplex* y, int* incy );
void     F77_zcopy  ( int* n, dcomplex* x, int* incx, dcomplex* y, int* incy );
// --- dot ---
float    F77_sdot   ( int* n, float*    x, int* incx, float*    y, int* incy );
double   F77_ddot   ( int* n, double*   x, int* incx, double*   y, int* incy );
scomplex F77_cdotu  ( int* n, scomplex* x, int* incx, scomplex* y, int* incy );
scomplex F77_cdotc  ( int* n, scomplex* x, int* incx, scomplex* y, int* incy );
dcomplex F77_zdotu  ( int* n, dcomplex* x, int* incx, dcomplex* y, int* incy );
dcomplex F77_zdotc  ( int* n, dcomplex* x, int* incx, dcomplex* y, int* incy );
// --- nrm2 ---
float    F77_snrm2  ( int* n, float*    x, int* incx );
double   F77_dnrm2  ( int* n, double*   x, int* incx );
float    F77_scnrm2 ( int* n, scomplex* x, int* incx );
double   F77_dznrm2 ( int* n, dcomplex* x, int* incx );
// --- scal ---
void     F77_sscal  ( int* n, float*    alpha, float*    y, int* incy );
void     F77_dscal  ( int* n, double*   alpha, double*   y, int* incy );
void     F77_cscal  ( int* n, scomplex* alpha, scomplex* y, int* incy );
void     F77_csscal ( int* n, float*    alpha, scomplex* y, int* incy );
void     F77_zscal  ( int* n, dcomplex* alpha, dcomplex* y, int* incy );
void     F77_zdscal ( int* n, double*   alpha, dcomplex* y, int* incy );
// --- swap ---
void     F77_sswap  ( int* n, float*    x, int* incx, float*    y, int* incy );
void     F77_dswap  ( int* n, double*   x, int* incx, double*   y, int* incy );
void     F77_cswap  ( int* n, scomplex* x, int* incx, scomplex* y, int* incy );
void     F77_zswap  ( int* n, dcomplex* x, int* incx, dcomplex* y, int* incy );

// --- Level-2 BLAS prototypes -------------------

// --- gemv ---
void     F77_sgemv  ( char* transa, int* m, int* n, float*    alpha, float*    a, int* lda, float*    x, int* incx, float*    beta, float*    y, int* incy );
void     F77_dgemv  ( char* transa, int* m, int* n, double*   alpha, double*   a, int* lda, double*   x, int* incx, double*   beta, double*   y, int* incy );
void     F77_cgemv  ( char* transa, int* m, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* x, int* incx, scomplex* beta, scomplex* y, int* incy );
void     F77_zgemv  ( char* transa, int* m, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* x, int* incx, dcomplex* beta, dcomplex* y, int* incy );
// --- ger ---
void     F77_sger   ( int* m, int* n, float*    alpha, float*    x, int* incx, float*    y, int* incy, float*    a, int* lda );
void     F77_dger   ( int* m, int* n, double*   alpha, double*   x, int* incx, double*   y, int* incy, double*   a, int* lda );
void     F77_cgerc  ( int* m, int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* y, int* incy, scomplex* a, int* lda );
void     F77_cgeru  ( int* m, int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* y, int* incy, scomplex* a, int* lda );
void     F77_zgerc  ( int* m, int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* y, int* incy, dcomplex* a, int* lda );
void     F77_zgeru  ( int* m, int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* y, int* incy, dcomplex* a, int* lda );
// --- hemv ---
void     F77_chemv  ( char* uplo, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* x, int* incx, scomplex* beta, scomplex* y, int* incy );
void     F77_zhemv  ( char* uplo, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* x, int* incx, dcomplex* beta, dcomplex* y, int* incy );
// --- her ---
void     F77_cher   ( char* uplo, int* n, float*    alpha, scomplex* x, int* incx, scomplex* a, int* lda );
void     F77_zher   ( char* uplo, int* n, double*   alpha, dcomplex* x, int* incx, dcomplex* a, int* lda );
// --- her2 ---
void     F77_cher2  ( char* uplo, int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* y, int* incy, scomplex* a, int* lda );
void     F77_zher2  ( char* uplo, int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* y, int* incy, dcomplex* a, int* lda );
// --- symv ---
void     F77_ssymv  ( char* uplo, int* n, float*    alpha, float*    a, int* lda, float*    x, int* incx, float*    beta, float*    y, int* incy );
void     F77_dsymv  ( char* uplo, int* n, double*   alpha, double*   a, int* lda, double*   x, int* incx, double*   beta, double*   y, int* incy );
// --- syr ---
void     F77_ssyr   ( char* uplo, int* n, float*    alpha, float*    x, int* incx, float*    a, int* lda );
void     F77_dsyr   ( char* uplo, int* n, double*   alpha, double*   x, int* incx, double*   a, int* lda );
// --- syr2 ---
void     F77_ssyr2  ( char* uplo, int* n, float*    alpha, float*    x, int* incx, float*    y, int* incy, float*    a, int* lda );
void     F77_dsyr2  ( char* uplo, int* n, double*   alpha, double*   x, int* incx, double*   y, int* incy, double*   a, int* lda );
// --- trmv ---
void     F77_strmv  ( char* uplo, char* transa, char* diag, int* n,  float*    a, int* lda, float*    y, int* incy );
void     F77_dtrmv  ( char* uplo, char* transa, char* diag, int* n,  double*   a, int* lda, double*   y, int* incy );
void     F77_ctrmv  ( char* uplo, char* transa, char* diag, int* n,  scomplex* a, int* lda, scomplex* y, int* incy );
void     F77_ztrmv  ( char* uplo, char* transa, char* diag, int* n,  dcomplex* a, int* lda, dcomplex* y, int* incy );
// --- trsv ---
void     F77_strsv  ( char* uplo, char* transa, char* diag, int* n,  float*    a, int* lda, float*    y, int* incy );
void     F77_dtrsv  ( char* uplo, char* transa, char* diag, int* n,  double*   a, int* lda, double*   y, int* incy );
void     F77_ctrsv  ( char* uplo, char* transa, char* diag, int* n,  scomplex* a, int* lda, scomplex* y, int* incy );
void     F77_ztrsv  ( char* uplo, char* transa, char* diag, int* n,  dcomplex* a, int* lda, dcomplex* y, int* incy );

// --- Level-3 BLAS prototypes -------------------

// --- gemm ---
void     F77_sgemm  ( char* transa, char* transb, int* m, int* n, int* k, float*    alpha, float*    a, int* lda, float*    b, int* ldb, float*    beta, float*    c, int* ldc );
void     F77_dgemm  ( char* transa, char* transb, int* m, int* n, int* k, double*   alpha, double*   a, int* lda, double*   b, int* ldb, double*   beta, double*   c, int* ldc );
void     F77_cgemm  ( char* transa, char* transb, int* m, int* n, int* k, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* beta, scomplex* c, int* ldc );
void     F77_zgemm  ( char* transa, char* transb, int* m, int* n, int* k, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* beta, dcomplex* c, int* ldc );
// --- hemm ---
void     F77_chemm  ( char* side, char* uplo, int* m, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* beta, scomplex* c, int* ldc );
void     F77_zhemm  ( char* side, char* uplo, int* m, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* beta, dcomplex* c, int* ldc );
// --- herk ---
void     F77_cherk  ( char* uplo, char* transa, int* n, int* k, float*  alpha, scomplex* a, int* lda, float*  beta, scomplex* c, int* ldc );
void     F77_zherk  ( char* uplo, char* transa, int* n, int* k, double* alpha, dcomplex* a, int* lda, double* beta, dcomplex* c, int* ldc );
// --- her2k ---
void     F77_cher2k ( char* uplo, char* transa, int* n, int* k, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb, float*  beta, scomplex* c, int* ldc );
void     F77_zher2k ( char* uplo, char* transa, int* n, int* k, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb, double* beta, dcomplex* c, int* ldc );
// --- symm ---
void     F77_ssymm  ( char* side, char* uplo, int* m, int* n, float*    alpha, float*    a, int* lda, float*    b, int* ldb, float*    beta, float*    c, int* ldc );
void     F77_dsymm  ( char* side, char* uplo, int* m, int* n, double*   alpha, double*   a, int* lda, double*   b, int* ldb, double*   beta, double*   c, int* ldc );
void     F77_csymm  ( char* side, char* uplo, int* m, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* beta, scomplex* c, int* ldc );
void     F77_zsymm  ( char* side, char* uplo, int* m, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* beta, dcomplex* c, int* ldc );
// --- syrk ---
void     F77_ssyrk  ( char* uplo, char* transa, int* n, int* k, float*    alpha, float*    a, int* lda, float*    beta, float*    c, int* ldc );
void     F77_dsyrk  ( char* uplo, char* transa, int* n, int* k, double*   alpha, double*   a, int* lda, double*   beta, double*   c, int* ldc );
void     F77_csyrk  ( char* uplo, char* transa, int* n, int* k, scomplex* alpha, scomplex* a, int* lda, scomplex* beta, scomplex* c, int* ldc );
void     F77_zsyrk  ( char* uplo, char* transa, int* n, int* k, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* beta, dcomplex* c, int* ldc );
// --- syr2k ---
void     F77_ssyr2k ( char* uplo, char* transa, int* n, int* k, float*    alpha, float*    a, int* lda, float*    b, int* ldb, float*    beta, float*    c, int* ldc );
void     F77_dsyr2k ( char* uplo, char* transa, int* n, int* k, double*   alpha, double*   a, int* lda, double*   b, int* ldb, double*   beta, double*   c, int* ldc );
void     F77_csyr2k ( char* uplo, char* transa, int* n, int* k, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* beta, scomplex* c, int* ldc );
void     F77_zsyr2k ( char* uplo, char* transa, int* n, int* k, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* beta, dcomplex* c, int* ldc );
// --- trmm ---
void     F77_strmm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, float*    alpha, float*    a, int* lda, float*    b, int* ldb );
void     F77_dtrmm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, double*   alpha, double*   a, int* lda, double*   b, int* ldb );
void     F77_ctrmm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb );
void     F77_ztrmm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb );
// --- trsm ---
void     F77_strsm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, float*    alpha, float*    a, int* lda, float*    b, int* ldb );
void     F77_dtrsm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, double*   alpha, double*   a, int* lda, double*   b, int* ldb );
void     F77_ctrsm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb );
void     F77_ztrsm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb );

// end blis_prototypes_blas.h
#endif

// End extern "C" construct block.
#ifdef __cplusplus
}
#endif

#endif
