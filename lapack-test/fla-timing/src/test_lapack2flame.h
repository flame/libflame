
#ifndef TEST_LAPACK2FLAME_H
#define TEST_LAPACK2FLAME_H

#include "FLAME.h"
#include "ctype.h"
#include "flops.h"
#include "float.h"

#define COMMENT_CHAR             '#'
#define SINGLE_BAR               "----------------------------------\n"
#define DOUBLE_BAR               "==================================\n"
#define MAX_NAME_LENGTH          256
#define MAX_STRING_LENGTH        4096
#define MAX_NUM_TESTITEMS        4096
#define MAX_NUM_PARAMS           4
#define MAX_TIME_VALUE           1e9
#define FLOPS_PER_UNIT_PERF      1e9
#define SECONDS_TO_SLEEP         5
#define NOT_DEFINED              0
#define PERCENT_UNIT             10
#define EPSILON_RESCALE          1000

// Test input/output
typedef struct
{
    char          testname[ MAX_NAME_LENGTH ];
    FLA_Datatype  datatype;
    FLA_Trans     trans[ MAX_NUM_PARAMS ];
    FLA_Uplo      uplos[ MAX_NUM_PARAMS ];
    dim_t         dims [ MAX_NUM_PARAMS ];
    dim_t         repeat;
} param_t;

typedef struct
{
    param_t       param;
    double        performance;
    double        residual;
} result_t;

typedef struct
{
    dim_t         size;
    result_t*     results;
} testitem_t;

// Utility macros
#define FP_PER_MUL(is_complex) ((is_complex) ?  (6.0) : (1.0))
#define FP_PER_ADD(is_complex) ((is_complex) ?  (2.0) : (1.0))

#define CHKERR(ierr)                                            \
  if ((ierr) != 0) {                                            \
    fprintf(stderr, "Error in %s : %d\n", __FILE__,__LINE__);   \
    return ierr;                                                \
  }

#define PROGRESS(name, i, size)                                         \
  if (size/PERCENT_UNIT) {                                              \
    if ((i+1) == size)                                                  \
      fprintf(stdout, "%-20s: [ %3.0f %%]      \n\n", name, ((double)(i+1)/size*100.0)); \
    else if (! ((i+1)%(size/PERCENT_UNIT)))                             \
      fprintf(stdout, "%-20s: [ %3.0f %%]      \n", name, ((double)(i+1)/size*100.0)); \
  } else {                                                              \
    if ((i+1) == size)                                                  \
      fprintf(stdout, "%-20s: [ %3.0f %%]      \n\n", name, ((double)(i+1)/size*100.0)); \
    else if ( (i+1) < size )                                            \
      fprintf(stdout, "%-20s: [ %3.0f %%]      \n", name, ((double)(i+1)/size*100.0)); \
  }

// FILE utility
int read_to_endl(FILE* in, char* str);
int move_to_endl(FILE* in);
int skip_comment(FILE* in);

// Test utility
int test_initialize(FLA_Datatype datatype);
int test_finalize();
int test_info( FILE* stream, const char* format, ...);

int flush_test_parameter( param_t* param );
int flush_test_result( result_t* result );

int get_io_filenames_from_list( FILE* fp, char* testkind, char* file_in, char* file_out );
int set_test_parameter( char* testname, FLA_Datatype datatype, dim_t repeat,
                        dim_t m, dim_t n, dim_t k,
                        FLA_Trans transa, FLA_Trans transb, FLA_Trans transc, FLA_Trans transd,
                        FLA_Uplo  uploa,  FLA_Uplo  uplob,  FLA_Uplo  uploc,  FLA_Uplo  uplod,
                        param_t *param );

int run( int (*func)(FILE*,param_t,result_t*),
         FILE* stream, param_t param, result_t* result );

int show_test_param( FILE* stream, param_t param );
int show_test_result( FILE* stream, result_t result );

int create_testitems( const char* name, testitem_t* testitem );
int free_testitems( testitem_t* testitem );

int write_test_result_to_file( testitem_t items, char* name );

// String utility
char* datatype2str( FLA_Datatype datatype );
FLA_Datatype str2datatype( char* name );

// My tests
int test_gemm  ( FILE* stream, param_t param, result_t *result);
int test_chol  ( FILE* stream, param_t param, result_t *result);
int test_trinv ( FILE* stream, param_t param, result_t *result);
int test_lu_piv( FILE* stream, param_t param, result_t *result);
int test_qr    ( FILE* stream, param_t param, result_t *result);
int test_lq    ( FILE* stream, param_t param, result_t *result);
int test_svd   ( FILE* stream, param_t param, result_t *result);
int test_appq  ( FILE* stream, param_t param, result_t *result);

#endif

