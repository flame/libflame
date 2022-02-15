/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <math.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>


#define PARAMETERS_FILENAME      "input.global.general"
#define OPERATIONS_FILENAME      "input.global.operations"
#define COMMENT_CHAR             '#'
#define MAX_BINARY_NAME_LENGTH   256
#define INPUT_BUFFER_SIZE        256
#define MAX_DT_STRING_LENGTH     32
#define MAX_STOR_STRING_LENGTH   32
#define MAX_FUNC_STRING_LENGTH   32
#define MAX_NUM_STORAGE          4
#define MAX_NUM_DATATYPES        4
#define FLOPS_PER_UNIT_PERF      1e9

#define DISABLE_ALL              0
#define SPECIFY                  1
#define DISABLE                  0
#define ENABLE                   1

#define MAX_PASS_STRING_LENGTH   32
#define TEST_FAIL_STRING     "FAILURE"
#define TEST_WARN_STRING     "MARGINAL PASS"
#define TEST_PASS_STRING     "PASS"

#define TEST_LAPACK    1

#define NUM_STORAGE_CHARS        3
#define STORAGE_SCHEME_CHARS     "crg"

#define ON_FAILURE_IGNORE_CHAR   'i'
#define ON_FAILURE_SLEEP_CHAR    's'
#define ON_FAILURE_ABORT_CHAR    'a'

#define SECONDS_TO_SLEEP         5


#if defined(_WIN32) && defined(ENABLE_ILP64) 
#define FS "ld"
#else
#define FS "d"
#endif

// Datatype
#define FLOAT             100
#define DOUBLE            101
#define COMPLEX           102
#define DOUBLE_COMPLEX    103
#define INT               104
#define CONSTANT          105

#if defined(ENABLE_ILP64)
typedef int64_t integer;
typedef uint64_t uinteger;
#else
typedef int integer;
typedef unsigned long uinteger;
#endif

typedef struct scomplex
{
	float real, imag;
} scomplex;

typedef struct dcomplex
{
	double real, imag;
} dcomplex;

typedef struct
{
	integer       n_repeats;
	integer       n_datatypes;
	char          datatype_char[ MAX_NUM_DATATYPES + 1 ];
	integer       datatype[ MAX_NUM_DATATYPES + 1 ];
	integer       p_first;
	integer       p_max;
	integer       p_inc;
	integer       p_nfact;
} test_params_t;

typedef struct
{
	int gesvd;
	int geqp3;
	int gerqf;
	int gerq2;
} test_ops_t;


typedef struct
{
	double failwarn_s;
	double warnpass_s;
	double failwarn_d;
	double warnpass_d;
	double failwarn_c;
	double warnpass_c;
	double failwarn_z;
	double warnpass_z;
} test_thresh_t;

typedef struct
{
	char *ops;
	void (*fp)(test_params_t);
}OPERATIONS;

// external declaration
extern int ilaver_(integer *vers_major__, integer *vers_minor__, integer *vers_patch__);

// Prototypes.
char* fla_test_get_string_for_result( double residual, integer datatype, test_thresh_t* thresh );
void fla_test_init_strings( void );
void fla_test_parse_command_line( int argc, char** argv );
void fla_test_output_op_struct( char* op_str, int op );
void fla_test_output_info( char* message, ... );
void fla_test_output_error( char* message, ... );
void fla_test_parse_message( FILE* output_stream, char* message, va_list args );
void fla_test_read_next_line( char* buffer, FILE* input_stream );
void fla_test_read_parameter_file( char* input_filename, test_params_t* params );
int fla_test_read_tests_for_op( FILE* input_stream, int* op, char* buffer );
void fla_test_lapack_suite( char* input_filename, test_params_t params, test_ops_t ops );

void fla_test_op_driver( char*            func_str,
							char*         impl_var_str,
							integer       n_pc,
							char**        pc_str,
							integer       n_matrices,
							test_params_t params,
							test_thresh_t thresh,
							void (*f_exp) (test_params_t, // params
										   integer,  // datatype
										   integer,  // p_cur
										   integer,  // pci (param combo counter)
										   integer,  // n_repeats
										   double*,       // perf
										   double*,	      // time
										   double* ) );   // residual

void fla_test_build_function_string( char*        func_base_str,
										char*        impl_var_str,
										unsigned int n_pc,
										char*        pc_str,
										char*        func_str );

void fill_string_with_n_spaces( char* str, unsigned int n_spaces );
void fla_test_sleep( void );
void fla_test_abort( void );
void fla_start_timer( void );
double fla_end_timer( void );