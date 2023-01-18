/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

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
#define FLA_TEST_FAIL_STRING     "FAILURE"
#define FLA_TEST_WARN_STRING     "MARGINAL PASS"
#define FLA_TEST_PASS_STRING     "PASS"

#define FLA_TEST_HIER_FRONT_END  (-1)
#define FLA_TEST_FLAT_FRONT_END  0
#define FLA_TEST_FLAT_UNB_VAR    1
#define FLA_TEST_FLAT_OPT_VAR    2
#define FLA_TEST_FLAT_BLK_VAR    3
#define FLA_TEST_FLAT_UNB_EXT    4
#define FLA_TEST_FLAT_BLK_EXT    5

#define NUM_STORAGE_CHARS        3
#define STORAGE_SCHEME_CHARS     "crg"

#define ON_FAILURE_IGNORE_CHAR   'i'
#define ON_FAILURE_SLEEP_CHAR    's'
#define ON_FAILURE_ABORT_CHAR    'a'

#define SECONDS_TO_SLEEP         5


typedef struct
{
	unsigned int  n_repeats;
	unsigned int  n_storage;
	char          storage[ MAX_NUM_STORAGE + 1 ];
	unsigned int  n_datatypes;
	char          datatype_char[ MAX_NUM_DATATYPES + 1 ];
	FLA_Datatype  datatype[ MAX_NUM_DATATYPES + 1 ];
	dim_t         b_flash;
	dim_t         b_alg_hier;
	dim_t         b_alg_flat;
	dim_t         p_first;
	dim_t         p_max;
	dim_t         p_inc;
	int           p_nfact;
	unsigned int  n_threads;
	char          reaction_to_failure;
} test_params_t;


typedef struct
{
	int flash_front;
	int fla_front;
	int fla_unb_vars;
	int fla_opt_vars;
	int fla_blk_vars;
	int fla_unb_ext;
	int fla_blk_ext;
} test_op_t;


typedef struct
{
	// BLAS level-3
	test_op_t gemm;
	test_op_t hemm;
	test_op_t herk;
	test_op_t her2k;
	test_op_t symm;
	test_op_t syrk;
	test_op_t syr2k;
	test_op_t trmm;
	test_op_t trsm;

	// LAPACK-level
	test_op_t chol;
	test_op_t lu_nopiv;
	test_op_t lu_piv;
	test_op_t lu_incpiv;
	test_op_t ldlt_nopiv_part;
	test_op_t lu_nopiv_i;
	test_op_t qrut;
	test_op_t qrutinc;
	test_op_t lqut;
	test_op_t apqut;
	test_op_t apqutinc;
	test_op_t caqrutinc;
	test_op_t apcaqutinc;
	test_op_t uddateut;
	test_op_t uddateutinc;
	test_op_t apqudut;
	test_op_t apqudutinc;
	test_op_t hessut;
	test_op_t tridiagut;
	test_op_t bidiagut;
	test_op_t eig_gest;
	test_op_t trinv;
	test_op_t spdinv;
	test_op_t sylv;
	test_op_t lyap;
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



// Prototypes.
char* libfla_test_get_string_for_datatype( FLA_Datatype datatype );
char* libfla_test_get_string_for_storage( char storage );
char* libfla_test_get_string_for_result( double         residual,
                                         FLA_Datatype   datatype,
                                         test_thresh_t* thresh );
void libfla_test_init_strings( void );

void libfla_test_fill_storage_strings( char** sc_str, unsigned int n_storage_runs,
                                                      unsigned int n_matrices );
void carryover( unsigned int* c, unsigned int n_matrices );

void libfla_test_parse_command_line( int argc, char** argv );

void libfla_test_output_op_struct( char* op_str, test_op_t op );
void libfla_test_output_op_struct_fla_only( char* op_str, test_op_t op );
void libfla_test_output_op_struct_flash_only( char* op_str, test_op_t op );
void libfla_test_output_op_struct_front_only( char* op_str, test_op_t op );
void libfla_test_output_op_struct_front_fla_only( char* op_str, test_op_t op );
void libfla_test_output_op_struct_blas3( char* op_str, test_op_t op );

void libfla_test_output_info( char* message, ... );
void libfla_test_output_error( char* message, ... );
void libfla_test_parse_message( FILE* output_stream, char* message, va_list args );

void libfla_test_read_next_line( char* buffer, FILE* input_stream );

void libfla_test_read_parameter_file( char* input_filename, test_params_t* params );

void libfla_test_read_operation_file( char* input_filename, test_ops_t* ops );
void libfla_test_read_tests_for_op( FILE* input_stream, test_op_t* op );
void libfla_test_read_tests_for_op_fla_only( FILE* input_stream, test_op_t* op );
void libfla_test_read_tests_for_op_flash_only( FILE* input_stream, test_op_t* op );
void libfla_test_read_tests_for_op_front_only( FILE* input_stream, test_op_t* op );
void libfla_test_read_tests_for_op_front_fla_only( FILE* input_stream, test_op_t* op );
void libfla_test_read_tests_for_op_blas3( FILE* input_stream, test_op_t* op );

void libfla_test_blas3_suite( FILE* output_stream, test_params_t params, test_ops_t ops );
void libfla_test_lapack_suite( FILE* output_stream, test_params_t params, test_ops_t ops );

void libfla_test_op_driver( char*         func_str,
                            char*         impl_var_str,
                            unsigned int  first_var,
                            unsigned int  last_var,
                            unsigned int  n_pc,
                            char**        pc_str,
                            unsigned int  n_matrices,
                            signed int    impl,
                            test_params_t params,
                            test_thresh_t thresh,
                            void (*f_exp) (test_params_t, // params
                                           unsigned int,  // var
                                           char*,         // sc_cur_str (current storage string)
                                           FLA_Datatype,  // datatype
                                           uinteger,  // p_cur
                                           unsigned int,  // pci (param combo counter)
                                           unsigned int,  // n_repeats
                                           signed int,    // impl
                                           double*,       // perf
                                           double*,	      // time
                                           double* ) );   // residual

void libfla_test_print_result_info( char  *func_param_str,
                                    char  *datatype_char,
                                    char  *sc_str,
                                    integer    p_cur,
                                    double perf,
				    double time_min,
                                    double residual,
                                    char  *pass_str,
                                    int    nfact );

void libfla_test_build_function_string( char*        func_base_str,
                                        signed int   impl,
                                        char*        impl_var_str,
                                        unsigned int var,
                                        unsigned int n_pc,
                                        char*        pc_str,
                                        char*        func_str );

void fill_string_with_n_spaces( char* str, unsigned int n_spaces );
void libfla_test_obj_create( FLA_Datatype dt, FLA_Trans trans, char storage, dim_t m, dim_t n, FLA_Obj* A );
void libfla_test_sleep( void );
void libfla_test_abort( void );
