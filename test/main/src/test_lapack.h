/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#ifndef TEST_LAPACK_H
#define TEST_LAPACK_H

#include <string.h>
#include <time.h>
#ifdef _WIN32
#include <Windows.h>
#endif
#include <stdarg.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#ifndef _WIN32
#include <unistd.h>
#include <sys/time.h>
#endif
#include "test_common.h"

#define OPERATIONS_FILENAME                "input.global.operations"
#define LINEAR_PARAMETERS_FILENAME         "Config/LIN_SLVR.dat"
#define SYM_EIG_PARAMETERS_FILENAME        "Config/EIG_PARAMS.dat"
#define SVD_PARAMETERS_FILENAME            "Config/SVD.dat"
#define NON_SYM_EIG_PARAMETERS_FILENAME    "Config/EIG_NSYM_PARAMS.dat"

#define COMMENT_CHAR             '#'
#define MAX_BINARY_NAME_LENGTH   256
#define INPUT_BUFFER_SIZE        256
#define MAX_DT_STRING_LENGTH     32
#define MAX_STOR_STRING_LENGTH   32
#define MAX_FUNC_STRING_LENGTH   20
#define MAX_NUM_STORAGE          4
#define MAX_NUM_DATATYPES        4
#define FLOPS_PER_UNIT_PERF      1e9

#define DISABLE_ALL              0
#define SPECIFY                  1
#define DISABLE                  0
#define ENABLE                   1

#define RECT_INPUT              0
#define SQUARE_INPUT            1

#define MAX_PASS_STRING_LENGTH   32

#define NUM_STORAGE_CHARS        3
#define STORAGE_SCHEME_CHARS     "crg"

#define NUM_SUB_TESTS (4)

#if defined(FLA_ENABLE_ILP64)
#ifdef _WIN32
#define FT_IS "lld"
#else
#define FT_IS "ld"
#endif
#else
#define FT_IS "d"
#endif

// API categories
#define LIN            (1)
#define EIG_SYM        (2)
#define EIG_NSYM       (3)
#define SVD            (4)

typedef struct Lin_solver_paramlist_t
{
    integer mode; // Any one of these:- 0: discrete, 1: Combinational, 2:Range with steps of increment of Matrix sizes
    integer num_ranges; // number of ranges to run
    integer m_range_start;
    integer m_range_end;
    integer m_range_step_size;
    integer n_range_start;
    integer n_range_end;
    integer n_range_step_size;
    integer num_tests;
    integer num_repeats;
    integer num_data_types;
    integer data_types[MAX_NUM_DATATYPES];
    char data_types_char[MAX_NUM_DATATYPES];
    integer matrix_layout;
    char Uplo;
    char transr; // Must be 'N' or 'T' or 'C'.
    integer m;
    integer n;  // The order of A; the number of rows in B
    integer nrhs; // number of rhight hand sides
    integer lda; //  leading dimension of the array a
    integer ldb; //  leading dimension of the array b
    integer ldab;  //  leading dimension of the array ab
    integer kl; // number of subdiagonals
    integer ku; // number of superdiagonals
    integer kd; // number of super or sub diagonals
    char diag; // flag to indicate unit diagonal

    // below params are used only by Lin solver driver APIs.
    char fact;  // Must be 'F', 'N', or 'E'.
    char equed; // Must be 'N', 'R'. 'C', 'B'
    char symm; // if symmetric 'S' or Hermitian 'H'
    float solver_threhold;// threshold to verify PASS/FAIL criteria
    char equed_porfsx; // Must be 'N', 'Y'.
    integer  n_err_bnds_porfsx;
    integer  nparams_porfsx;
    char norm_gbcon; // norm param for gbcon API
    integer kl_gbcon; // number of subdiagonals
    integer ku_gbcon; // number of superdiagonals
    integer ldab_gbcon;  //  leading dimension of the array ab
} Lin_solver_paramlist;

/* struct to hold eigen parameters */
typedef struct EIG_paramlist_t
{
    integer mode; // Any one of these:- 0: discrete, 1: Combinational, 2:Range with steps of increment of Matrix sizes
    integer num_ranges; // number of ranges to run
    integer m_range_start;
    integer m_range_end;
    integer m_range_step_size;
    integer n_range_start;
    integer n_range_end;
    integer n_range_step_size;
    integer num_tests;
    integer num_repeats;
    integer num_data_types;
    integer data_types[MAX_NUM_DATATYPES];
    char data_types_char[MAX_NUM_DATATYPES];
    integer matrix_layout;
    char trans; // Must be 'N' or 'T' or 'C'.
    char uplo; // Must be 'U' or 'L'
    char job; // Must be 'N', 'P', 'S' or 'B'
    char jobz; //Must be 'N' or 'V'
    char vect; // Vector must be 'Q' or  'P'
    integer m;   //
    integer n;  // The order of A; the number of rows in B
    integer p; //
    integer nrhs; // number of rhight hand sides
    integer lda; //  leading dimension of the array a
    integer ldb; //  leading dimension of the array b
    integer nb;  //  leading dimension of the array ab
    integer ldt; // number of subdiagonals
    integer k;
    integer isgn;
    char compz;
    integer kb;
    integer itype;
    char vect_rd;
    char side;
    char job_seqr; // Must be 'E', 'S'
    char eigsrc; // Must be 'Q' or  'N'.
    char initv; // Must be 'N' or 'U'.
    char norm;
    char diag;
    char storev;
    integer tsize;
    integer threshold_value; // threshold value for EIG
}EIG_paramlist;


/* struct to hold eigen parameters */
typedef struct EIG_Non_symmetric_paramlist_t
{
    integer mode; // Any one of these:- 0: discrete, 1: Combinational, 2:Range with steps of increment of Matrix sizes
    integer num_ranges; // number of ranges to run
    integer m_range_start;
    integer m_range_end;
    integer m_range_step_size;
    integer n_range_start;
    integer n_range_end;
    integer n_range_step_size;
    integer num_repeats;
    integer num_tests;
    integer num_data_types;
    integer data_types[MAX_NUM_DATATYPES];
    char data_types_char[MAX_NUM_DATATYPES];
    integer matrix_layout;
    integer n;  // The order of A; the number of rows in B
    char howmny; // Must be 'A' or 'B' or 'S'.
    char initv; // Must be 'N' or 'U'.
    char job_seqr; // Must be 'E', 'S'
    char eigsrc; // Must be 'Q' or  'N'.
    char side;  // Must be 'R' or 'L' or 'B'.
    char job; //  Must be 'E' or 'V' or 'B'.
    char howmny_trsna; // Must be 'A' or 'S'.

    /* used params for trsen API  */
    char job_trsen;// must bie 'N' or 'E' or 'V' or 'B'.
    char compq; // Must be 'V' or 'N' .

    /* used params for trsyl API  */
    char trana_real; //  Must be 'N' or 'T' or 'C'.
    char trana_complex; //  Must be 'N' or 'T' or 'C'.
    char tranb_real; //  Must be 'N' or 'T' or 'C'.
    char tranb_complex; //  Must be 'N' or 'T' or 'C'.
    integer isgn; //  +1 or -1

    /* Thresholds for the APIs  */
    float gghrd_threshold; // threshold for the gghrd API
    float ggbal_threshold; // threshold for the ggbal API
    float GenNonSymEigProblem_threshold; // threshold for the ggbal API

    char compq_hgeqz; // Must be 'I' or 'V' or 'N'
    char compz_hgeqz; // Must be 'I' or 'V' or 'N'

    char side_tgevc; // Must be 'R', 'L', or 'B'.
    char jobvsl; // Must be 'N', or 'V'.
    char jobvsr; // Must be 'N', or 'V'.
    char sort; // Must be 'N', or 'S'.
    char sense_ggesx;// must be 'N' or 'E' or 'V' or 'B'.
    char balance_ggevx;// must be 'N' or 'P' or 'S' or 'B'.
    char sense_ggevx;// must be 'N' or 'E' or 'V' or 'B'.
    char sort_gees; // Must be 'N', or 'S'.
    integer  wantz; // Must be 1 or 0
    integer  wantq; // Must be 1 or 0
    integer  tgsen_ijob; // Must be between 0 to 5
    char unmhr_trans; // Must be N or C
}EIG_Non_symmetric_paramlist;


/* struct to hold SVD parameters */
typedef struct SVD_paramlist_t
{
    integer mode; // Any one of these:- 0: discrete, 1: Combinational, 2:Range with steps of increment of Matrix sizes
    integer num_ranges; // number of ranges to run
    integer m_range_start;
    integer m_range_end;
    integer m_range_step_size;
    integer n_range_start;
    integer n_range_end;
    integer n_range_step_size;
    integer num_repeats;
    integer num_tests;
    integer num_data_types;
    integer data_types[MAX_NUM_DATATYPES];
    char data_types_char[MAX_NUM_DATATYPES];
    integer matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobu; // Must be 'U' or 'N'.
    char jobv; // Must be 'V' or 'N'.
    char jobq; // Must be 'Q' or 'N'.
    integer m; // The number of rows of the matrix A
    integer p; // The number of rows of the matrix B
    integer n; // The number of columns of the matrices A and B
    float tola;
    float tolb;
    char jobu_gesvd; // Must be 'A', 'S', 'O', or 'N'.
    char jobvt_gesvd; // Must be 'A', 'S', 'O', or 'N'.

    char joba_gejsv; //  Must be 'C', 'E', 'F', 'G', 'A', or 'R'.
    char jobu_gejsv; // Must be 'U', 'F', 'W', or 'N'.
    char jobv_gejsv; // Must be 'V', 'J', 'W', or 'N'.
    char jobr_gejsv; // Must be 'N' or 'R'.
    char jobt_gejsv; // Must be 'T' or 'N'.
    char jobp_gejsv; //  Must be 'P' or 'N'.
    integer m_gejsv; // The number of rows of the matrix A
    integer n_gejsv; // The number of rows of the matrix B

    /* Parameters for 'gesvj' API  */
    char joba_gesvj; //  Must be 'L', 'U' or 'G'.
    char jobu_gesvj; // Must be 'U', 'C' or 'N'.
    char jobv_gesvj; // Must be 'V', 'A' or 'N'.
    integer m_gesvj; // The number of rows of the matrix A
    integer n_gesvj; // The number of rows of the matrix B
    integer mv_gesvj;
    float ctol_gesvj; // convergence of threshold

    /* Parameters for 'gesvdx' API  */
    char jobu_gesvdx; //  Must be 'V', or 'N'.
    char jobvt_gesvdx; // Must be 'V', or 'N'.
    char range_gesvdx; // Must be 'A', 'V', 'I'.
    integer il, iu; // the indices of the smallest and largest singular values.
    float vl, vu; //  the lower and upper bounds of the interval.

    /* Parameters for 'gesvdq' API  */
        char joba_gesvdq; //  Must be 'A', 'H', 'M' , 'E'
        char jobu_gesvdq; // Must be 'A', 'S', 'R' , 'N'
        char jobv_gesvdq; // Must be 'A', 'V', 'R' , 'N'.

    /* Thresholds for the APIs  */
    float svd_threshold; // threshold for the gghrd API

}SVD_paramlist;


typedef struct
{
    integer       n_repeats;
    integer       n_datatypes;
    char          *datatype_char;
    integer       *datatype;
    integer       p_first;
    integer       p_max;
    integer       p_inc;
    integer       p_nfact;

    struct SVD_paramlist_t svd_paramslist[NUM_SUB_TESTS];
    struct EIG_Non_symmetric_paramlist_t eig_non_sym_paramslist[NUM_SUB_TESTS];
    struct EIG_paramlist_t eig_sym_paramslist[NUM_SUB_TESTS];
    struct Lin_solver_paramlist_t lin_solver_paramslist[NUM_SUB_TESTS];

} test_params_t;

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
    void (*fp)(test_params_t *);
}OPERATIONS;

// external declaration
extern integer ilaver_(integer *vers_major__, integer *vers_minor__, integer *vers_patch__);

// Prototypes.
char* fla_test_get_string_for_result( double residual, integer datatype, double thresh );
void fla_test_init_strings( void );
void fla_test_parse_command_line( integer argc, char** argv );
void fla_test_output_op_struct( char* op_str, integer op );
void fla_test_output_info( char* message, ... );
void fla_test_output_error( char* message, ... );
void fla_test_parse_message( FILE* output_stream, char* message, va_list args );
void fla_test_read_next_line( char* buffer, FILE* input_stream );
integer fla_test_read_tests_for_op( FILE* input_stream, integer* op, char* buffer );

/*Read Linear API parameters from config file */
void fla_test_read_linear_param( const char* input_filename, test_params_t* params );

/*Function to read EIG parametes from config file */
void fla_test_read_sym_eig_params( const char* input_filename, test_params_t* params );

/*Function to read EIG non-symmetric parametes from config file */
void fla_test_read_non_sym_eig_params( const char* input_filename, test_params_t* params );

/*Function to read SVD parametes from config file */
void fla_test_read_svd_params ( const char* input_filename, test_params_t* params );

void fla_test_lapack_suite( char* input_filename, test_params_t *params );

void fla_test_op_driver( char*            func_str,
                         integer       square_inp,
                         test_params_t *params,
                         integer       api_type,
                         void (*f_exp) (test_params_t *, // params
                                        integer,        // datatype
                                        integer,        // p_cur
                                        integer,        // q_cur
                                        integer,        // pci (param combo counter)
                                        integer,        // n_repeats
                                        double*,        // perf
                                        double*,        // time
                                        double* ) );    // residual

void fla_test_build_function_string( char*        func_base_str,
                                        char*        impl_var_str,
                                        char*        func_str );
void fill_string_with_n_spaces( char* str, integer n_spaces );
double fla_test_clock(void);
void fla_test_get_time_unit(char * scale , double * time);

#endif // TEST_LAPACK_H
