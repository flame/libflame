/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "ctype.h"
#include "test_lapack.h"
#include "test_routines.h"

// Global variables.
int n_threads = 1;
char fla_test_binary_name[ MAX_BINARY_NAME_LENGTH + 1 ];
char fla_test_pass_string[ MAX_PASS_STRING_LENGTH + 1 ];
char fla_test_warn_string[ MAX_PASS_STRING_LENGTH + 1 ];
char fla_test_fail_string[ MAX_PASS_STRING_LENGTH + 1 ];
char fla_test_incomplete_string [ MAX_PASS_STRING_LENGTH + 1];
char fla_test_invalid_string [ MAX_PASS_STRING_LENGTH + 1];
char fla_test_storage_format_string[ 200 ];
char fla_test_stor_chars[ NUM_STORAGE_CHARS + 1 ];
double ref_time_sec = 0.0;

integer total_tests;
integer total_failed_tests;
integer total_incomplete_tests;
integer tests_passed[4];
integer tests_failed[4];
integer tests_incomplete[4];

#define SKIP_EXTRA_LINE_READ() \
        eol = fgetc(fp); \
        if(eol != '\n') \
            fscanf(fp, "%*[^\n]\n")


#define CHECK_LINE_SKIP() \
        eol = fgetc(fp); \
        if((eol == '\r') || (eol == '\n')){ \
            num_ranges = ( (i+1) < num_ranges)? (i+1):num_ranges; \
            break; \
        } \

#if AOCL_FLA_SET_PROGRESS_ENABLE == 1
int aocl_fla_progress(const char* const api,const integer lenapi,const  integer* const progress,const integer* const current_thread,const integer* const total_threads)
{
  printf( "In AOCL FLA  Progress thread  %"FT_IS", at API  %s, progress  %"FT_IS" total threads= %"FT_IS"\n", *current_thread, api, *progress,*total_threads );
  return 0;

}
#endif

#if AOCL_FLA_SET_PROGRESS_ENABLE == 2
int test_progress(const char* const api,const integer lenapi,const integer* const progress,const integer* const current_thread,const integer* const total_threads)
{
  printf( "In AOCL Progress thread  %"FT_IS", at API  %s, progress %"FT_IS" total threads= %"FT_IS" \n", *current_thread, api, *progress,*total_threads );
  return 0;

}
#endif


int  main( int argc, char** argv )
{
    test_params_t params;
    integer vers_major, vers_minor, vers_patch;

    ilaver_(&vers_major, &vers_minor, &vers_patch);

    printf(" LAPACK version: %"FT_IS".%"FT_IS".%"FT_IS" \n", vers_major, vers_minor, vers_patch);

    /* Initialize some strings. */
    fla_test_init_strings();

    /* Check if multithread variable is set or not*/
    char* str = getenv("FLA_TEST_NUM_THREADS");

    if( str != NULL){

        n_threads = atoi(str);

        if(n_threads < 1){
            fla_test_output_error("Number of threads should be greater than or equal to 1 \n");
        }
    }

    /* Check for Command line requests */
    if ( argc > 1 )
    {
        fla_test_execute_cli_api(argc, argv, &params);
    }
    else
    {
        /* Copy the binary name to a global string so we can use it later. */
        strncpy( fla_test_binary_name, argv[0], MAX_BINARY_NAME_LENGTH );

        /* Read Linear API parameters from config file */
        fla_test_read_linear_param( LINEAR_PARAMETERS_FILENAME, &params );

        /* Read eigen parameters from config file */
        fla_test_read_sym_eig_params( SYM_EIG_PARAMETERS_FILENAME, &params );
        fla_test_read_non_sym_eig_params( NON_SYM_EIG_PARAMETERS_FILENAME, &params );

        /* Read SVD parameters from config file */
        fla_test_read_svd_params ( SVD_PARAMETERS_FILENAME, &params  );

        #if AOCL_FLA_SET_PROGRESS_ENABLE == 2
            aocl_fla_set_progress(test_progress);
        #endif

        /* Test the LAPACK-level operations. */
        fla_test_lapack_suite( OPERATIONS_FILENAME, &params );
    }

    return 0;
}


/* This function reads the operation file to execute selected LAPACK APIs*/
void fla_test_lapack_suite( char* input_filename, test_params_t *params )
{
    char buffer[ INPUT_BUFFER_SIZE ];
    integer check_flag;
    integer i, j, op, test_api_count, g_id, test_group_count;
    FILE* input_stream;

    total_tests = 0;
    total_failed_tests = 0;
    total_incomplete_tests = 0;
    tests_passed[0] = tests_passed[1] = tests_passed[2] = tests_passed[3] = 0;
    tests_failed[0] = tests_failed[1] = tests_failed[2] = tests_failed[3] = 0;
    tests_incomplete[0] = tests_incomplete[1] = tests_incomplete[2] = tests_incomplete[3] = 0;

    test_api_count = sizeof(API_test_functions) / sizeof(API_test_functions[0]);
    test_group_count = sizeof(API_test_group) / sizeof(API_test_group[0]);

    fla_test_output_info( "\n" );
    fla_test_output_info( "--- LAPACK-level operation tests ---------------------\n" );
    fla_test_output_info( "\n" );

    // Attempt to open input file corresponding to input_filename as
    // read-only/binary.
    input_stream = fopen( input_filename, "rb" );

    // Check for success.
    if ( input_stream == NULL )
    {
        fla_test_output_error( "Failed to open input file %s. Check existence and permissions.\n", input_filename );
    }

    // Check for '2' option in input config and 3 for group testing
    check_flag = fla_test_check_run_only(input_stream, &op,  buffer);

    if(check_flag == 3)
    {
        for(i =0; i < test_group_count; i++)
        {
            fla_test_read_tests_for_op( input_stream, &op, buffer);
            if(op == 1)
            {
                g_id = fla_test_get_group_id(buffer);
                if(g_id >= 0)
                {
                    for( j = 0; j < test_api_count; j++)
                    {
                        if(API_test_functions[j].type == g_id)
                        {
                            API_test_functions[j].fp(1, NULL, params);
                        }
                    }
                }
            }
        }
    }
    else
    {
        while(fla_test_read_tests_for_op(input_stream, &op, buffer))
        {
            if(op == check_flag)
            {
                // Check if the specified API is supported in test suite
                for( i=0; i < test_api_count; i++)
                {
                    if(!strcmp(API_test_functions[i].ops, buffer))
                    {
                        API_test_functions[i].fp(1, NULL, params);
                    }
                }

                op = 0;
            }
        }
    }
    fclose( input_stream );
    fla_test_print_summary();
}

void fla_test_print_summary()
{
    fla_test_output_info("\n\nResults Summary:\n\n");
    fla_test_output_info("%2sDATATYPE%13s No. of Tests%6s Passed%9s Failed%9s Incomplete\n", "", "", "", "", "");
    fla_test_output_info( "=====================================================================================\n" );
    fla_test_output_info("%2sFLOAT%15s %8d%8s %8d%6s %8d%15s %d\n", "", "", tests_passed[0] + tests_failed[0] + tests_incomplete[0], "", tests_passed[0], "", tests_failed[0], "", tests_incomplete[0]);
    fla_test_output_info("%2sDOUBLE%14s %8d%8s %8d%6s %8d%15s %d\n", "", "", tests_passed[1] + tests_failed[1]  + tests_incomplete[0], "", tests_passed[1], "", tests_failed[1], "", tests_incomplete[1]);
    fla_test_output_info("%2sCOMPLEX%13s %8d%8s %8d%6s %8d%15s %d\n", "", "", tests_passed[2] + tests_failed[2]  + tests_incomplete[0], "", tests_passed[2], "", tests_failed[2], "", tests_incomplete[2]);
    fla_test_output_info("%2sDOUBLE COMPLEX%6s %8d%8s %8d%6s %8d%15s %d\n", "", "", tests_passed[3] + tests_failed[3]  + tests_incomplete[0], "", tests_passed[3], "", tests_failed[3], "", tests_incomplete[3]);

    if(total_failed_tests > 0)
        printf("\n\nThere are failed tests, Please look at output log for more details\n\n");
    
    if(total_incomplete_tests > 0)
        printf("\n\nThere are some incomplete tests, Please look at the parameters passed to the API\n\n");
}


void fla_test_output_op_struct( char* op_str, integer op )
{
    fla_test_output_info( "%s LAPACK  %"FT_IS"\n", op_str, op );
}

/* THis function checks if operations file specified any run only option */
integer fla_test_check_run_only( FILE* input_stream, integer* op, char* buffer)
{
    integer run_only_flag = 1, i;
    integer test_group_count = sizeof(API_test_group) / sizeof(API_test_group[0]);

    /*checks if operations file specified any sub-group of API is enable*/
    for(i = 0; i < test_group_count; i++)
    {
        fla_test_read_tests_for_op( input_stream, op, buffer);
        if((*op) == 1 )
        {
                run_only_flag = 3;
                fseek( input_stream, 0, SEEK_SET );
                return run_only_flag;
        }
    }

    while(fla_test_read_tests_for_op( input_stream, op, buffer))
    {
        run_only_flag = fla_max(run_only_flag, *op);
    }

    fseek( input_stream, 0, SEEK_SET );

    return run_only_flag;
}

/* This function reads group_id for given group name*/
integer fla_test_get_group_id(char *buffer)
{
    integer i;
    integer test_group_count = sizeof(API_test_group) / sizeof(API_test_group[0]);

    for(i = 0; i < test_group_count; i++)
    {
        if(!strcmp(buffer, API_test_group[i]))
            return i;
    }

    return -1;
}

/* This functiom extract enable option and API name from operation file */
integer fla_test_read_tests_for_op( FILE* input_stream, integer* op, char* buffer )
{
    char temp[ INPUT_BUFFER_SIZE ];

    // We want to read at least one line, so we use a do-while loop.
    do
    {
        // Read the next line into a temporary buffer and check success.
        if ( fgets( temp, INPUT_BUFFER_SIZE-1, input_stream ) == NULL )
            return 0;
    }
    // We continue to read lines into buffer until the line is neither
    // commented nor blank.
    while ( temp[0] == COMMENT_CHAR || temp[0] == '\n' ||
            temp[0] == ' '          || temp[0] == '\t' );

    // Save the string in temp, up to first white space character, into buffer.
    sscanf( temp, "%"FT_IS" %s", op, buffer );

    return 1;

}


/* This function reads operation file line by line*/
void fla_test_read_next_line( char* buffer, FILE* input_stream )
{
    char temp[ INPUT_BUFFER_SIZE ];

    // We want to read at least one line, so we use a do-while loop.
    do
    {
        // Read the next line into a temporary buffer and check success.
        if ( fgets( temp, INPUT_BUFFER_SIZE-1, input_stream ) == NULL )
        {
            if ( feof( input_stream ) )
                fla_test_output_error( "Error reading input file: encountered unexpected EOF." );
            else
                fla_test_output_error( "Error (non-EOF) reading input file." );
        }
    }
    // We continue to read lines into buffer until the line is neither
    // commented nor blank.
    while ( temp[0] == COMMENT_CHAR || temp[0] == '\n' ||
            temp[0] == ' '          || temp[0] == '\t' );

    // Save the string in temp, up to first white space character, into buffer.
    sscanf( temp, "%s ", buffer );
}


/* This function reads parameters needed for Linear solver APIs
   from the config settings file 'LIN_SLVR.dat' and saves in the
   'lin_solver_paramslist' structure array   */
void fla_test_read_linear_param ( const char *file_name, test_params_t* params )
{
    FILE *fp;
    integer i, j;
    char line[20];
    char *str;
    char eol;
    integer num_tests, ndata_types;
    integer num_ranges;

    str = &line[0];
    fp = fopen( file_name, "r");
    if (fp == NULL){
    printf("Error: Lin solver config file missing. Exiting.. \n");
    exit(-1);
    }

    /* Read the number of Tests */
    fscanf(fp, "%s", &line[0]);
    fscanf(fp, "%"FT_IS"", &num_tests);

    num_ranges = num_tests;
    for (i=0; i<NUM_SUB_TESTS; i++){
        params->lin_solver_paramslist[i].num_tests = num_tests;
    }

    fscanf(fp, "%s", &line[0]); // Range_start
    for (i=0; i < NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].m_range_start) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_end
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].m_range_end) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_step_size
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].m_range_step_size) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_start

    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].n_range_start) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_end
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].n_range_end) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_step_size
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].n_range_step_size) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // leading dimension for A
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].lda) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // leading dimension for B
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].ldb) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // leading dimension LDAB
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].ldab) );
        CHECK_LINE_SKIP ();
    }

    for (i=0; i<NUM_SUB_TESTS; i++){
        params->lin_solver_paramslist[i].num_ranges = num_ranges;
    }

    fscanf(fp, "%s", &line[0]); // Numer of repeats
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].num_repeats) );
    }

    ndata_types = NUM_SUB_TESTS;
    fscanf(fp, "%s", &line[0]); // Datatypes
    str = &line[0];
    for( i = 0; i < NUM_SUB_TESTS; i++ )
    {
        fscanf(fp, "%s", str);
        for( j = 0; j < NUM_SUB_TESTS; j++ )
        {
            params->lin_solver_paramslist[j].data_types_char[i] = *str;
            params->lin_solver_paramslist[j].data_types[i] = get_datatype(*str);
        }
        eol = fgetc(fp);
        if((eol == '\r') || (eol == '\n')){
            ndata_types = ( (i+1) < ndata_types)? (i+1):ndata_types;
            break;
        }
    }
    for (i=0; i<NUM_SUB_TESTS; i++){
        params->lin_solver_paramslist[i].num_data_types = ndata_types;
    }

    fscanf(fp, "%s", &line[0]); // Matrix Layout (row or col major)
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].matrix_layout) );
        CHECK_LINE_SKIP ();
    }

    str = &line[0];
    fscanf(fp, "%s", str);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->lin_solver_paramslist[i].transr = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->lin_solver_paramslist[i].Uplo = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].nrhs) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i = 0; i < NUM_SUB_TESTS; i++) {
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].ncolm));
        CHECK_LINE_SKIP();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].kl) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].ku) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].kd) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->lin_solver_paramslist[i].diag = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->lin_solver_paramslist[i].fact = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->lin_solver_paramslist[i].equed = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->lin_solver_paramslist[i].symm = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%f", &(params->lin_solver_paramslist[i].solver_threshold) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->lin_solver_paramslist[i].equed_porfsx = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].n_err_bnds_porfsx) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].nparams_porfsx) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->lin_solver_paramslist[i].norm_gbcon = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].kl_gbcon) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].ku_gbcon) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->lin_solver_paramslist[i].ldab_gbcon) );
        CHECK_LINE_SKIP ();
    }
    fclose(fp);

}


/* This function reads parameters needed for Eigen APIs
   from the config settings file 'EIG_PARAMS.dat' and saves in the
   'params->eig_sym_paramslist' structure array   */
void fla_test_read_sym_eig_params( const char *file_name , test_params_t* params )
{
    FILE *fp;
    integer i, j;
    char line[20], eol;
    char *str, c[20];
    integer num_tests;
    integer num_ranges;
    integer ndata_types = NUM_SUB_TESTS;

    str = &c[0];
    fp = fopen( file_name, "r");
    if (fp == NULL){
    printf("Error: Symmetric EIG params config file missing. Exiting.. \n");
    exit(-1);
    }

    /* Read the number of Ranges */
    fscanf(fp, "%s", &line[0]);
    fscanf(fp, "%"FT_IS"", &num_tests);
    for (i=0; i<NUM_SUB_TESTS; i++){
        params->eig_sym_paramslist[i].num_tests = num_tests;
    }

    num_ranges = num_tests;
    fscanf(fp, "%s", &line[0]); // Range_start

    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].m_range_start) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_end
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].m_range_end) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_step_size
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].m_range_step_size) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_end
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].n_range_start) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_end
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].n_range_end) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_step_size
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].n_range_step_size) );
        CHECK_LINE_SKIP ();
    }

    for (i=0; i<NUM_SUB_TESTS; i++){
        params->eig_sym_paramslist[i].num_ranges = num_ranges;
    }

    fscanf(fp, "%s", &line[0]); // number of repeats
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].num_repeats) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    str = &line[0];

    ndata_types = NUM_SUB_TESTS;
    for( i = 0; i < NUM_SUB_TESTS; i++ )
    {
        fscanf(fp, "%s", str); // num data types
        for( j = 0; j < NUM_SUB_TESTS; j++ )
        {
            params->eig_sym_paramslist[j].data_types_char[i] = *str;
            params->eig_sym_paramslist[j].data_types[i] = get_datatype(*str);
        }
        eol = fgetc(fp);
        if((eol == '\r') || (eol == '\n')){
            ndata_types = ( (i+1) < ndata_types)? (i+1):ndata_types;
            break;
        }
    }

    for (i=0; i<NUM_SUB_TESTS; i++){
        params->eig_sym_paramslist[i].num_data_types = ndata_types;
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].matrix_layout) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", str);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].trans = *str;
        CHECK_LINE_SKIP ();
    }
    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].uplo = *str;
        CHECK_LINE_SKIP ();
    }
    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].job = *str;
        CHECK_LINE_SKIP ();
    }
    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].jobz = *str;
        CHECK_LINE_SKIP ();
    }
    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].job_seqr = *str;
        CHECK_LINE_SKIP ();
    }
    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].vect = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].nrhs) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].lda) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].ldb) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].ldz) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].nb) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].ldt) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].k) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].isgn) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].compz = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].compz_hseqr = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].kb) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].itype) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].vect_rd = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].side = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].eigsrc = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].initv = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].norm = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].diag = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_sym_paramslist[i].storev = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].tsize) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].ilo) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].ihi) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_sym_paramslist[i].threshold_value) );
        CHECK_LINE_SKIP ();
    }

    fclose(fp);

}

/* This function reads parameters needed for Non symmetric Eigen APIs
   from the config settings file 'EIG_NSYM_PARAMS.dat' and saves in the
   'params->eig_non_sym_paramslist' structure array   */
/* This function reads parameters needed for Non symmetric Eigen APIs
   from the config settings file 'EIG_NSYM_PARAMS.dat' and saves in the
   'params->eig_non_sym_paramslist' structure array   */

void fla_test_read_non_sym_eig_params( const char *file_name , test_params_t* params )
{
    FILE *fp;
    integer i, j;
    char line[20], eol;
    char *str = &line[0];
    integer num_tests;
    integer num_ranges;
    integer ndata_types = NUM_SUB_TESTS;
    fp = fopen( file_name, "r");
    if (fp == NULL){
    printf("Error: EIG non symmetric API params config file missing. Exiting.. \n");
    exit(-1);
    }

    /* Read the number of Tests */
    fscanf(fp, "%s", &line[0]);
    fscanf(fp, "%"FT_IS"", &num_tests);
    for (i=0; i<NUM_SUB_TESTS; i++){
        params->eig_non_sym_paramslist[i].num_tests = num_tests;
    }

    num_ranges = num_tests;

    fscanf(fp, "%s", &line[0]); // Range_start
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].m_range_start) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_end
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].m_range_end) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_step_size
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].m_range_step_size) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_end
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].n_range_start) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_end
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].n_range_end) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_step_size
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].n_range_step_size) );
        CHECK_LINE_SKIP ();
    }

    for (i=0; i<NUM_SUB_TESTS; i++){
        params->eig_non_sym_paramslist[i].num_ranges = num_ranges;
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].lda) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].ldb) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].ldvl) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].ldvr) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // number of repeats
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].num_repeats) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for( i = 0; i < NUM_SUB_TESTS; i++ )
    {
        fscanf(fp, "%s", str); // num data types
        for( j = 0; j < NUM_SUB_TESTS; j++ )
        {
            params->eig_non_sym_paramslist[j].data_types_char[i] = *str;
            params->eig_non_sym_paramslist[j].data_types[i] = get_datatype(*str);
        }
        eol = fgetc(fp);
        if((eol == '\r') || (eol == '\n')){
            ndata_types = ( (i+1) < ndata_types)? (i+1):ndata_types;
            break;
        }
    }

    for (i=0; i<NUM_SUB_TESTS; i++){
        params->eig_non_sym_paramslist[i].num_data_types = ndata_types;
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].matrix_layout) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].howmny = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].initv = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].job_seqr = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].eigsrc = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].initv = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].job = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].howmny_trsna = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].job_trsen = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].compq = *str;
        CHECK_LINE_SKIP ();
    }

    /* Reading config params for 'trsyl' API  */
    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].trana_real = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].trana_complex = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].tranb_real= *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].tranb_complex = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].isgn) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%f", &(params->eig_non_sym_paramslist[i].gghrd_threshold) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%f", &(params->eig_non_sym_paramslist[i].ggbal_threshold) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%f", &(params->eig_non_sym_paramslist[i].GenNonSymEigProblem_threshold) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].compq_hgeqz = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].compz_hgeqz = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].side_tgevc = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].jobvsl = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].jobvsr = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].sort = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].sense_ggesx = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].balance_ggevx = *str;
        CHECK_LINE_SKIP ();
    }
    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].sense_ggevx = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].sort_gees = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].wantz) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].wantq) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->eig_non_sym_paramslist[i].tgsen_ijob) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->eig_non_sym_paramslist[i].unmhr_trans = *str;
        CHECK_LINE_SKIP ();
    }


    fclose(fp);
}



/* This function reads parameters needed for SVD APIs
   from the config settings file 'SVD.dat' and saves in the
   'params->svd_paramslist' structure array   */
void fla_test_read_svd_params ( const char *file_name, test_params_t* params )
{
    FILE *fp;
    integer i, j;
    char line[25];
    char *str;
    char eol;
    integer num_tests;
    integer ndata_types;
    integer num_ranges;

    str = &line[0];
    fp = fopen( file_name, "r");
    if (fp == NULL){
    printf("Error: SVD config file missing. Exiting.. \n");
    exit(-1);
    }

    fscanf(fp, "%s", &line[0]);
    fscanf(fp, "%"FT_IS"", &num_tests);
    for (i=0; i<NUM_SUB_TESTS; i++){
        params->svd_paramslist[i].num_tests = num_tests;
    }

    num_ranges = num_tests;
    fscanf(fp, "%s", &line[0]); // Range_start
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].m_range_start) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_end
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].m_range_end) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_step_size
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].m_range_step_size) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_start
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].n_range_start) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_end
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].n_range_end) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]); // Range_step_size
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].n_range_step_size) );
        CHECK_LINE_SKIP ();
    }

    for (i=0; i<NUM_SUB_TESTS; i++){
        params->svd_paramslist[i].num_ranges = num_ranges;
    }

    fscanf(fp, "%s", &line[0]); // number of repeats
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].num_repeats) );
        CHECK_LINE_SKIP ();
    }

    ndata_types = NUM_SUB_TESTS;
    fscanf(fp, "%s", &line[0]);
    for( i = 0; i < NUM_SUB_TESTS; i++ )
    {
        fscanf(fp, "%s", str); // num data types
        for( j = 0; j < NUM_SUB_TESTS; j++ )
        {
            params->svd_paramslist[j].data_types_char[i] = *str;
            params->svd_paramslist[j].data_types[i] = get_datatype(*str);
        }
        eol = fgetc(fp);
        if((eol == '\r') || (eol == '\n')){
            ndata_types = ( (i+1) < ndata_types)? (i+1):ndata_types;
            break;
        }
    }

    for (i=0; i<NUM_SUB_TESTS; i++){
        params->svd_paramslist[i].num_data_types = ndata_types;
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].matrix_layout) );
        CHECK_LINE_SKIP ();
    }

    str = &line[0];
    fscanf(fp, "%s", str);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobu = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobv = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobq = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].m) );
        CHECK_LINE_SKIP ();
    }


    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].p) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].n) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%f", &(params->svd_paramslist[i].tola) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%f", &(params->svd_paramslist[i].tolb) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%f", &(params->svd_paramslist[i].svd_threshold) );
        CHECK_LINE_SKIP ();
    }

    str = &line[0];
    fscanf(fp, "%s", str);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobu_gesvd = *str;
        CHECK_LINE_SKIP ();
    }

    str = &line[0];
    fscanf(fp, "%s", str);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobvt_gesvd = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].joba_gejsv = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobu_gejsv = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobv_gejsv = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobr_gejsv = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobt_gejsv = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobp_gejsv = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].m_gejsv) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].n_gejsv) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].joba_gesvj = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobu_gesvj = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobv_gesvj = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].m_gesvj) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].n_gesvj) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].mv_gesvj) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%f", &(params->svd_paramslist[i].ctol_gesvj) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", &(params->svd_paramslist[i].jobu_gesvdx) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", &(params->svd_paramslist[i].jobvt_gesvdx) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", &(params->svd_paramslist[i].range_gesvdx) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].il) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].iu) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%f", &(params->svd_paramslist[i].vl) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%f", &(params->svd_paramslist[i].vu) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].joba_gesvdq = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobu_gesvdq = *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%s", str);
        params->svd_paramslist[i].jobv_gesvdq= *str;
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].lda) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].ldu) );
        CHECK_LINE_SKIP ();
    }

    fscanf(fp, "%s", &line[0]);
    for (i=0; i<NUM_SUB_TESTS; i++){
        fscanf(fp, "%"FT_IS"", &(params->svd_paramslist[i].ldvt) );
        CHECK_LINE_SKIP ();
    }

    fclose(fp);
}


void fla_test_output_info( char* message, ... )
{
    FILE* output_stream = stdout;
    va_list args;

    // Initialize variable argument environment.
    va_start( args, message );

    // Parse the received message and print its components.
    fla_test_parse_message( output_stream, message, args );

    // Shutdown variable argument environment and clean up stack.
    va_end( args );

    // Flush the output stream.
    fflush( output_stream );
}


void fla_test_output_error( char* message, ... )
{
    FILE* output_stream = stderr;
    va_list args;

    fprintf( output_stream, "%s: *** error ***: ", fla_test_binary_name );

    // Initialize variable argument environment.
    va_start( args, message );

    // Parse the received message and print its components.
    fla_test_parse_message( output_stream, message, args );

    // Shutdown variable argument environment and clean up stack.
    va_end( args );

    // Flush the output stream.
    fflush( output_stream );

    // Exit.
    exit(1);
}


void fla_test_parse_message( FILE* output_stream, char* message, va_list args )
{
    integer           c, cf;
    char          format_spec[20];
    uinteger  the_uint;
    integer           the_int;
    double        the_double;
    char*         the_string;
    char          the_char;

    // Begin looping over message to insert variables wherever there are
    // format specifiers.
    for ( c = 0; message[c] != '\0'; )
    {
        if ( message[c] != '%' )
        {
            fprintf( output_stream, "%c", message[c] );
            c += 1;
        }
        else if ( message[c] == '%' && message[c+1] == '%' ) // handle escaped '%' chars.
        {
            fprintf( output_stream, "%c", message[c] );
            c += 2;
        }
        else
        {
            // Save the format string if there is one.
            format_spec[0] = '%';
            for ( c += 1, cf = 1; strchr( "udefsc", message[c] ) == NULL; ++c, ++cf )
            {
                format_spec[cf] = message[c];
            }

            // Add the final type specifier, and null-terminate the string.
            format_spec[cf] = message[c];
            format_spec[cf+1] = '\0';

            // Switch based on type, since we can't predict what will
            // va_args() will return.
            switch ( message[c] )
            {
                case 'u':
                the_uint = va_arg( args, uinteger );
                fprintf( output_stream, format_spec, the_uint );
                break;

                case 'd':
                the_int = va_arg( args, integer );
                fprintf( output_stream, format_spec, the_int );
                break;

                case 'e':
                the_double = va_arg( args, double );
                fprintf( output_stream, format_spec, the_double );
                break;

                case 'f':
                the_double = va_arg( args, double );
                fprintf( output_stream, format_spec, the_double );
                break;

                case 's':
                the_string = va_arg( args, char* );
                fprintf( output_stream, format_spec, the_string );
                break;

                case 'c':
                the_char = va_arg( args, integer );
                fprintf( output_stream, "%c", the_char );
                break;
            }

            // Move to next character past type specifier.
            c += 1;
        }
    }
}


void fla_test_execute_cli_api( integer argc, char** argv, test_params_t *params )
{
    integer i, test_api_count;
    char s_name[MAX_FUNC_STRING_LENGTH];

    test_api_count = sizeof(API_test_functions) / sizeof(API_test_functions[0]);
    fla_test_output_info( "%2sAPI%13s DATA_TYPE%6s SIZE%9s FLOPS%9s TIME%9s ERROR%9s STATUS\n", "", "", "", "", "", "", "" );
    fla_test_output_info( "%1s=====%12s===========%4s========%7s=======%7s========%6s==========%6s========\n", "", "", "", "", "", "", "" );

    /* Check if the specified API is supported in test suite */
    strcpy(s_name, argv[1]);
    for( i = 0; i < strlen(s_name); i++)
    {
        s_name[i] = tolower(s_name[i]);
    }
    for( i = 0; i < test_api_count; i++)
    {
        if(!strcmp(API_test_functions[i].ops, s_name))
        {
            API_test_functions[i].fp(argc, argv, params);
            break;
        }
    }

    if(total_tests == 0)
        printf("\nNo test was run, give valid arguments\n");
}


char* fla_test_get_string_for_result( double residual, integer datatype, double thresh )
{
    char* r_val;

    if ( datatype == FLOAT )
    {
        if      ( residual == DBL_MAX )   r_val = fla_test_incomplete_string;
        else if ( residual == DBL_MIN )   r_val = fla_test_invalid_string;
        else if ( residual > thresh )     r_val = fla_test_fail_string;
        else                              r_val = fla_test_pass_string;
    }
    else if ( datatype == DOUBLE )
    {
        if      ( residual == DBL_MAX )   r_val = fla_test_incomplete_string;
        else if ( residual == DBL_MIN )   r_val = fla_test_invalid_string;
        else if ( residual > thresh )     r_val = fla_test_fail_string;
        else                              r_val = fla_test_pass_string;
    }
    else if ( datatype == COMPLEX )
    {
        if      ( residual == DBL_MAX )   r_val = fla_test_incomplete_string;
        else if ( residual == DBL_MIN )   r_val = fla_test_invalid_string;
        else if ( residual > thresh )     r_val = fla_test_fail_string;
        else                              r_val = fla_test_pass_string;
    }
    else
    {
        if      ( residual == DBL_MAX )   r_val = fla_test_incomplete_string;
        else if ( residual == DBL_MIN )   r_val = fla_test_invalid_string;
        else if ( residual > thresh )     r_val = fla_test_fail_string;
        else                              r_val = fla_test_pass_string;
    }

    return r_val;
}


void fla_test_init_strings( void )
{
    sprintf( fla_test_pass_string, "PASS" );
    sprintf( fla_test_warn_string, "MARGINAL" );
    sprintf( fla_test_fail_string, "FAIL" );
    sprintf( fla_test_incomplete_string, "INCOMPLETE" );
    sprintf( fla_test_invalid_string, "INVALID_LDA" );
    sprintf( fla_test_storage_format_string, "Row(r) and General(g) storage format is not supported by External LAPACK interface" );
    sprintf( fla_test_stor_chars, STORAGE_SCHEME_CHARS );
}


void fla_test_op_driver( char*         func_str,
                         integer       sqr_inp,
                         test_params_t *params,
                         integer       api_type,
                         void (*f_exp) (test_params_t *,  // params
                                           integer,          // datatype
                                           integer,          // p_cur
                                           integer,          // q_cur
                                           integer,          // pci (param combo counter)
                                           integer,          // n_repeats
                                           double*,          // perf
                                           double*,          //time
                                           double* ) )       // residual
{
    integer n_datatypes = params->n_datatypes;
    integer n_repeats, ith;
    integer num_ranges, range_loop_counter;
    integer p_first, p_max, p_inc;
    integer q_first, q_max, q_inc;
    integer dt, p_cur, q_cur;
    char    datatype_char;
    integer datatype;
    double thresh;
    double perf_max_val, time_min_val, residual_max_val;
    double* perf = (double*)malloc(n_threads*sizeof(double));
    double* time = (double*)malloc(n_threads*sizeof(double));
    double* residual = (double*)malloc(n_threads*sizeof(double));

    fla_test_output_info( "%2sAPI%13s DATA_TYPE%6s SIZE%9s FLOPS%9s TIME%9s ERROR%9s STATUS\n", "", "", "", "", "", "", "" );
    fla_test_output_info( "%1s=====%12s===========%4s========%7s=======%7s========%6s==========%6s========\n", "", "", "", "", "", "", "" );
    switch (api_type)
    {
        case LIN:
            num_ranges          = params->lin_solver_paramslist[0].num_ranges;
            break;

        case EIG_SYM:
            num_ranges          = params->eig_sym_paramslist[0].num_ranges;
            break;

        case EIG_NSYM:
            num_ranges          = params->eig_non_sym_paramslist[0].num_ranges;
            break;

        case SVD:
            num_ranges          = params->svd_paramslist[0].num_ranges;
            break;

        default:
            fla_test_output_error( "Invalid API type. Exiting...\n" );
            return;
    }

    // Loop over the requested ranges.
    for ( range_loop_counter = 0; range_loop_counter < num_ranges; ++range_loop_counter )
    {
        switch (api_type)
        {
            case LIN:
                p_first               = params->lin_solver_paramslist[range_loop_counter].m_range_start;
                p_max                 = params->lin_solver_paramslist[range_loop_counter].m_range_end;
                p_inc                 = params->lin_solver_paramslist[range_loop_counter].m_range_step_size;
                q_first               = params->lin_solver_paramslist[range_loop_counter].n_range_start;
                q_max                 = params->lin_solver_paramslist[range_loop_counter].n_range_end;
                q_inc                 = params->lin_solver_paramslist[range_loop_counter].n_range_step_size;
                thresh                = params->lin_solver_paramslist[range_loop_counter].solver_threshold;
                params->datatype      = params->lin_solver_paramslist[range_loop_counter].data_types;
                params->datatype_char = params->lin_solver_paramslist[range_loop_counter].data_types_char;
                n_repeats             = params->lin_solver_paramslist[range_loop_counter].num_repeats;
                n_datatypes           = params->lin_solver_paramslist[range_loop_counter].num_data_types;
                break;

            case EIG_SYM:
                p_first               = params->eig_sym_paramslist[range_loop_counter].m_range_start;
                p_max                 = params->eig_sym_paramslist[range_loop_counter].m_range_end;
                p_inc                 = params->eig_sym_paramslist[range_loop_counter].m_range_step_size;
                q_first               = p_first;
                q_max                 = p_max;
                q_inc                 = p_inc;
                thresh                = params->eig_sym_paramslist[range_loop_counter].threshold_value;
                params->datatype      = params->eig_sym_paramslist[range_loop_counter].data_types;
                params->datatype_char = params->eig_sym_paramslist[range_loop_counter].data_types_char;
                n_repeats             = params->eig_sym_paramslist[range_loop_counter].num_repeats;
                n_datatypes           = params->eig_sym_paramslist[range_loop_counter].num_data_types;
                break;

            case EIG_NSYM:
                p_first               = params->eig_non_sym_paramslist[range_loop_counter].m_range_start;
                p_max                 = params->eig_non_sym_paramslist[range_loop_counter].m_range_end;
                p_inc                 = params->eig_non_sym_paramslist[range_loop_counter].m_range_step_size;
                q_first               = params->eig_non_sym_paramslist[range_loop_counter].n_range_start;
                q_max                 = params->eig_non_sym_paramslist[range_loop_counter].n_range_end;
                q_inc                 = params->eig_non_sym_paramslist[range_loop_counter].n_range_step_size;
                thresh                = params->eig_non_sym_paramslist[range_loop_counter].gghrd_threshold;
                params->datatype      = params->eig_non_sym_paramslist[range_loop_counter].data_types;
                params->datatype_char = params->eig_non_sym_paramslist[range_loop_counter].data_types_char;
                n_repeats             = params->eig_non_sym_paramslist[range_loop_counter].num_repeats;
                n_datatypes           = params->eig_non_sym_paramslist[range_loop_counter].num_data_types;
                break;

            case SVD:
                p_first               = params->svd_paramslist[range_loop_counter].m_range_start;
                p_max                 = params->svd_paramslist[range_loop_counter].m_range_end;
                p_inc                 = params->svd_paramslist[range_loop_counter].m_range_step_size;
                q_first               = params->svd_paramslist[range_loop_counter].n_range_start;
                q_max                 = params->svd_paramslist[range_loop_counter].n_range_end;
                q_inc                 = params->svd_paramslist[range_loop_counter].n_range_step_size;
                thresh                = params->svd_paramslist[range_loop_counter].svd_threshold;
                params->datatype      = params->svd_paramslist[range_loop_counter].data_types;
                params->datatype_char = params->svd_paramslist[range_loop_counter].data_types_char;
                n_repeats             = params->svd_paramslist[range_loop_counter].num_repeats;
                n_datatypes           = params->svd_paramslist[range_loop_counter].num_data_types;
                break;
            default:
                return;
        }

        /* Loop over the requested datatypes. */
        for ( dt = 0; dt < n_datatypes; ++dt )
        {
            datatype              = params->datatype[dt];
            datatype_char         = params->datatype_char[dt];

            /* Loop over the requested problem sizes */
            for ( p_cur = p_first, q_cur = q_first; (p_cur <= p_max && q_cur <= q_max); p_cur += p_inc, q_cur += q_inc )
            {
                if ( n_threads > 1)
                {
#pragma omp parallel num_threads( n_threads )
#pragma omp for
                    for ( ith = 0; ith < n_threads; ith++ )
                    {
                        f_exp(params, datatype, p_cur, q_cur, range_loop_counter, n_repeats, (perf+ith), (time+ith), (residual+ith));
                    }

                    get_max(DOUBLE, (void*)residual, (void*)&residual_max_val, n_threads);
                    get_min(DOUBLE, (void*)time, (void*)&time_min_val, n_threads);
                    get_max(DOUBLE, (void*)perf, (void*)&perf_max_val, n_threads);

                    fla_test_print_status(func_str, datatype_char, sqr_inp, p_cur, q_cur, residual_max_val, thresh, time_min_val, perf_max_val);
                }
                else
                {
                    f_exp(params, datatype, p_cur, q_cur, range_loop_counter, n_repeats, perf, time, residual);
                    fla_test_print_status(func_str, datatype_char, sqr_inp, p_cur, q_cur, *residual, thresh, *time, *perf);
                }

            }
        }

        fla_test_output_info( "\n" );
    }

    free(perf);
    free(residual);
    free(time);
}

void fla_test_print_status(char* func_str,
                           char datatype_char,
                           integer sqr_inp,
                           integer p_cur,
                           integer q_cur,
                           double residual,
                           double thresh,
                           double time,
                           double perf)
{
    char func_param_str[64];
    char blank_str[32];
    char scale[3] = "";
    integer n_spaces, datatype;
    char*   pass_str;

    datatype = get_datatype(datatype_char);

    pass_str = fla_test_get_string_for_result( residual, datatype, thresh );

    fla_test_build_function_string( func_str, NULL, func_param_str ); 

    n_spaces = MAX_FUNC_STRING_LENGTH - strlen( func_param_str );
    fill_string_with_n_spaces( blank_str, n_spaces );

    fla_test_get_time_unit(scale , &time);

    if ( sqr_inp == SQUARE_INPUT )
    {
        if (pass_str[0] == 'P' || pass_str[0] == 'F')
            fla_test_output_info(" %s%s  %c  %10"FT_IS" x %-9"FT_IS" %-10.2lf  %6.2lf %-7s  %-7.2le   %10s\n",
                                func_param_str, blank_str,
                                datatype_char,
                                p_cur, p_cur, perf, time, scale, residual, pass_str );
        else
           fla_test_output_info(" %s%s  %c  %10"FT_IS" x %-9"FT_IS" %-10.2lf  %6.2lf %-7s  %-7.2le   %12s\n",
                                func_param_str, blank_str,
                                datatype_char,
                                p_cur, p_cur, perf, time, scale, residual, pass_str ); 
    }
    else
    {
        if (pass_str[0] == 'P' || pass_str[0] == 'F')
            fla_test_output_info(" %s%s  %c  %10"FT_IS" x %-9"FT_IS" %-10.2lf  %6.2lf %-7s  %-7.2le   %10s\n",
                                func_param_str, blank_str,
                                datatype_char,
                                p_cur, q_cur, perf, time, scale, residual, pass_str );
        else
            fla_test_output_info(" %s%s  %c  %10"FT_IS" x %-9"FT_IS" %-10.2lf  %6.2lf %-7s  %-7.2le   %12s\n",
                                func_param_str, blank_str,
                                datatype_char,
                                p_cur, q_cur, perf, time, scale, residual, pass_str );
    }

    total_tests++;
    tests_passed[datatype - FLOAT] += (pass_str[0] == 'P');
    tests_failed[datatype - FLOAT] += (pass_str[0] == 'F');
    tests_incomplete[datatype - FLOAT] += (pass_str[0] == 'I');
    total_failed_tests += (pass_str[0] == 'F');
    total_incomplete_tests += (pass_str[0] == 'I');
}


void fla_test_build_function_string( char*        func_base_str,
                                        char*     impl_var_str,
                                        char*     func_str )
{
    sprintf( func_str, "%s", func_base_str );
}


void fill_string_with_n_spaces( char* str, integer n_spaces )
{
    integer i;

    for ( i = 0; i < n_spaces; ++i )
        sprintf( &str[i], " " );
}


/* This function is used to clock the execution time of APIs*/
double fla_test_clock()
{
#ifdef _WIN32
    LARGE_INTEGER clock_freq = {0};
    LARGE_INTEGER clock_val;

    QueryPerformanceFrequency( &clock_freq );
    QueryPerformanceCounter( &clock_val );

    return ( ( double) clock_val.QuadPart / ( double) clock_freq.QuadPart );
#else
    double the_time, norm_sec;
    struct timespec ts;

    clock_gettime( CLOCK_MONOTONIC, &ts );

    if ( ref_time_sec == 0.0 )
        ref_time_sec = ( double ) ts.tv_sec;

    norm_sec = ( double ) ts.tv_sec - ref_time_sec;

    the_time = norm_sec + ts.tv_nsec * 1.0e-9;

    return the_time;
#endif
}


/* This function sets the unit of time*/
void fla_test_get_time_unit(char * scale , double * time)
{
    if(*time >= 1)
    {
        scale[0] = 's';
        scale[1] = '\0';
        return ;
    }

    if ((*time < 1) && (*time > 0.001))
    {
        scale[0]='m';
        *time *= 1000;
    }
    else if ((*time <0.001) && (*time > 0.000001))
    {
        scale[0]='u';
        *time *= 1000000;
    }
    else if ((*time < 0.000001) && (*time > 0.000000001))
    {
        scale[0]='n';
        *time *= 1000000000;
    }
    else if ((*time < 0.000000001) && (*time > 0.000000000001))
    {
        scale[0]='p';
        *time *= 1000000000000;
    }

    scale[1] = 's';
    scale[2] = '\0';
}
