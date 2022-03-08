/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_libflame.h"
#include "test_routines.h"


// Global variables.
char fla_test_binary_name[ MAX_BINARY_NAME_LENGTH + 1 ];
char fla_test_pass_string[ MAX_PASS_STRING_LENGTH + 1 ];
char fla_test_warn_string[ MAX_PASS_STRING_LENGTH + 1 ];
char fla_test_fail_string[ MAX_PASS_STRING_LENGTH + 1 ];
char fla_test_storage_format_string[ 200 ];
char fla_test_stor_chars[ NUM_STORAGE_CHARS + 1 ];
clock_t timer;


int main( int argc, char** argv )
{
    test_params_t params;
    test_ops_t    ops;
    integer vers_major, vers_minor, vers_patch;

    ilaver_(&vers_major, &vers_minor, &vers_patch);

    printf(" LibFlame version: %"FS".%"FS".%"FS" \n", vers_major, vers_minor, vers_patch);

    // Initialize some strings.
    fla_test_init_strings();

    // Parse the command line parameters.
    fla_test_parse_command_line( argc, argv );

    /*Read Linear API parameters from config file */
    fla_test_read_linear_param( LINEAR_PARAMETERS_FILENAME, &params );

    /*Read eigen parameters from config file */
    fla_test_read_sym_eig_params( SYM_EIG_PARAMETERS_FILENAME, &params );
    fla_test_read_non_sym_eig_params( NON_SYM_EIG_PARAMETERS_FILENAME, &params );

    /*Read SVD parameters from config file */
    fla_test_read_svd_params ( SVD_PARAMETERS_FILENAME, &params  );

    // Test the LAPACK-level operations.
    fla_test_lapack_suite( OPERATIONS_FILENAME, &params, ops );

    return 0;
}


void fla_test_lapack_suite( char* input_filename, test_params_t *params, test_ops_t ops )
{
    char buffer[ INPUT_BUFFER_SIZE ];
	char  temp[ INPUT_BUFFER_SIZE ];
	integer i;
    integer op;
    FILE* input_stream;


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

    // Read the number of repeats.
    fla_test_read_next_line( buffer, input_stream );
    sscanf( buffer, "%u ", &(params->n_repeats) );

    // Read the datatypes to test. We should have at most four: 's', 'd', 'c',
    // and 'z'.
    fla_test_read_next_line( buffer, input_stream );
    sscanf( buffer, "%s ", temp );

    params->n_datatypes = strlen( temp );

    if ( params->n_datatypes > MAX_NUM_DATATYPES )
    {
        fla_test_output_error( "Detected too many datatype requests (%u) in input file.\n", params->n_datatypes );
    }

    for( i = 0; i < params->n_datatypes; ++i )
    {
        if ( temp[i] == 's' ) params->datatype[i] = FLOAT;
        else if ( temp[i] == 'd' ) params->datatype[i] = DOUBLE;
        else if ( temp[i] == 'c' ) params->datatype[i] = COMPLEX;
        else if ( temp[i] == 'z' ) params->datatype[i] = DOUBLE_COMPLEX;

        params->datatype_char[i] = temp[i];
    }

    fla_test_read_next_line( buffer, input_stream );

    while(fla_test_read_tests_for_op( input_stream, &op, buffer))
    {
        if(op)
        {
            for( i=0; i < test_api_count; i++)
            {
                if(!strcmp(API_test_functions[i].ops, buffer))
                {
                    API_test_functions[i].fp(params);
                }
            }

            op = 0;
        }
    }

    fclose( input_stream );
}


void fla_test_output_op_struct( char* op_str, int op )
{
    fla_test_output_info( "%s LAPACK  %d\n", op_str, op );
}


int fla_test_read_tests_for_op( FILE* input_stream, int* op, char* buffer )
{
    char temp[ INPUT_BUFFER_SIZE ];

    // We want to read at least one line, so we use a do-while loop.
    do
    {
        // Read the next line into a temporary buffer and check success.
        if ( fgets( temp, INPUT_BUFFER_SIZE-1, input_stream ) == NULL )
        {
            return 0;
        }
    }
    // We continue to read lines into buffer until the line is neither
    // commented nor blank.
    while ( temp[0] == COMMENT_CHAR || temp[0] == '\n' ||
            temp[0] == ' '          || temp[0] == '\t' );


    // Save the string in temp, up to first white space character, into buffer.
    sscanf( temp, "%d %s", op, buffer );

    return 1;

}


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
   integer i;
   char line[20];
   char *str;
   integer num_tests;
   integer mode;
   integer num_ranges;

   str = &line[0];
   fp = fopen( file_name, "r");
   if (fp == NULL){
    printf("Error: Lin solver config file missing. Exiting.. \n");
    exit(-1);
   }

   /* Read the mode */
   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &mode);
   fscanf(fp, "%*[^\n]\n");

   /* Read the number of Ranges */
   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &num_ranges);

   fscanf(fp, "%s", &line[0]); // Range_start

   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].m_range_start) );
      params->lin_solver_paramslist[i].num_ranges = num_ranges;
      params->lin_solver_paramslist[i].mode = mode;
   }

   fscanf(fp, "%s", &line[0]); // Range_end
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].m_range_end) );
   }

   fscanf(fp, "%s", &line[0]); // Range_step_size
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].m_range_step_size) );
   }

   fscanf(fp, "%s", &line[0]); // Range_start

   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].n_range_start) );
   }

   fscanf(fp, "%s", &line[0]); // Range_end
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].n_range_end) );
   }

   fscanf(fp, "%s", &line[0]); // Range_step_size
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].n_range_step_size) );
   }

   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &num_tests);
   for (i=0; i<num_tests; i++){
      params->lin_solver_paramslist[i].num_tests = num_tests;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].matrix_layout) );
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->lin_solver_paramslist[i].transr = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->lin_solver_paramslist[i].Uplo = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].m) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].n) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].nrhs) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].lda) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].ldb) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].ldab) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].kl) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].ku) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].kd) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->lin_solver_paramslist[i].diag = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->lin_solver_paramslist[i].fact = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->lin_solver_paramslist[i].equed = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->lin_solver_paramslist[i].symm = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &(params->lin_solver_paramslist[i].solver_threhold) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->lin_solver_paramslist[i].equed_porfsx = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].n_err_bnds_porfsx) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].nparams_porfsx) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->lin_solver_paramslist[i].norm_gbcon = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].kl_gbcon) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].ku_gbcon) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->lin_solver_paramslist[i].ldab_gbcon) );
   }

   fclose(fp);

}

/* This function reads parameters needed for Eigen APIs
   from the config settings file 'EIG_PARAMS.dat' and saves in the
   'params->eig_sym_paramslist' structure array   */
void fla_test_read_sym_eig_params( const char *file_name , test_params_t* params )
{
   FILE *fp;
   integer i;
   char line[20];
   char *str, c[20];
   integer num_tests;
   integer mode;
   integer num_ranges;

   str = &c[0];
   fp = fopen( file_name, "r");
   if (fp == NULL){
    printf("Error: Symmetric EIG params config file missing. Exiting.. \n");
    exit(-1);
   }
   /* Read the mode */
   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &mode);
   fscanf(fp, "%*[^\n]\n");

   /* Read the number of Ranges */
   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &num_ranges);

   fscanf(fp, "%s", &line[0]); // Range_start

      for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].n_range_start) );
      params->eig_sym_paramslist[i].num_ranges = num_ranges;
      params->eig_sym_paramslist[i].mode=mode;
   }

   fscanf(fp, "%s", &line[0]); // Range_end
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].n_range_end) );
   }

   fscanf(fp, "%s", &line[0]); // Range_step_size
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].n_range_step_size) );
   }

   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &num_tests);
   for (i=0; i<num_tests; i++){
      params->eig_sym_paramslist[i].num_tests = num_tests;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].matrix_layout) );
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].trans = *str;
   }
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].uplo = *str;
   }
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].job = *str;
   }
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].jobz = *str;
   }
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].vect = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].m) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].n) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].p) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].nrhs) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].lda) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].ldb) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].nb) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].ldt) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].k) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].isgn) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].compz = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].kb) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].itype) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].vect_rd = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].side = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].job_seqr = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].eigsrc = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].initv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].norm = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].diag = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_sym_paramslist[i].storev = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].tsize) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_sym_paramslist[i].threshold_value) );
   }


   fclose(fp);

}

/* This function reads parameters needed for Non symmetric Eigen APIs
   from the config settings file 'EIG_NSYM_PARAMS.dat' and saves in the
   'params->eig_non_sym_paramslist' structure array   */

void fla_test_read_non_sym_eig_params( const char *file_name , test_params_t* params )
{
   FILE *fp;
   integer i;
   char line[20];
   char *str, c[20];
   integer num_tests;
   integer mode;
   integer num_ranges;
   str = c;
   fp = fopen( file_name, "r");
   if (fp == NULL){
    printf("Error: EIG non symmetric API params config file missing. Exiting.. \n");
    exit(-1);
   }
/* Read the mode */
   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &mode);
   fscanf(fp, "%*[^\n]\n");

   /* Read the number of Ranges */
   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &num_ranges);

   fscanf(fp, "%s", &line[0]); // Range_start

      for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->eig_non_sym_paramslist[i].n_range_start) );
      params->eig_non_sym_paramslist[i].num_ranges = num_ranges;
      params->eig_non_sym_paramslist[i].mode=mode;
   }

   fscanf(fp, "%s", &line[0]); // Range_end
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->eig_non_sym_paramslist[i].n_range_end) );
   }

   fscanf(fp, "%s", &line[0]); // Range_step_size
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->eig_non_sym_paramslist[i].n_range_step_size) );
   }

   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &num_tests);
   for (i=0; i<num_tests; i++){
      params->eig_non_sym_paramslist[i].num_tests = num_tests;
   }
   fscanf(fp, "%s", &line[0]);//n
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_non_sym_paramslist[i].n) );
   }
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].howmny = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].initv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].job_seqr = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].eigsrc = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].initv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].job = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].howmny_trsna = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].job_trsen = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].compq = *str;
   }

   /* Reading config params for 'trsyl' API  */
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].trana_real = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].trana_complex = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].tranb_real= *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].tranb_complex = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_non_sym_paramslist[i].isgn) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &(params->eig_non_sym_paramslist[i].gghrd_threshold) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &(params->eig_non_sym_paramslist[i].ggbal_threshold) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &(params->eig_non_sym_paramslist[i].GenNonSymEigProblem_threshold) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].compq_hgeqz = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].compz_hgeqz = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].side_tgevc = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].jobvsl = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].jobvsr = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].sort = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].sense_ggesx = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].balance_ggevx = *str;
   }
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].sense_ggevx = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].sort_gees = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_non_sym_paramslist[i].wantz) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_non_sym_paramslist[i].wantq) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->eig_non_sym_paramslist[i].tgsen_ijob) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->eig_non_sym_paramslist[i].unmhr_trans = *str;
   }


   fclose(fp);
}

/* This function reads parameters needed for SVD APIs
   from the config settings file 'SVD.dat' and saves in the
   'params->svd_paramslist' structure array   */
void fla_test_read_svd_params ( const char *file_name, test_params_t* params )
{
   FILE *fp;
   integer i;
   char line[25];
   char *str, c[20];
   integer num_tests;
   integer mode;
   integer num_ranges;

   str = c;
   fp = fopen( file_name, "r");
   if (fp == NULL){
    printf("Error: SVD config file missing. Exiting.. \n");
    exit(-1);
   }

   /* Read the mode */
   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &mode);
   fscanf(fp, "%*[^\n]\n");

   /* Read the number of Ranges */
   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &num_ranges);

   fscanf(fp, "%s", &line[0]); // Range_start

   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].m_range_start) );
      params->svd_paramslist[i].num_ranges = num_ranges;
      params->svd_paramslist[i].mode = mode;
   }

   fscanf(fp, "%s", &line[0]); // Range_end
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].m_range_end) );
   }

   fscanf(fp, "%s", &line[0]); // Range_step_size
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].m_range_step_size) );
   }

   fscanf(fp, "%s", &line[0]); // Range_start

   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].n_range_start) );
   }

   fscanf(fp, "%s", &line[0]); // Range_end
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].n_range_end) );
   }

   fscanf(fp, "%s", &line[0]); // Range_step_size
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].n_range_step_size) );
   }

   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &num_tests);
   for (i=0; i<num_tests; i++){
      params->svd_paramslist[i].num_tests = num_tests;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].matrix_layout) );
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobu = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobq = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].m) );
   }


   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].p) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].n) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &(params->svd_paramslist[i].tola) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &(params->svd_paramslist[i].tolb) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &(params->svd_paramslist[i].svd_threshold) );
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobu_gesvd = *str;
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobvt_gesvd = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].joba_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobu_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobv_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobr_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobt_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobp_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].m_gejsv) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].n_gejsv) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].joba_gesvj = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobu_gesvj = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobv_gesvj = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].m_gesvj) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].n_gesvj) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].mv_gesvj) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &(params->svd_paramslist[i].ctol_gesvj) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", &(params->svd_paramslist[i].jobu_gesvdx) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", &(params->svd_paramslist[i].jobvt_gesvdx) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", &(params->svd_paramslist[i].range_gesvdx) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].il) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &(params->svd_paramslist[i].iu) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &(params->svd_paramslist[i].vl) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &(params->svd_paramslist[i].vu) );
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].joba_gesvdq = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobu_gesvdq = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      params->svd_paramslist[i].jobv_gesvdq= *str;
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
    int           c, cf;
    char          format_spec[8];
    unsigned int  the_uint;
    int           the_int;
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
                the_uint = va_arg( args, unsigned int );
                fprintf( output_stream, format_spec, the_uint );
                break;

                case 'd':
                the_int = va_arg( args, int );
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
                the_char = va_arg( args, int );
                fprintf( output_stream, "%c", the_char );
                break;
            }

            // Move to next character past type specifier.
            c += 1;
        }
    }
}



void fla_test_parse_command_line( int argc, char** argv )
{
    if ( argc > 1 )
    {
        fprintf( stderr, "Too many command line arguments.\n" );
        exit(1);
    }
    // Copy the binary name to a global string so we can use it later.
    strncpy( fla_test_binary_name, argv[0], MAX_BINARY_NAME_LENGTH );
}



char* fla_test_get_string_for_result( double residual, integer datatype, test_thresh_t* thresh )
{
    char* r_val;

    if ( datatype == FLOAT )
    {
        if      ( residual > thresh->failwarn_s ) r_val = fla_test_fail_string;
        else if ( residual > thresh->warnpass_s ) r_val = fla_test_warn_string;
        else                                      r_val = fla_test_pass_string;
    }
    else if ( datatype == DOUBLE )
    {
        if      ( residual > thresh->failwarn_d ) r_val = fla_test_fail_string;
        else if ( residual > thresh->warnpass_d ) r_val = fla_test_warn_string;
        else                                      r_val = fla_test_pass_string;
    }
    else if ( datatype == COMPLEX )
    {
        if      ( residual > thresh->failwarn_c ) r_val = fla_test_fail_string;
        else if ( residual > thresh->warnpass_c ) r_val = fla_test_warn_string;
        else                                      r_val = fla_test_pass_string;
    }
    else
    {
        if      ( residual > thresh->failwarn_z ) r_val = fla_test_fail_string;
        else if ( residual > thresh->warnpass_z ) r_val = fla_test_warn_string;
        else                                      r_val = fla_test_pass_string;
    }

    return r_val;
}



void fla_test_init_strings( void )
{
    sprintf( fla_test_pass_string, "PASS" );
    sprintf( fla_test_warn_string, "MARGINAL" );
    sprintf( fla_test_fail_string, "FAILURE" );
    sprintf( fla_test_storage_format_string, "Row(r) and General(g) storage format is not supported by External LAPACK interface" );
    sprintf( fla_test_stor_chars, STORAGE_SCHEME_CHARS );
}


void fla_test_op_driver( char*            func_str,
                            char*         impl_var_str,
                            integer       n_pc,
                            char**        pc_str,
                            integer       n_matrices,
                            test_params_t *params,
                            test_thresh_t thresh,
                            integer       api_type,
                            void (*f_exp) (test_params_t *, // params
                                           integer,           // datatype
                                           integer,      // p_cur
                                           integer,  // pci (param combo counter)
                                           integer,  // n_repeats
                                           double*,       // perf
                                           double*,       //time
                                           double* ) )    // residual
{
    integer n_datatypes        = params->n_datatypes;
    integer n_repeats           = params->n_repeats;
    integer p_first, p_max, p_inc;
    integer q_first, q_max, q_inc;
    integer dt, p_cur, q_cur, pci;
    char         datatype_char;
    integer      datatype;
    double       perf, time, residual;
    char*        pass_str;
    char         blank_str[32];
    char         func_param_str[64];
    unsigned int n_spaces;


    fla_test_output_info( "%3sAPI%28s DATA_TYPE%4s SIZE%1s FLOPS%2s TIME(s)%6s ERROR%5s STATUS\n", "", "", "", "", "", "", "" );
    fla_test_output_info( "%3s====%28s==========%4s====%1s=======%2s========%5s==========%2s========\n", "", "", "", "", "", "", "" );

    switch (api_type)
    {
        // Supporting only one range of matrix sizes. So referring only index '0' of paramslist[] array
        case LIN:
           p_first             = params->lin_solver_paramslist[0].m_range_start;
           p_max               = params->lin_solver_paramslist[0].m_range_end;
           p_inc               = params->lin_solver_paramslist[0].m_range_step_size;

           q_first             = params->lin_solver_paramslist[0].n_range_start;
           q_max               = params->lin_solver_paramslist[0].n_range_end;
           q_inc               = params->lin_solver_paramslist[0].n_range_step_size;
           break;

        case EIG_SYM:
           p_first             = params->eig_sym_paramslist[0].m_range_start;
           p_max               = params->eig_sym_paramslist[0].m_range_end;
           p_inc               = params->eig_sym_paramslist[0].m_range_step_size;

           q_first             = params->eig_sym_paramslist[0].n_range_start;
           q_max               = params->eig_sym_paramslist[0].n_range_end;
           q_inc               = params->eig_sym_paramslist[0].n_range_step_size;
           break;

        case EIG_NSYM:
           p_first             = params->eig_non_sym_paramslist[0].m_range_start;
           p_max               = params->eig_non_sym_paramslist[0].m_range_end;
           p_inc               = params->eig_non_sym_paramslist[0].m_range_step_size;

           q_first             = params->eig_non_sym_paramslist[0].n_range_start;
           q_max               = params->eig_non_sym_paramslist[0].n_range_end;
           q_inc               = params->eig_non_sym_paramslist[0].n_range_step_size;
           break;

        case SVD:
           p_first             = params->svd_paramslist[0].m_range_start;
           p_max               = params->svd_paramslist[0].m_range_end;
           p_inc               = params->svd_paramslist[0].m_range_step_size;

           q_first             = params->svd_paramslist[0].n_range_start;
           q_max               = params->svd_paramslist[0].n_range_end;
           q_inc               = params->svd_paramslist[0].n_range_step_size;
           break;

        default: break;
    }


    // Loop over the requested datatypes.
    for ( dt = 0; dt < n_datatypes; ++dt )
    {
        datatype      = params->datatype[dt];
        datatype_char = params->datatype_char[dt];

        // Loop over the requested problem sizes.
        for ( p_cur = p_first, q_cur = q_first; (p_cur <= p_max && q_cur <= q_max); p_cur += p_inc, q_cur += q_inc )
        {
            // Loop over the operation's parameter combinations.
            for ( pci = 0; pci < n_pc; ++pci )
            {
                f_exp( params,
                        datatype,
                        p_cur, pci, n_repeats,
                        &perf, &time, &residual );

                pass_str = fla_test_get_string_for_result( residual, datatype, &thresh );

                // Output the results. Use different formats depending on
                // whether the results are from a front-end or variant.
                fla_test_build_function_string( func_str,
                                                    impl_var_str,
                                                    n_pc, pc_str[pci],
                                                    func_param_str );

                n_spaces = MAX_FUNC_STRING_LENGTH - strlen( func_param_str );
                fill_string_with_n_spaces( blank_str, n_spaces );

                fla_test_output_info( "   %s%s  %c  %12u  %6.3lf  %6.10lf  %9.2le   %s\n",
                                            func_param_str, blank_str,
                                            datatype_char,
                                            p_cur, perf, time, residual, pass_str );
            }
        }

        fla_test_output_info( "\n" );
    }
}




void fla_test_build_function_string( char*        func_base_str,
                                        char*        impl_var_str,
                                        unsigned int n_pc,
                                        char*        pc_str,
                                        char*        func_str )
{

    sprintf( func_str, "%s", func_base_str );


    sprintf( &func_str[strlen(func_str)], "_%s", impl_var_str );

    if ( n_pc > 1 )
        sprintf( &func_str[strlen(func_str)], ":%s", pc_str );


}


void fill_string_with_n_spaces( char* str, unsigned int n_spaces )
{
    unsigned int i;

    for ( i = 0; i < n_spaces; ++i )
        sprintf( &str[i], " " );
}


void fla_test_sleep( void )
{
    int i;

    fla_test_output_info( "Resuming in " );
    for ( i = SECONDS_TO_SLEEP; i > 0; --i )
    {
        fla_test_output_info( "%d ", i );
#ifdef _WIN32
        Sleep(1);
#else
        sleep(1);
#endif
    }
    fla_test_output_info( "\n" );
}


void fla_test_abort( void )
{
    abort();
}


void fla_start_timer( void )
{
    timer = clock();
}


double fla_end_timer( void )
{
    timer = clock() - timer;
    return ((double)timer)/CLOCKS_PER_SEC;
}