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

	// Read the main test suite parameters.
	fla_test_read_parameter_file( PARAMETERS_FILENAME, &params );

	// Test the LAPACK-level operations.
	fla_test_lapack_suite( OPERATIONS_FILENAME, params, ops );

	return 0;
}


void fla_test_lapack_suite( char* input_filename, test_params_t params, test_ops_t ops )
{
	char buffer[ INPUT_BUFFER_SIZE ];
	int op;
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

	while(fla_test_read_tests_for_op( input_stream, &op, buffer))
	{
		if(op)
		{
			for(int i=0; i < test_api_count; i++)
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


void fla_test_read_parameter_file( char* input_filename, test_params_t* params )
{
	FILE* input_stream;
	char  buffer[ INPUT_BUFFER_SIZE ];
	char  temp[ INPUT_BUFFER_SIZE ];
	int   i;

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

	// Read the initial problem size to test.
	fla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%" FS , &(params->p_first) );

	// Read the maximum problem size to test.
	fla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%" FS , &(params->p_max) );

	// Read the problem size increment to test.
	fla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%" FS , &(params->p_inc) );

	// Read the partial number of matrix size for incomplete factorization.
	fla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &(params->p_nfact) );


	// Close the file.
	fclose( input_stream );

	fla_test_output_info( "\n" );
	fla_test_output_info( "--- test suite parameters ----------------------------\n" );
	fla_test_output_info( "\n" );
	fla_test_output_info( "n_repeats            %u\n", params->n_repeats );
	fla_test_output_info( "n_datatypes          %u\n", params->n_datatypes );
	fla_test_output_info( "datatype[0]          %d (%c)\n", params->datatype[0], params->datatype_char[0] );
	for( i = 1; i < params->n_datatypes; ++i )
		fla_test_output_info( "        [%d]          %d (%c)\n", i, params->datatype[i], params->datatype_char[i] );
	fla_test_output_info( "p_first              %u\n", params->p_first );
	fla_test_output_info( "p_max                %u\n", params->p_max );
	fla_test_output_info( "p_inc                %u\n", params->p_inc );
	fla_test_output_info( "p_nfact              %d\n", params->p_nfact );
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
							test_params_t params,
							test_thresh_t thresh,
							void (*f_exp) (test_params_t, // params
										   int,           // datatype
										   integer,      // p_cur
										   integer,  // pci (param combo counter)
										   integer,  // n_repeats
										   double*,       // perf
										   double*,       //time
										   double* ) )    // residual
{
	integer n_datatypes        = params.n_datatypes;
	integer p_first             = params.p_first;
	integer p_max               = params.p_max;
	integer p_inc               = params.p_inc;
	integer n_repeats           = params.n_repeats;
	integer dt, p_cur, pci;
	char         datatype_char;
	integer      datatype;
	double       perf, time, residual;
	char*        pass_str;
	char         blank_str[32];
	char         func_param_str[64];
	unsigned int n_spaces;


	fla_test_output_info( "%3sAPI%28s DATA_TYPE%4s SIZE%1s FLOPS%2s TIME(s)%6s ERROR%5s STATUS\n", "", "", "", "", "", "", "" );
	fla_test_output_info( "%3s====%28s==========%4s====%1s=======%2s========%5s==========%2s========\n", "", "", "", "", "", "", "" );


	// Loop over the requested datatypes.
	for ( dt = 0; dt < n_datatypes; ++dt )
	{
		datatype      = params.datatype[dt];
		datatype_char = params.datatype_char[dt];

		// Loop over the requested problem sizes.
		for ( p_cur = p_first; p_cur <= p_max; p_cur += p_inc )
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