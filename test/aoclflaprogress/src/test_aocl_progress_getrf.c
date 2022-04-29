#include "test_aocl_progress.h"

#define PROCESS_ENABLE
#ifdef PROCESS_ENABLE
#define BUFLEN 16
int aocl_fla_progress(char* api,integer lenapi,integer *progress,integer *current_thread,integer *total_threads)
{
  char buf[BUFLEN];
  if( lenapi >= BUFLEN ) lenapi = BUFLEN-1;
  strncpy( buf, api, lenapi );
  buf[lenapi] = '\0';
  printf( "In AOCL FLA  Progress thread  %"FS", at API  %s, progress  %"FS" total threads= %"FS"\n", *current_thread, buf, *progress,*total_threads );
  return 0;

}
#endif
// uncomment to enable AOCL_FLA_SET_Progress Printing
#define SET_PROCESS_ENABLE
#ifdef SET_PROCESS_ENABLE
#define BUFLEN 16
int test_progress(char* api,integer lenapi,integer *progress,integer *current_thread,integer *total_threads)
{
char buf[BUFLEN];
  if( lenapi >= BUFLEN ) lenapi = BUFLEN-1;
  strncpy( buf, api, lenapi );
  buf[lenapi] = '\0';
  printf( "In AOCL Progress thread  %"FS", at API  %s, progress %"FS" total threads= %"FS" \n", *current_thread, buf, *progress,*total_threads );
  return 0;

}
#endif

int main()
{
    integer  m, n,lda,info;
    void *A;
    void  *ipiv;
    FILE* input_stream;
	char  buffer[ INPUT_BUFFER_SIZE ];
	char *input_filename= "input.global.general";
	// Attempt to open input file corresponding to input_filename as
	// read-only/binary.
	input_stream = fopen(input_filename, "rb" );

	// Check for success.
	if ( input_stream == NULL )
	{
		printf( "Failed to open input file %s. Check existence and permissions.\n", input_filename );
	}

	// Read the m
	fla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%"FS"", &m);

	// Read the n
	fla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%"FS"", &n);
    
    lda = m;
    A = (double *)malloc(m*n*sizeof(double));
    ipiv = (integer *)malloc(m*sizeof(integer));
    /* Initialization of input Buffers */
    fla_test_init_buffer_rand( A, m*n );
    printf("in main m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
    #ifdef SET_PROCESS_ENABLE
      	aocl_fla_set_progress(test_progress);
    #endif 
	dgetrf_(&m,&n,A,&lda,ipiv,&info);

    return 0;

}

void  fla_test_init_buffer_rand( double *buf1, integer size)
{
   integer j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
    }
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
				printf( "Error reading input file: encountered unexpected EOF.\n " );
			else
				printf( "Error (non-EOF) reading input file. \n" );
		}
	}
	// We continue to read lines into buffer until the line is neither
	// commented nor blank.
	while ( temp[0] == COMMENT_CHAR || temp[0] == '\n' ||
			temp[0] == ' '          || temp[0] == '\t' );

	// Save the string in temp, up to first white space character, into buffer.
	sscanf( temp, "%s ", buffer );
}
