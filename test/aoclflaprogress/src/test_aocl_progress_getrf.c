#include "test_aocl_progress.h"
#define TEST_SGETRF
#define TEST_DGETRF
#define TEST_CGETRF
#define TEST_ZGETRF
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

	#ifdef TEST_SGETRF
    A = (float *)malloc(m*n*sizeof(float));
    ipiv = (integer *)malloc(m*sizeof(integer));
    /* Initialization of input Buffers */
    printf("Testing SGETRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
	fla_test_init_buffer_rand_s( A,m,n,lda);
    #ifdef SET_PROCESS_ENABLE
      	aocl_fla_set_progress(test_progress);
    #endif 
	sgetrf_(&m,&n,A,&lda,ipiv,&info);
	free(A);
	free(ipiv);
	#endif
	#ifdef TEST_DGETRF
    A = (double *)malloc(m*n*sizeof(double));
    ipiv = (integer *)malloc(m*sizeof(integer));
    /* Initialization of input Buffers */
    printf("Testing DGETRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
	fla_test_init_buffer_rand_d( A,m,n,lda);
    #ifdef SET_PROCESS_ENABLE
      	aocl_fla_set_progress(test_progress);
    #endif 
	dgetrf_(&m,&n,A,&lda,ipiv,&info);
	free(A);
	free(ipiv);
	#endif

    #ifdef TEST_CGETRF
    A = (scomplex *)malloc(m*n*sizeof(scomplex));
    ipiv = (integer *)malloc(m*sizeof(integer));
    /* Initialization of input Buffers */
    printf("Testing CGETRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
	fla_test_init_buffer_rand_c( A,m,n,lda);
    #ifdef SET_PROCESS_ENABLE
      	aocl_fla_set_progress(test_progress);
    #endif 
	cgetrf_(&m,&n,A,&lda,ipiv,&info);
	free(A);
	free(ipiv);
	#endif
	#ifdef TEST_ZGETRF
    A = (dcomplex *)malloc(m*n*sizeof(dcomplex));
    ipiv = (integer *)malloc(m*sizeof(integer));
    /* Initialization of input Buffers */
    printf("Testing ZGETRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
	fla_test_init_buffer_rand_z( A,m,n,lda);
    #ifdef SET_PROCESS_ENABLE
      	aocl_fla_set_progress(test_progress);
    #endif 
	zgetrf_(&m,&n,A,&lda,ipiv,&info);
	free(A);
	free(ipiv);
	#endif

    return 0;

}

void  fla_test_init_buffer_rand_d( double *buf1, integer M,integer N,integer LDA)
{
   int i, j;

   for( i = 0; i < N; i++ )
   {
      for( j = 0; j < M; j++ )
      {
         buf1[i * LDA + j] = DRAND();
      }
   }

   return;

}

void  fla_test_init_buffer_rand_z( dcomplex *buf1, integer M,integer N,integer LDA)
{
   
   int i, j;

   for( i = 0; i < N; i++ )
   {
      for( j = 0; j < M; j++ )
      {
         buf1[i * LDA + j].real = DRAND();
         buf1[i * LDA + j].imag = DRAND();
      }
   }

   return;
}

void  fla_test_init_buffer_rand_s( float *buf1, integer M,integer N,integer LDA)
{
   int i, j;

   for( i = 0; i < N; i++ )
   {
      for( j = 0; j < M; j++ )
      {
         buf1[i * LDA + j] = SRAND();
      }
   }

   return;

}

void  fla_test_init_buffer_rand_c( scomplex *buf1, integer M,integer N,integer LDA)
{
   
   int i, j;

   for( i = 0; i < N; i++ )
   {
      for( j = 0; j < M; j++ )
      {
         buf1[i * LDA + j].real = SRAND();
         buf1[i * LDA + j].imag = SRAND();
      }
   }

   return;
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
