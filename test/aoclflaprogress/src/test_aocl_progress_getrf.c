#include "test_aocl_progress.h"
//#define TEST_SGETRF
//#define TEST_DGETRF
//#define TEST_CGETRF
//#define TEST_ZGETRF
#define TEST_SGBTRF
#define TEST_DGBTRF
#define TEST_CGBTRF
#define TEST_ZGBTRF

#define PROCESS_ENABLE
#ifdef PROCESS_ENABLE

int aocl_fla_progress(const char* const api,const integer lenapi,const  integer* const progress,const integer* const current_thread,const integer* const total_threads)
{
  printf( "In AOCL FLA  Progress thread  %"FS", at API  %s, progress  %"FS" total threads= %"FS"\n", *current_thread, api, *progress,*total_threads );
  return 0;

}
#endif
// uncomment to enable AOCL_FLA_SET_Progress Printing
#define SET_PROCESS_ENABLE
#ifdef SET_PROCESS_ENABLE

int test_progress(const char* const api,const integer lenapi,const integer * const progress,const integer *const current_thread,const integer *const total_threads)
{
  printf( "In AOCL Progress thread  %"FS", at API  %s, progress %"FS" total threads= %"FS" \n", *current_thread, api, *progress,*total_threads );
  return 0;

}
#endif

int main()
{
    integer  m, n,lda,kl,ku,loop_count,info;
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

	// Read the kl
        fla_test_read_next_line( buffer, input_stream );
        sscanf( buffer, "%"FS"", &kl);

        // Read the ku
        fla_test_read_next_line( buffer, input_stream );
        sscanf( buffer, "%"FS"", &ku);

        // Read the loop_count
        fla_test_read_next_line( buffer, input_stream );
        sscanf( buffer, "%"FS"", &loop_count);

	#ifdef SET_PROCESS_ENABLE
                aocl_fla_set_progress(test_progress);
        #endif


     	#ifdef TEST_SGBTRF
    	for(integer i =0;i<loop_count;++i){
    		A = (float *)malloc(m*n*sizeof(float));
    		ipiv = (integer *)malloc(m*sizeof(integer));
    		/* Initialization of input Buffers */
    		printf("Testing SGBTRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
        	fla_test_init_buffer_rand_s( A,m,n,lda);
        	sgbtrf_(&m,&n,&kl,&ku,A,&lda,ipiv,&info);
        	free(A);
        	free(ipiv);
    	}
        #endif

        #ifdef TEST_DGBTRF
 	for(integer i =0;i<loop_count;++i){
        	A = (double *)malloc(m*n*sizeof(double));
    		ipiv = (integer *)malloc(m*sizeof(integer));
    		/* Initialization of input Buffers */
    		printf("Testing DGBTRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
        	fla_test_init_buffer_rand_d( A,m,n,lda);
       		dgbtrf_(&m,&n,&kl,&ku,A,&lda,ipiv,&info);
        	free(A);
        	free(ipiv);
	}
        #endif

    	#ifdef TEST_CGBTRF
	for(integer i =0;i<loop_count;++i){
    		A = (scomplex *)malloc(m*n*sizeof(scomplex));
    		ipiv = (integer *)malloc(m*sizeof(integer));
    		/* Initialization of input Buffers */
    		printf("Testing CGBTRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
        	fla_test_init_buffer_rand_c( A,m,n,lda);
        	cgbtrf_(&m,&n,&kl,&ku,A,&lda,ipiv,&info);
        	free(A);
        	free(ipiv);
	}
        #endif

	#ifdef TEST_ZGBTRF
        for(integer i =0;i<loop_count;++i){
                A = (dcomplex *)malloc(m*n*sizeof(dcomplex));
                ipiv = (integer *)malloc(m*sizeof(integer));
                /* Initialization of input Buffers */
                printf("Testing ZGBTRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
                fla_test_init_buffer_rand_z( A,m,n,lda);
                zgbtrf_(&m,&n,&kl,&ku,A,&lda,ipiv,&info);
                free(A);
                free(ipiv);
        }
        #endif



	#ifdef TEST_SGETRF
	for(integer i =0;i<loop_count;++i){
		A = (float *)malloc(m*n*sizeof(float));
    		ipiv = (integer *)malloc(m*sizeof(integer));
    		/* Initialization of input Buffers */
    		printf("Testing SGETRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
		fla_test_init_buffer_rand_s( A,m,n,lda);
		sgetrf_(&m,&n,A,&lda,ipiv,&info);
		free(A);
		free(ipiv);
	}
	#endif

	#ifdef TEST_DGETRF
        for(integer i =0;i<loop_count;++i){
                A = (double *)malloc(m*n*sizeof(double));
                ipiv = (integer *)malloc(m*sizeof(integer));
                /* Initialization of input Buffers */
                printf("Testing DGETRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
                fla_test_init_buffer_rand_d( A,m,n,lda);
                dgetrf_(&m,&n,A,&lda,ipiv,&info);
                free(A);
                free(ipiv);
        }
        #endif

	#ifdef TEST_CGETRF
        for(integer i =0;i<loop_count;++i){
                A = (scomplex *)malloc(m*n*sizeof(scomplex));
		ipiv = (integer *)malloc(m*sizeof(integer));
                /* Initialization of input Buffers */
                printf("Testing CGETRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
                fla_test_init_buffer_rand_c( A,m,n,lda);
                cgetrf_(&m,&n,A,&lda,ipiv,&info);
                free(A);
                free(ipiv);
        }
        #endif

	#ifdef TEST_ZGETRF
        for(integer i =0;i<loop_count;++i){
                A = (dcomplex *)malloc(m*n*sizeof(dcomplex));
                ipiv = (integer *)malloc(m*sizeof(integer));
                /* Initialization of input Buffers */
                printf("Testing ZGETRF API values are m = %"FS" ,n=%"FS",lda=%"FS"  \n",m,n,lda);
                fla_test_init_buffer_rand_z( A,m,n,lda);
                zgetrf_(&m,&n,A,&lda,ipiv,&info);
                free(A);
                free(ipiv);
        }
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
