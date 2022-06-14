#include "FLAME.h"
#include <string.h>

#define INPUT_BUFFER_SIZE  256
#define COMMENT_CHAR       '#'

#if defined(FLA_ENABLE_ILP64)
#define FS "llu"
#else
#define FS "lu"
#endif

#define DRAND()  ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;
#define SRAND()  ( float ) ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;

void  fla_test_read_next_line( char* buffer, FILE* input_stream );
void  fla_test_init_buffer_rand_s( float *buf1, integer M,integer N,integer LDA);
void  fla_test_init_buffer_rand_d( double *buf1, integer M,integer N,integer LDA);
void  fla_test_init_buffer_rand_c( scomplex *buf1, integer M,integer N,integer LDA);
void  fla_test_init_buffer_rand_z( dcomplex *buf1, integer M,integer N,integer LDA);
