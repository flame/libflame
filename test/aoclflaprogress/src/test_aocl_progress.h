#include "FLAME.h"
#include <string.h>

#define INPUT_BUFFER_SIZE  256
#define COMMENT_CHAR       '#'

#if defined(FLA_ENABLE_ILP64)
#define FS "llu"
#else
#define FS "lu"
#endif

void fla_test_read_next_line( char* buffer, FILE* input_stream );
void  fla_test_init_buffer_rand( double *buf1, integer size);
