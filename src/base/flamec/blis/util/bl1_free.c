
#include "blis1.h"

#ifdef BLIS1_ENABLE_USE_OF_FLA_MALLOC
  #include "FLAME.h"
  #define BLIS1_FREE FLA_free
#else
  #define BLIS1_FREE free
#endif

void bl1_vfree( void* p )
{
	free( ( void* ) p );
}

void bl1_ifree( int* p )
{
	free( ( int* ) p );
}

void bl1_sfree( float* p )
{
	free( ( void* ) p );
}

void bl1_dfree( double* p )
{
	free( ( void* ) p );
}

void bl1_cfree( scomplex* p )
{
	free( ( void* ) p );
}

void bl1_zfree( dcomplex* p )
{
	free( ( void* ) p );
}

