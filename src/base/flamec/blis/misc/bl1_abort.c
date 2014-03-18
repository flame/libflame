
#include "blis1.h"

void bl1_abort( void )
{
	abort();
}

void bl1_abort_msg( char* message )
{
	fprintf( stderr, "BLIS: %s\n", message );
	fprintf( stderr, "BLIS: Aborting.\n" );
	bl1_abort();
}

