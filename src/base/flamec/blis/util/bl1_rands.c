
#include "blis1.h"

void bl1_srands( float* alpha )
{
	*alpha = ( float ) ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;
}

void bl1_drands( double* alpha )
{
	*alpha = ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0 ) ) - 1.0;
}

void bl1_crands( scomplex* alpha )
{
	bl1_srands( &(alpha->real) );
	bl1_srands( &(alpha->imag) );
}

void bl1_zrands( dcomplex* alpha )
{
	bl1_drands( &(alpha->real) );
	bl1_drands( &(alpha->imag) );
}

