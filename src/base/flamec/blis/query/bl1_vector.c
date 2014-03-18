
#include "blis1.h"

int bl1_vector_dim( int m, int n )
{
	if ( m == 1 ) return n;
	else          return m;
}

int bl1_vector_inc( trans1_t trans, int m, int n, int rs, int cs )
{
	if ( bl1_does_notrans( trans ) )
	{
		if ( m == 1 ) return cs;
		else          return rs;
	}
	else // if ( bl1_does_trans( trans ) )
	{
		if ( m == 1 ) return rs;
		else          return cs;
	}
}
