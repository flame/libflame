
#include "blis1.h"

void bl1_set_dims_with_trans( trans1_t trans, int m, int n, int* m_new, int* n_new )
{
	if ( bl1_does_trans( trans ) )
	{
		*m_new = n;
		*n_new = m;
	}
	else
	{
		*m_new = m;
		*n_new = n;
	}
}

void bl1_set_dim_with_side( side1_t side, int m, int n, int* dim_new )
{
	if ( bl1_is_left( side ) )
	{
		*dim_new = m;
	}
	else // if ( bl1_is_right( side ) )
	{
		*dim_new = n;
	}
}

