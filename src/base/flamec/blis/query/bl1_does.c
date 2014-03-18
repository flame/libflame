
#include "blis1.h"

int bl1_does_trans( trans1_t trans )
{
	return ( trans == BLIS1_TRANSPOSE ||
	         trans == BLIS1_CONJ_TRANSPOSE );
}

int bl1_does_notrans( trans1_t trans )
{
	return ( trans == BLIS1_NO_TRANSPOSE ||
	         trans == BLIS1_CONJ_NO_TRANSPOSE );
}

int bl1_does_conj( trans1_t trans )
{
	return ( trans == BLIS1_CONJ_NO_TRANSPOSE ||
	         trans == BLIS1_CONJ_TRANSPOSE );
}

