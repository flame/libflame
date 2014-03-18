
#include "blis1.h"

void bl1_sinverts( conj1_t conj, float* alpha )
{
	float  one = 1.0F;

	*alpha = one / *alpha;
}

void bl1_dinverts( conj1_t conj, double* alpha )
{
	double one = 1.0;

	*alpha = one / *alpha;
}

void bl1_cinverts( conj1_t conj, scomplex* alpha )
{
	float  temp;
	float  s, xr_s, xi_s;

	s           = bl1_fmaxabs( alpha->real, alpha->imag ); \
	xr_s        = alpha->real / s;
	xi_s        = alpha->imag / s;
	temp        = xr_s * alpha->real + xi_s * alpha->imag;

	alpha->real =  xr_s / temp;
	alpha->imag = -xi_s / temp;

	if ( bl1_is_conj( conj ) )
		bl1_cconjs( alpha );
}

void bl1_zinverts( conj1_t conj, dcomplex* alpha )
{
	double temp;
	double s, xr_s, xi_s;

	s           = bl1_fmaxabs( alpha->real, alpha->imag ); \
	xr_s        = alpha->real / s;
	xi_s        = alpha->imag / s;
	temp        = xr_s * alpha->real + xi_s * alpha->imag;

	alpha->real =  xr_s / temp;
	alpha->imag = -xi_s / temp;

	if ( bl1_is_conj( conj ) )
		bl1_zconjs( alpha );
}

