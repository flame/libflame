
#include "FLAME.h"

float FLA_random_float()
{
  return ( float )( ( ( double ) rand() / ( double ) RAND_MAX ) * 2.0 - 1.0 );
}

double FLA_random_double()
{
  return ( ( double ) rand() / ( double ) RAND_MAX ) * 2.0 - 1.0;
}

scomplex FLA_random_scomplex()
{
  scomplex z;

  z.real = FLA_random_float();
  z.imag = FLA_random_float();

  return z;
}

dcomplex FLA_random_dcomplex()
{
  dcomplex z;

  z.real = FLA_random_double();
  z.imag = FLA_random_double();

  return z;
}

