
#include "FLAME.h"

FLA_Error FLA_Mult_add( FLA_Obj alpha, FLA_Obj beta, FLA_Obj gamma )
{
  FLA_Datatype datatype;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Mult_add_check( alpha, beta, gamma );

  datatype = FLA_Obj_datatype( gamma );
  
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );
    float *buff_gamma = ( float * ) FLA_FLOAT_PTR( gamma );

    *buff_gamma = *buff_gamma + *buff_alpha * *buff_beta;

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );
    double *buff_gamma = ( double * ) FLA_DOUBLE_PTR( gamma );

    *buff_gamma = *buff_gamma + *buff_alpha * *buff_beta;

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );
    scomplex *buff_beta  = ( scomplex * ) FLA_COMPLEX_PTR( beta );
    scomplex *buff_gamma = ( scomplex * ) FLA_COMPLEX_PTR( gamma );
    scomplex  alphabeta;

    alphabeta.real = buff_alpha->real * buff_beta->real -
                     buff_alpha->imag * buff_beta->imag;

    alphabeta.imag = buff_alpha->real * buff_beta->imag +
                     buff_alpha->imag * buff_beta->real;

    buff_gamma->real = buff_gamma->real + alphabeta.real;
    buff_gamma->imag = buff_gamma->imag + alphabeta.imag;

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    dcomplex *buff_beta  = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );
    dcomplex *buff_gamma = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( gamma );
    dcomplex  alphabeta;

    alphabeta.real = buff_alpha->real * buff_beta->real -
                     buff_alpha->imag * buff_beta->imag;

    alphabeta.imag = buff_alpha->real * buff_beta->imag +
                     buff_alpha->imag * buff_beta->real;

    buff_gamma->real = buff_gamma->real + alphabeta.real;
    buff_gamma->imag = buff_gamma->imag + alphabeta.imag;

    break;
  }

  }

  return FLA_SUCCESS;
}

