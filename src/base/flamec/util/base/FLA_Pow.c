
#include "FLAME.h"

FLA_Error FLA_Pow( FLA_Obj base, FLA_Obj exp, FLA_Obj btoe )
{
  FLA_Datatype datatype;
  int          r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Pow_check( base, exp, btoe );

  datatype = FLA_Obj_datatype( base );
  
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_base = ( float * ) FLA_FLOAT_PTR( base );
    float *buff_exp  = ( float * ) FLA_FLOAT_PTR( exp );
    float *buff_btoe = ( float * ) FLA_FLOAT_PTR( btoe );

    *buff_btoe = ( float ) pow( *buff_base, *buff_exp );
    
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_base = ( double * ) FLA_DOUBLE_PTR( base );
    double *buff_exp  = ( double * ) FLA_DOUBLE_PTR( exp );
    double *buff_btoe = ( double * ) FLA_DOUBLE_PTR( btoe );

    *buff_btoe = ( double ) pow( *buff_base, *buff_exp );
    
    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_base = ( scomplex * ) FLA_COMPLEX_PTR( base );
    scomplex *buff_exp  = ( scomplex * ) FLA_COMPLEX_PTR( exp );
    scomplex *buff_btoe = ( scomplex * ) FLA_COMPLEX_PTR( btoe );

    buff_btoe->real = ( float ) pow( buff_base->real, buff_exp->real );
    buff_btoe->imag = 0.0;
    
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_base = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( base );
    dcomplex *buff_exp  = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( exp );
    dcomplex *buff_btoe = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( btoe );

    buff_btoe->real = ( double ) pow( buff_base->real, buff_exp->real );
    buff_btoe->imag = 0.0;
    
    break;
  }

  }

  return r_val;
}

