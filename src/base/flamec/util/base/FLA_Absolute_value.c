
#include "FLAME.h"

FLA_Error FLA_Absolute_value( FLA_Obj alpha )
{
  FLA_Datatype datatype;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Absolute_value_check( alpha );

  datatype = FLA_Obj_datatype( alpha );
  
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    *buff_alpha = ( float ) fabs( ( double ) *buff_alpha );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    *buff_alpha = fabs( *buff_alpha );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );

    buff_alpha->real = ( float ) sqrt( ( double ) buff_alpha->real * buff_alpha->real + 
                                                  buff_alpha->imag * buff_alpha->imag );
    buff_alpha->imag = 0.0F; 


    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );

    buff_alpha->real = sqrt( buff_alpha->real * buff_alpha->real + 
                             buff_alpha->imag * buff_alpha->imag );
    buff_alpha->imag = 0.0; 

    break;
  }

  }

  return FLA_SUCCESS;
}

