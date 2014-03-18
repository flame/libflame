
#include "FLAME.h"

FLA_Error FLA_Asum_external( FLA_Obj x, FLA_Obj asum_x )
{
  FLA_Datatype datatype;
  FLA_Datatype dt_asum;
  int          num_elem;
  int          inc_x;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Asum_check( x, asum_x );

  if ( FLA_Obj_has_zero_dim( x ) )
  {
    FLA_Set( FLA_ZERO, asum_x );
    return FLA_SUCCESS;
  }

  dt_asum  = FLA_Obj_datatype( asum_x );
  datatype = FLA_Obj_datatype( x );

  inc_x    = FLA_Obj_vector_inc( x );
  num_elem = FLA_Obj_vector_dim( x );


  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_x      = ( float * ) FLA_FLOAT_PTR( x );
    float *buff_asum_x = ( float * ) FLA_FLOAT_PTR( asum_x );

    bl1_sasum( num_elem,
               buff_x, inc_x,
               buff_asum_x );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_x      = ( double * ) FLA_DOUBLE_PTR( x );
    double *buff_asum_x = ( double * ) FLA_DOUBLE_PTR( asum_x );

    bl1_dasum( num_elem,
               buff_x, inc_x,
               buff_asum_x );

    break;
  }

  case FLA_COMPLEX:
  {
    if ( dt_asum == FLA_FLOAT )
    {
      scomplex *buff_x      = ( scomplex * ) FLA_COMPLEX_PTR( x );
      float    *buff_asum_x = ( float    * ) FLA_FLOAT_PTR( asum_x );

      bl1_casum( num_elem,
                 buff_x, inc_x,
                 buff_asum_x );
    }
    else if ( dt_asum == FLA_COMPLEX )
    {
      scomplex *buff_x      = ( scomplex * ) FLA_COMPLEX_PTR( x );
      scomplex *buff_asum_x = ( scomplex * ) FLA_COMPLEX_PTR( asum_x );

      bl1_casum( num_elem,
                 buff_x, inc_x,
                 &(buff_asum_x->real) );
      buff_asum_x->imag = 0.0F;
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    if ( dt_asum == FLA_DOUBLE )
    {
      dcomplex *buff_x      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );
      double   *buff_asum_x = ( double   * ) FLA_DOUBLE_PTR( asum_x );

      bl1_zasum( num_elem,
                 buff_x, inc_x,
                 buff_asum_x );
    }
    else if ( dt_asum == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_x      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );
      dcomplex *buff_asum_x = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( asum_x );

      bl1_zasum( num_elem,
                 buff_x, inc_x,
                 &(buff_asum_x->real) );
      buff_asum_x->imag = 0.0;
    }

    break;
  }

  }

  return FLA_SUCCESS;
}

