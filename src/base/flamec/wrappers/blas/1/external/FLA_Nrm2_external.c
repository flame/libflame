
#include "FLAME.h"

FLA_Error FLA_Nrm2_external( FLA_Obj x, FLA_Obj norm_x )
{
  FLA_Datatype datatype;
  int          num_elem;
  int          inc_x;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Nrm2_check( x, norm_x );

  if ( FLA_Obj_has_zero_dim( x ) )
  {
    FLA_Set( FLA_ZERO, norm_x );
    return FLA_SUCCESS;
  }

  datatype = FLA_Obj_datatype( x );

  inc_x    = FLA_Obj_vector_inc( x );
  num_elem = FLA_Obj_vector_dim( x );


  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_x      = ( float * ) FLA_FLOAT_PTR( x );
    float *buff_norm_x = ( float * ) FLA_FLOAT_PTR( norm_x );

    bl1_snrm2( num_elem,
               buff_x, inc_x,
               buff_norm_x );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_x      = ( double * ) FLA_DOUBLE_PTR( x );
    double *buff_norm_x = ( double * ) FLA_DOUBLE_PTR( norm_x );

    bl1_dnrm2( num_elem,
               buff_x, inc_x,
               buff_norm_x );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_x      = ( scomplex * ) FLA_COMPLEX_PTR( x );
    float    *buff_norm_x = ( float    * ) FLA_FLOAT_PTR( norm_x );

    bl1_cnrm2( num_elem,
               buff_x, inc_x,
               buff_norm_x );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_x      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );
    double   *buff_norm_x = ( double   * ) FLA_DOUBLE_PTR( norm_x );

    bl1_znrm2( num_elem,
               buff_x, inc_x,
               buff_norm_x );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

