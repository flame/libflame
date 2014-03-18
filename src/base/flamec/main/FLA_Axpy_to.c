
#include "FLAME.h"

FLA_Error FLA_Axpy_buffer_to_object( FLA_Trans trans, FLA_Obj alpha, dim_t m, dim_t n, void* X_buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj Y )
{
  FLA_Obj  X;
  FLA_Obj  YTL, YTR, 
           YBL, Yij;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Axpy_buffer_to_object_check( trans, alpha, m, n, X_buffer, rs, cs, i, j, Y );

  FLA_Part_2x2( Y,  &YTL, &YTR,
                    &YBL, &Yij,     i, j, FLA_TL );

  FLA_Obj_create_without_buffer( FLA_Obj_datatype( Y ), m, n, &X );
  FLA_Obj_attach_buffer( X_buffer, rs, cs, &X );

  FLA_Axpyt_external( trans, alpha, X, Yij );

  FLA_Obj_free_without_buffer( &X );

  return FLA_SUCCESS;
}



FLA_Error FLA_Axpy_object_to_buffer( FLA_Trans trans, FLA_Obj alpha, dim_t i, dim_t j, FLA_Obj X, dim_t m, dim_t n, void* Y_buffer, dim_t rs, dim_t cs )
{
  FLA_Obj  Y;
  FLA_Obj  XTL, XTR, 
           XBL, Xij;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Axpy_object_to_buffer_check( trans, alpha, i, j, X, m, n, Y_buffer, rs, cs );

  FLA_Part_2x2( X,  &XTL, &XTR,
                    &XBL, &Xij,     i, j, FLA_TL );

  FLA_Obj_create_without_buffer( FLA_Obj_datatype( X ), m, n, &Y );
  FLA_Obj_attach_buffer( Y_buffer, rs, cs, &Y );

  FLA_Axpyt_external( trans, alpha, Xij, Y );

  FLA_Obj_free_without_buffer( &Y );

  return FLA_SUCCESS;
}

