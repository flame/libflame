
#include "FLAME.h"
#include "FLA_Syrk_ln_omp.h"

FLA_Error FLA_Syrk_ln_omp1t_var5_rc( FLA_Obj A, FLA_Obj C, int nb_alg )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;
  FLA_Obj MyC;

  int b;
  
  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  #pragma intel omp parallel taskq
  {
  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    nb_alg = FLA_Obj_width( A )/omp_get_num_threads() + 1;
    b = min( FLA_Obj_width( AR ), nb_alg );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    #pragma intel omp task captureprivate(A1) private(MyC)
    {
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &MyC );
    FLA_Obj_set_to_zero( MyC );
    
    /* MyC := A1 * A1' */
    FLA_Syrk( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_ONE, A1, FLA_ZERO, MyC );

    /* C := MyC */
    //REF_Axpy_sync_circular( FLA_ONE, MyC, C );

    FLA_Obj_free( &MyC );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );
  }
  }

  return FLA_SUCCESS;
}

