
#include "FLAME.h"

FLA_Error FLA_SA_LU_task( FLA_Obj U,
                          FLA_Obj D, FLA_Obj p, FLA_Obj L, dim_t nb_alg,
                                                           fla_lu_t* cntl )
{
  FLA_Error info;

  info = FLA_SA_LU_blk( U,
                        D, p, L, nb_alg );

  return info;
}

