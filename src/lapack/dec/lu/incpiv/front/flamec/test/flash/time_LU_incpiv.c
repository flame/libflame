#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_LU_incpiv( FLA_Obj A, FLA_Obj p );
void time_LU(
              int is_pivoting, int type, int nrepeats, int m, int n, dim_t nb_alg, dim_t nb_flash,
              FLA_Obj A, FLA_Obj p, FLA_Obj b, FLA_Obj b_ref, FLA_Obj b_norm,
              double *dtime, double *diff, double *gflops );


void time_LU(
              int is_pivoting, int type, int nrepeats, int m, int n, dim_t nb_alg, dim_t nb_flash,
              FLA_Obj A, FLA_Obj p, FLA_Obj b, FLA_Obj b_ref, FLA_Obj b_norm,
              double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj A_save, AH_save;
  FLA_Obj AH, pH, bH, LH;

  //FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  //FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, A, &A_flat );

  FLASH_LU_incpiv_create_hier_matrices( A, 1, &nb_flash, nb_alg,
                                        &AH, &pH, &LH );
  FLASH_Obj_create_hier_copy_of_flat( b, 1, &nb_flash, &bH );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, AH, &AH_save );

  FLA_Copy_external( A, A_save );
  FLASH_Copy( AH, AH_save );


  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLA_Copy_external( A_save, A );
    FLASH_Copy( AH_save, AH );

    *dtime = FLA_Clock();

    switch( is_pivoting ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_LU_incpiv( A, p );
        break;
      case FLA_ALG_FRONT:
        FLASH_LU_incpiv( AH, pH, LH );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );
  }

  if ( type == FLA_ALG_REFERENCE )
  {
    *diff = 0.0;
  }
  else
  {
    FLASH_FS_incpiv( AH, pH, LH, bH );
    FLASH_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                AH, bH );

    FLASH_Obj_flatten( bH, b );

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE,
                       A_save, b, FLA_MINUS_ONE, b_ref );

    FLA_Nrm2_external( b_ref, b_norm );
    FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, b_norm,
                               1, 1, diff, 1, 1 );
  }

  *gflops = 2.0 / 3.0 * m * m * n /
            dtime_old / 1e9;

  *dtime = dtime_old;

  FLASH_Copy( A_save, A );

  FLA_Obj_free( &A_save );
  FLASH_Obj_free( &AH_save );

  FLASH_Obj_free( &AH );
  FLASH_Obj_free( &pH );
  FLASH_Obj_free( &bH );
  FLASH_Obj_free( &LH );
}

