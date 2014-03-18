
#include "FLAME.h"

FLA_Error REF_Svd_uv( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
                      double* dtime_bred, double* dtime_bsvd, double* dtime_appq,
                      double* dtime_qrfa, double* dtime_gemm )

{
  *dtime_bred = 1;
  *dtime_bsvd = 1;
  *dtime_appq = 1;
  *dtime_qrfa = 1;
  *dtime_gemm = 1;

  return FLA_Svd_external( FLA_SVD_VECTORS_ALL, FLA_SVD_VECTORS_ALL, A, s, U, V );
}
