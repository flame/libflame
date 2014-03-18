
#include "FLAME.h"

void FLA_Gemm_pack_C( FLA_Trans transC, FLA_Obj C, FLA_Obj *packed_C )
{
  *packed_C = C;
}

void FLA_Gemm_unpack_andor_scale_C( FLA_Trans transC, FLA_Obj alpha, 
				    FLA_Obj C, FLA_Obj *packed_C )
{
}

void FLA_Gemm_pack_andor_scale_B( FLA_Trans transB, FLA_Obj alpha,
				  FLA_Obj B, FLA_Obj *packed_B )
{
  *packed_B = B;
}

void FLA_Gemm_release_pack_B( FLA_Trans transB, FLA_Obj *packed_B )
{
}

void FLA_Gemm_pack_andor_scale_A( FLA_Trans transA, FLA_Obj alpha,
				  FLA_Obj A, FLA_Obj *packed_A )
{
  *packed_A = A;
}

void FLA_Gemm_release_pack_A( FLA_Trans transA, FLA_Obj *packed_A )
{
}

void FLA_Gemm_kernel( FLA_Obj alpha, FLA_Obj packed_A, 
		      FLA_Obj packed_B, FLA_Obj packed_C )
{
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, packed_A, packed_B,
	    FLA_ONE, packed_C );
}
