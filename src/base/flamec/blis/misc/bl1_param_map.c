
#include "blis1.h"

// --- BLIS to BLAS/LAPACK mappings --------------------------------------------

void bl1_param_map_to_netlib_trans( trans1_t blis_trans, void* blas_trans )
{
	if ( bl1_is_notrans( blis_trans ) || bl1_is_conjnotrans( blis_trans ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasNoTrans;
#else
		*( ( char*                 ) blas_trans ) = 'N';
#endif
	}
	else if ( bl1_is_trans( blis_trans ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasTrans;
#else
		*( ( char*                 ) blas_trans ) = 'T';
#endif
	}
	else if ( bl1_is_conjtrans( blis_trans ))
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasConjTrans;
#else
		*( ( char*                 ) blas_trans ) = 'C';
#endif
	}
	else
	{
		bl1_abort_msg( "Invalid BLIS trans value to map." );
	}
}

void bl1_param_map_to_netlib_uplo( uplo1_t blis_uplo, void* blas_uplo )
{
	if ( bl1_is_lower( blis_uplo ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_UPLO* ) blas_uplo ) = CblasLower;
#else
		*( ( char*            ) blas_uplo ) = 'L';
#endif
	}
	else if ( bl1_is_upper( blis_uplo ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_UPLO* ) blas_uplo ) = CblasUpper;
#else
		*( ( char*            ) blas_uplo ) = 'U';
#endif
	}
	else
	{
		bl1_abort_msg( "Invalid BLIS uplo value to map." );
	}
}

void bl1_param_map_to_netlib_side( side1_t blis_side, void* blas_side )
{
	if ( bl1_is_left( blis_side ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_SIDE* ) blas_side ) = CblasLeft;
#else
		*( ( char*            ) blas_side ) = 'L';
#endif
	}
	else if ( bl1_is_right( blis_side ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_SIDE* ) blas_side ) = CblasRight;
#else
		*( ( char*            ) blas_side ) = 'R';
#endif
	}
	else
	{
		bl1_abort_msg( "Invalid BLIS side value to map." );
	}
}

void bl1_param_map_to_netlib_diag( diag1_t blis_diag, void* blas_diag )
{
	if ( bl1_is_nonunit_diag( blis_diag ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_DIAG* ) blas_diag ) = CblasNonUnit;
#else
		*( ( char*            ) blas_diag ) = 'N';
#endif
	}
	else if ( bl1_is_unit_diag( blis_diag ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_DIAG* ) blas_diag ) = CblasUnit;
#else
		*( ( char*            ) blas_diag ) = 'U';
#endif
	}
	else
	{
		bl1_abort_msg( "Invalid BLIS diag value to map." );
	}
}

