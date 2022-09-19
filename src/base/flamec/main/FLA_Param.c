/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#ifdef FLA_ENABLE_HIP
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>
#endif

// --- FLAME to BLAS/LAPACK mappings -------------------------------------------

void FLA_Param_map_flame_to_netlib_trans( FLA_Trans trans, void* blas_trans )
{
	if ( trans == FLA_NO_TRANSPOSE )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasNoTrans;
#else
		*( ( char*                 ) blas_trans ) = 'N';
#endif
	}
	else if ( trans == FLA_TRANSPOSE )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasTrans;
#else
		*( ( char*                 ) blas_trans ) = 'T';
#endif
	}
	else if ( trans == FLA_CONJ_TRANSPOSE )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasConjTrans;
#else
		*( ( char*                 ) blas_trans ) = 'C';
#endif
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_TRANS );
	}
}

#ifdef FLA_ENABLE_HIP
rocblas_operation FLA_Param_map_flame_to_rocblas_trans( FLA_Trans trans, FLA_Bool is_real )
{
        if ( trans == FLA_NO_TRANSPOSE )
        {
                return rocblas_operation_none;
        } else if ( trans == FLA_TRANSPOSE )
        {
                return rocblas_operation_transpose;
        }
        else if ( trans == FLA_CONJ_TRANSPOSE && is_real)
        {
                return rocblas_operation_transpose;
        }
        else if ( trans == FLA_CONJ_NO_TRANSPOSE && is_real )
        {
                return rocblas_operation_none;
        }
        else if ( trans == FLA_CONJ_TRANSPOSE && !is_real)
        {
                return rocblas_operation_conjugate_transpose;
        }
        else if ( trans == FLA_CONJ_NO_TRANSPOSE && !is_real)
        {
                // not supported by rocBLAS
                fprintf( stderr, "FLA_CONJ_NO_TRANSPOSE not supported by rocBLAS.\n" );
                FLA_Check_error_code( FLA_INVALID_TRANS );
                return rocblas_operation_none; // to silence warning
        }
        else
        {
                FLA_Check_error_code( FLA_INVALID_TRANS );
                return rocblas_operation_none; // to silence warning
        }
}
#endif

void FLA_Param_map_flame_to_netlib_uplo( FLA_Uplo uplo, void* blas_uplo )
{
	if ( uplo == FLA_LOWER_TRIANGULAR )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_UPLO* ) blas_uplo ) = CblasLower;
#else
		*( ( char*            ) blas_uplo ) = 'L';
#endif
	}
	else if ( uplo == FLA_UPPER_TRIANGULAR )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_UPLO* ) blas_uplo ) = CblasUpper;
#else
		*( ( char*            ) blas_uplo ) = 'U';
#endif
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_UPLO );
	}
}

#ifdef FLA_ENABLE_HIP
rocblas_fill FLA_Param_map_flame_to_rocblas_uplo( FLA_Uplo uplo )
{
	if ( uplo == FLA_LOWER_TRIANGULAR )
	{
		return rocblas_fill_lower;
	}
	else if ( uplo == FLA_UPPER_TRIANGULAR )
	{
		return rocblas_fill_upper;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_UPLO );
		return rocblas_fill_lower; // to silence warning
	}
}
#endif

void FLA_Param_map_flame_to_netlib_side( FLA_Side side, void* blas_side )
{
	if ( side == FLA_LEFT )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_SIDE* ) blas_side ) = CblasLeft;
#else
		*( ( char*            ) blas_side ) = 'L';
#endif
	}
	else if ( side == FLA_RIGHT )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_SIDE* ) blas_side ) = CblasRight;
#else
		*( ( char*            ) blas_side ) = 'R';
#endif
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_SIDE );
	}
}

#ifdef FLA_ENABLE_HIP
rocblas_side FLA_Param_map_flame_to_rocblas_side( FLA_Side side )
{
	if ( side == FLA_LEFT )
	{
		return rocblas_side_left;
	}
	else if ( side == FLA_RIGHT )
	{
		return rocblas_side_right;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_SIDE );
		return rocblas_side_left; // to silence warning
	}
}
#endif

void FLA_Param_map_flame_to_netlib_diag( FLA_Diag diag, void* blas_diag )
{
	if ( diag == FLA_NONUNIT_DIAG )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_DIAG* ) blas_diag ) = CblasNonUnit;
#else
		*( ( char*            ) blas_diag ) = 'N';
#endif
	}
	else if ( diag == FLA_UNIT_DIAG )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_DIAG* ) blas_diag ) = CblasUnit;
#else
		*( ( char*            ) blas_diag ) = 'U';
#endif
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_DIAG );
	}
}

#ifdef FLA_ENABLE_HIP
rocblas_diagonal FLA_Param_map_flame_to_rocblas_diag( FLA_Diag diag )
{
	if ( diag == FLA_NONUNIT_DIAG )
	{
		return rocblas_diagonal_non_unit;
	}
	else if ( diag == FLA_UNIT_DIAG )
	{
		return rocblas_diagonal_unit;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_DIAG );
		return rocblas_diagonal_non_unit; // to silence warning
	}
}
#endif

void FLA_Param_map_flame_to_netlib_direct( FLA_Direct direct, void* lapack_direct )
{
	if ( direct == FLA_FORWARD )
	{
		*( ( char* ) lapack_direct ) = 'F';
	}
	else if ( direct == FLA_BACKWARD )
	{
		*( ( char* ) lapack_direct ) = 'B';
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_DIRECT );
	}
}

void FLA_Param_map_flame_to_netlib_storev( FLA_Store storev, void* lapack_storev )
{
	if ( storev == FLA_COLUMNWISE )
	{
		*( ( char* ) lapack_storev ) = 'C';
	}
	else if ( storev == FLA_ROWWISE )
	{
		*( ( char* ) lapack_storev ) = 'R';
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_STOREV );
	}
}

void FLA_Param_map_flame_to_netlib_evd_type( FLA_Evd_type evd_type, void* lapack_evd_type )
{
	if ( evd_type == FLA_EVD_WITHOUT_VECTORS )
	{
		*( ( char* ) lapack_evd_type ) = 'N';
	}
	else if ( evd_type == FLA_EVD_WITH_VECTORS )
	{
		*( ( char* ) lapack_evd_type ) = 'V';
	}
	else if ( evd_type == FLA_EVD_OF_TRIDIAG_WITH_VECTORS )
	{
		*( ( char* ) lapack_evd_type ) = 'I';
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_EVD_TYPE );
	}
}

#ifdef FLA_ENABLE_HIP
rocblas_evect FLA_Param_map_flame_to_rocblas_evd_type( FLA_Evd_type evd_type )
{
        if ( evd_type == FLA_EVD_WITHOUT_VECTORS )
        {
                return rocblas_evect_none;
        }
        else if ( evd_type == FLA_EVD_WITH_VECTORS )
        {
                return rocblas_evect_original;
        }
        else if ( evd_type == FLA_EVD_OF_TRIDIAG_WITH_VECTORS )
        {
                return rocblas_evect_tridiagonal;
        }
        else
        {
                FLA_Check_error_code( FLA_INVALID_EVD_TYPE );
        }
}
#endif

void FLA_Param_map_flame_to_netlib_svd_type( FLA_Svd_type svd_type, void* lapack_svd_type )
{
	if      ( svd_type == FLA_SVD_VECTORS_ALL )
	{
		*( ( char* ) lapack_svd_type ) = 'A';
	}
	else if ( svd_type == FLA_SVD_VECTORS_MIN_COPY )
	{
		*( ( char* ) lapack_svd_type ) = 'S';
	}
	else if ( svd_type == FLA_SVD_VECTORS_MIN_OVERWRITE )
	{
		*( ( char* ) lapack_svd_type ) = 'O';
	}
	else if ( svd_type == FLA_SVD_VECTORS_NONE )
	{
		*( ( char* ) lapack_svd_type ) = 'N';
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_SVD_TYPE );
	}
}

#ifdef FLA_ENABLE_HIP
rocblas_svect     FLA_Param_map_flame_to_rocblas_svd_type( FLA_Svd_type svd_type )
{
        if      ( svd_type == FLA_SVD_VECTORS_ALL )
        {
                return rocblas_svect_all;
        }
        else if ( svd_type == FLA_SVD_VECTORS_MIN_COPY )
        {
                return rocblas_svect_singular;
        }
        else if ( svd_type == FLA_SVD_VECTORS_MIN_OVERWRITE )
        {
                return rocblas_svect_overwrite;
        }
        else if ( svd_type == FLA_SVD_VECTORS_NONE )
        {
                return rocblas_svect_none;
        }
        else
        {
                FLA_Check_error_code( FLA_INVALID_SVD_TYPE );
        }
}

#endif

void FLA_Param_map_flame_to_netlib_machval( FLA_Machval machval, void* blas_machval )
{
	if      ( machval == FLA_MACH_EPS )
	{
		*( ( char* ) blas_machval ) = 'E';
	}
	else if ( machval == FLA_MACH_SFMIN )
	{
		*( ( char* ) blas_machval ) = 'S';
	}
	else if ( machval == FLA_MACH_BASE )
	{
		*( ( char* ) blas_machval ) = 'B';
	}
	else if ( machval == FLA_MACH_PREC )
	{
		*( ( char* ) blas_machval ) = 'P';
	}
	else if ( machval == FLA_MACH_NDIGMANT )
	{
		*( ( char* ) blas_machval ) = 'N';
	}
	else if ( machval == FLA_MACH_RND )
	{
		*( ( char* ) blas_machval ) = 'R';
	}
	else if ( machval == FLA_MACH_EMIN )
	{
		*( ( char* ) blas_machval ) = 'M';
	}
	else if ( machval == FLA_MACH_RMIN )
	{
		*( ( char* ) blas_machval ) = 'U';
	}
	else if ( machval == FLA_MACH_EMAX )
	{
		*( ( char* ) blas_machval ) = 'L';
	}
	else if ( machval == FLA_MACH_RMAX )
	{
		*( ( char* ) blas_machval ) = 'O';
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_MACHVAL );
	}
}

// --- FLAME to BLIS mappings --------------------------------------------------

void FLA_Param_map_flame_to_blis_trans( FLA_Trans trans, trans1_t* blis_trans )
{
	if ( trans == FLA_NO_TRANSPOSE )
	{
		*blis_trans = BLIS1_NO_TRANSPOSE;
	}
	else if ( trans == FLA_TRANSPOSE )
	{
		*blis_trans = BLIS1_TRANSPOSE;
	}
	else if ( trans == FLA_CONJ_NO_TRANSPOSE )
	{
		*blis_trans = BLIS1_CONJ_NO_TRANSPOSE;
	}
	else if ( trans == FLA_CONJ_TRANSPOSE )
	{
		*blis_trans = BLIS1_CONJ_TRANSPOSE;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_TRANS );
	}
}

void FLA_Param_map_flame_to_blis_conj( FLA_Conj conj, conj1_t* blis_conj )
{
	if ( conj == FLA_NO_CONJUGATE )
	{
		*blis_conj = BLIS1_NO_CONJUGATE;
	}
	else if ( conj == FLA_CONJUGATE )
	{
		*blis_conj = BLIS1_CONJUGATE;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_CONJ );
	}
}

void FLA_Param_map_flame_to_blis_uplo( FLA_Uplo uplo, uplo1_t* blis_uplo )
{
	if ( uplo == FLA_LOWER_TRIANGULAR )
	{
		*blis_uplo = BLIS1_LOWER_TRIANGULAR;
	}
	else if ( uplo == FLA_UPPER_TRIANGULAR )
	{
		*blis_uplo = BLIS1_UPPER_TRIANGULAR;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_UPLO );
	}
}

void FLA_Param_map_flame_to_blis_side( FLA_Side side, side1_t* blis_side )
{
	if ( side == FLA_LEFT )
	{
		*blis_side = BLIS1_LEFT;
	}
	else if ( side == FLA_RIGHT )
	{
		*blis_side = BLIS1_RIGHT;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_SIDE );
	}
}

void FLA_Param_map_flame_to_blis_diag( FLA_Diag diag, diag1_t* blis_diag )
{
	if ( diag == FLA_NONUNIT_DIAG )
	{
		*blis_diag = BLIS1_NONUNIT_DIAG;
	}
	else if ( diag == FLA_UNIT_DIAG )
	{
		*blis_diag = BLIS1_UNIT_DIAG;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_DIAG );
	}
}

// --- BLAS/LAPACK to FLAME mappings -------------------------------------------

void FLA_Param_map_netlib_to_flame_trans( char* trans, FLA_Trans* flame_trans )
{
	if      ( *trans == 'n' || *trans == 'N' )
		*flame_trans = FLA_NO_TRANSPOSE;
	else if ( *trans == 't' || *trans == 'T' )
		*flame_trans = FLA_TRANSPOSE;
	else if ( *trans == 'c' || *trans == 'C' )
		*flame_trans = FLA_CONJ_TRANSPOSE;
	else
		FLA_Check_error_code( FLA_INVALID_TRANS );
}

void FLA_Param_map_netlib_to_flame_uplo( char* uplo, FLA_Uplo* flame_uplo )
{
	if      ( *uplo == 'l' || *uplo == 'L' )
		*flame_uplo = FLA_LOWER_TRIANGULAR;
	else if ( *uplo == 'u' || *uplo == 'U' )
		*flame_uplo = FLA_UPPER_TRIANGULAR;
	else
		FLA_Check_error_code( FLA_INVALID_UPLO );
}

void FLA_Param_map_netlib_to_flame_side( char* side, FLA_Side* flame_side )
{
	if      ( *side == 'l' || *side == 'L' )
		*flame_side = FLA_LEFT;
	else if ( *side == 'r' || *side == 'R' )
		*flame_side = FLA_RIGHT;
	else
		FLA_Check_error_code( FLA_INVALID_SIDE );
}

void FLA_Param_map_netlib_to_flame_diag( char* diag, FLA_Diag* flame_diag )
{
	if      ( *diag == 'n' || *diag == 'N' )
		*flame_diag = FLA_NONUNIT_DIAG;
	else if ( *diag == 'u' || *diag == 'U' )
		*flame_diag = FLA_UNIT_DIAG;
	else
		FLA_Check_error_code( FLA_INVALID_DIAG );
}

void FLA_Param_map_netlib_to_flame_inv( int* itype, FLA_Inv* flame_inv )
{
	if      ( *itype == 1 )
		*flame_inv = FLA_INVERSE;
	else if ( *itype == 2 || *itype == 3 )
		*flame_inv = FLA_NO_INVERSE;
	else
		FLA_Check_error_code( FLA_INVALID_INVERSE );
}

void FLA_Param_map_netlib_to_flame_svd_type( char* svd, FLA_Svd_type* flame_svd )
{
        if      ( *svd == 'A' || *svd == 'a' ) 
                *flame_svd = FLA_SVD_VECTORS_ALL;
	else if ( *svd == 'S' || *svd == 's' ) 
                *flame_svd = FLA_SVD_VECTORS_MIN_COPY;
	else if ( *svd == 'O' || *svd == 'o' ) 
                *flame_svd = FLA_SVD_VECTORS_MIN_OVERWRITE;
	else if ( *svd == 'N' || *svd == 'n' ) 
                *flame_svd = FLA_SVD_VECTORS_NONE;
	else
		FLA_Check_error_code( FLA_INVALID_SVD_TYPE );
}


// --- BLIS to FLAME mappings --------------------------------------------------

void FLA_Param_map_blis_to_flame_trans( trans1_t trans, FLA_Trans* flame_trans )
{
	if      ( bl1_is_notrans( trans ) )
		*flame_trans = FLA_NO_TRANSPOSE;
	else if ( bl1_is_trans( trans ) )
		*flame_trans = FLA_TRANSPOSE;
	else if ( bl1_is_conjnotrans( trans ) )
		*flame_trans = FLA_CONJ_NO_TRANSPOSE;
	else if ( bl1_is_conjtrans( trans ) )
		*flame_trans = FLA_CONJ_TRANSPOSE;
	else
		FLA_Check_error_code( FLA_INVALID_TRANS );
}

void FLA_Param_map_blis_to_flame_uplo( uplo1_t uplo, FLA_Uplo* flame_uplo )
{
	if      ( bl1_is_lower( uplo ) )
		*flame_uplo = FLA_LOWER_TRIANGULAR;
	else if ( bl1_is_upper( uplo ) )
		*flame_uplo = FLA_UPPER_TRIANGULAR;
	else
		FLA_Check_error_code( FLA_INVALID_UPLO );
}

void FLA_Param_map_blis_to_flame_side( side1_t side, FLA_Side* flame_side )
{
	if      ( bl1_is_left( side ) )
		*flame_side = FLA_LEFT;
	else if ( bl1_is_right( side ) )
		*flame_side = FLA_RIGHT;
	else
		FLA_Check_error_code( FLA_INVALID_SIDE );
}

void FLA_Param_map_blis_to_flame_diag( diag1_t diag, FLA_Diag* flame_diag )
{
	if      ( bl1_is_nonunit_diag( diag ) )
		*flame_diag = FLA_NONUNIT_DIAG;
	else if ( bl1_is_unit_diag( diag ) )
		*flame_diag = FLA_UNIT_DIAG;
	else if ( bl1_is_zero_diag( diag ) )
		*flame_diag = FLA_ZERO_DIAG;
	else
		FLA_Check_error_code( FLA_INVALID_DIAG );
}

// --- FLAME char to FLAME mappings --------------------------------------------

void FLA_Param_map_char_to_flame_trans( char* trans, FLA_Trans* flame_trans )
{
	if      ( *trans == 'n' || *trans == 'N' )
		*flame_trans = FLA_NO_TRANSPOSE;
	else if ( *trans == 't' || *trans == 'T' )
		*flame_trans = FLA_TRANSPOSE;
	else if ( *trans == 'c' || *trans == 'C' )
		*flame_trans = FLA_CONJ_NO_TRANSPOSE;
	else if ( *trans == 'h' || *trans == 'H' )
		*flame_trans = FLA_CONJ_TRANSPOSE;
	else
		FLA_Check_error_code( FLA_INVALID_TRANS );
}

void FLA_Param_map_char_to_flame_uplo( char* uplo, FLA_Uplo* flame_uplo )
{
	if      ( *uplo == 'l' || *uplo == 'L' )
		*flame_uplo = FLA_LOWER_TRIANGULAR;
	else if ( *uplo == 'u' || *uplo == 'U' )
		*flame_uplo = FLA_UPPER_TRIANGULAR;
	else
		FLA_Check_error_code( FLA_INVALID_UPLO );
}

void FLA_Param_map_char_to_flame_side( char* side, FLA_Side* flame_side )
{
	if      ( *side == 'l' || *side == 'L' )
		*flame_side = FLA_LEFT;
	else if ( *side == 'r' || *side == 'R' )
		*flame_side = FLA_RIGHT;
	else
		FLA_Check_error_code( FLA_INVALID_SIDE );
}

void FLA_Param_map_char_to_flame_diag( char* diag, FLA_Diag* flame_diag )
{
	if      ( *diag == 'n' || *diag == 'N' )
		*flame_diag = FLA_NONUNIT_DIAG;
	else if ( *diag == 'u' || *diag == 'U' )
		*flame_diag = FLA_UNIT_DIAG;
	else
		FLA_Check_error_code( FLA_INVALID_DIAG );
}

void FLA_Param_map_char_to_flame_direct( char* direct, FLA_Direct* flame_direct )
{
	if      ( *direct == 'b' || *direct == 'B' )
		*flame_direct = FLA_BACKWARD;
	else if ( *direct == 'f' || *direct == 'F' )
		*flame_direct = FLA_FORWARD;
	else
		FLA_Check_error_code( FLA_INVALID_DIRECT );
}

void FLA_Param_map_char_to_flame_storev( char* storev, FLA_Direct* flame_storev )
{
	if      ( *storev == 'c' || *storev == 'C' )
		*flame_storev = FLA_COLUMNWISE;
	else if ( *storev == 'r' || *storev == 'R' )
		*flame_storev = FLA_ROWWISE;
	else
		FLA_Check_error_code( FLA_INVALID_STOREV );
}

void FLA_Param_map_char_to_flame_inv( char* inv, FLA_Inv* flame_inv )
{
	if      ( *inv == 'i' || *inv == 'I' )
		*flame_inv = FLA_INVERSE;
	else if ( *inv == 'n' || *inv == 'N' )
		*flame_inv = FLA_NO_INVERSE;
	else
		FLA_Check_error_code( FLA_INVALID_INVERSE );
}

