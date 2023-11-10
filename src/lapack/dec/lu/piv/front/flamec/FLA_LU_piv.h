/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
    Modifications Copyright (c) 2021-2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLA_LU_piv_vars.h"

FLA_Error FLA_LU_piv_internal( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl );
integer   FLA_LU_piv_small_s_var0( integer *m, integer *n, real *a, integer *lda,
                                   integer *ipiv, integer *info );
integer   FLA_LU_piv_small_s_var1( integer *m, integer *n, real *a, integer *lda,
                                   integer *ipiv, integer *info );
integer   FLA_LU_piv_small_d_var0( integer *m, integer *n, doublereal *a, integer *lda,
                                   integer *ipiv, integer *info );
integer   FLA_LU_piv_small_d_var1( integer *m, integer *n, doublereal *a, integer *lda,
                                   integer *ipiv, integer *info );
integer   FLA_LU_piv_small_d_var2( integer *m, integer *n, doublereal *a, integer *lda,
                                   integer *ipiv, integer *info );
int   FLA_LU_piv_small_z_var0( integer *m, integer *n, dcomplex *a, integer *lda,
                                   integer *ipiv, integer *info);
int   FLA_LU_piv_z_var0(integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info);
int   FLA_LU_piv_z_parallel( integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info);
integer   FLA_LU_piv_small_s_var0( integer *m, integer *n, real *a, integer *lda,
                                   integer *ipiv, integer *info );
integer   FLA_LU_piv_small_s_var1( integer *m, integer *n, real *a, integer *lda,
                                   integer *ipiv, integer *info );


FLA_Error FLA_LU_piv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj B, FLA_Obj X );

FLA_Error FLASH_LU_piv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj B, FLA_Obj X );

integer lapack_cgetf2(integer *m, integer *n, scomplex *a, integer *lda,
	 integer *ipiv, integer *info);
integer lapack_cgetrf(integer *m, integer *n, scomplex *a, integer *lda,
	 integer *ipiv, integer *info);
integer lapack_dgetf2(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);
integer lapack_dgetrf(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);
integer lapack_sgetf2(integer *m, integer *n, real *a, integer *lda,
	integer *ipiv, integer *info);
integer lapack_sgetrf(integer *m, integer *n, real *a, integer *lda,
	integer *ipiv, integer *info);
integer lapack_zgetf2(integer *m, integer *n, dcomplex *a,
	integer *lda, integer *ipiv, integer *info);
integer lapack_zgetrf(integer *m, integer *n, dcomplex *a,
	integer *lda, integer *ipiv, integer *info);
