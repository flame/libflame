/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
    Modifications Copyright (c) 2021-2026 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLA_LU_piv_vars.h"

void FLA_get_optimum_params_getrf(fla_dim_t m, fla_dim_t n, fla_dim_t *nb, int *n_threads);
FLA_Error FLA_LU_piv_internal( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl );
fla_dim_t   FLA_LU_piv_small_s_var0( fla_dim_t *m, fla_dim_t *n, real *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv, fla_dim_t *info );
fla_dim_t   FLA_LU_piv_small_s_var1( fla_dim_t *m, fla_dim_t *n, real *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv, fla_dim_t *info );
fla_dim_t   FLA_LU_piv_small_d_var0( fla_dim_t *m, fla_dim_t *n, doublereal *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv, fla_dim_t *info );
fla_dim_t   FLA_LU_piv_small_d_var1( fla_dim_t *m, fla_dim_t *n, doublereal *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv, fla_dim_t *info );
fla_dim_t   FLA_LU_piv_small_d_var2( fla_dim_t *m, fla_dim_t *n, doublereal *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv, fla_dim_t *info );
fla_dim_t   FLA_LU_piv_c_var0(fla_dim_t* m, fla_dim_t* n, scomplex* a, fla_dim_t* lda, 
                                   aocl_int_t* ipiv, fla_dim_t* info);
fla_dim_t   FLA_LU_piv_c_parallel( fla_dim_t* m, fla_dim_t* n, scomplex* a, fla_dim_t* lda, 
                                   aocl_int_t* ipiv, fla_dim_t* info);
int   FLA_LU_piv_small_z_var0( fla_dim_t *m, fla_dim_t *n, dcomplex *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv, fla_dim_t *info);
int   FLA_LU_piv_z_var0(fla_dim_t *m, fla_dim_t *n, dcomplex *a, fla_dim_t *lda, aocl_int_t *ipiv, fla_dim_t *info);
int   FLA_LU_piv_z_parallel( fla_dim_t *m, fla_dim_t *n, dcomplex *a, fla_dim_t *lda, aocl_int_t *ipiv, fla_dim_t *info);
int   FLA_LU_piv_z_var1_parallel(fla_dim_t* m, fla_dim_t* n, dcomplex* a, fla_dim_t* lda, aocl_int_t* ipiv, fla_dim_t* info);
fla_dim_t   FLA_LU_piv_small_s_var0( fla_dim_t *m, fla_dim_t *n, real *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv, fla_dim_t *info );
fla_dim_t   FLA_LU_piv_small_s_var1( fla_dim_t *m, fla_dim_t *n, real *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv, fla_dim_t *info );
int FLA_LU_piv_d_parallel( fla_dim_t *m, fla_dim_t *n, doublereal *a, fla_dim_t *lda, aocl_int_t *ipiv, fla_dim_t *info);

int FLA_LU_piv_s_parallel( fla_dim_t *m, fla_dim_t *n, real *a, fla_dim_t *lda, aocl_int_t *ipiv, fla_dim_t *info);

FLA_Error FLA_LU_piv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj B, FLA_Obj X );

FLA_Error FLASH_LU_piv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj B, FLA_Obj X );

fla_dim_t lapack_cgetf2(fla_dim_t *m, fla_dim_t *n, scomplex *a, fla_dim_t *lda,
	 aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_cgetrf(fla_dim_t *m, fla_dim_t *n, scomplex *a, fla_dim_t *lda,
	 aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_dgetrf(fla_dim_t *m, fla_dim_t *n, doublereal *a, fla_dim_t *
	lda, aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_sgetf2(fla_dim_t *m, fla_dim_t *n, real *a, fla_dim_t *lda,
	aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_sgetrf(fla_dim_t *m, fla_dim_t *n, real *a, fla_dim_t *lda,
	aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_zgetf2(fla_dim_t *m, fla_dim_t *n, dcomplex *a,
	fla_dim_t *lda, aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_zgetrf(fla_dim_t *m, fla_dim_t *n, dcomplex *a,
	fla_dim_t *lda, aocl_int_t *ipiv, fla_dim_t *info);
