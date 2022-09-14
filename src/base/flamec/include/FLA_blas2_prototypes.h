/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifdef FLA_ENABLE_HIP
#include <rocblas/rocblas.h>
#endif

// --- top-level wrapper prototypes --------------------------------------------

FLA_Error FLA_Gemv( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Gemvc( FLA_Trans transa, FLA_Conj conjx, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Ger( FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Gerc( FLA_Conj conjx, FLA_Conj conjy, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Hemv( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Hemvc( FLA_Uplo uplo, FLA_Conj conja, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Her( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Herc( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Her2( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Her2c( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Symv( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Syr( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Syr2( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Trmv( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trmvsx( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Trsv( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trsvsx( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );


// --- task wrapper prototypes -------------------------------------------------

FLA_Error FLA_Gemv_task( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );
FLA_Error FLA_Trsv_task( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );

FLA_Error FLA_Gemv_h_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );
FLA_Error FLA_Gemv_n_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );
FLA_Error FLA_Gemv_t_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );

FLA_Error FLA_Trsv_lc_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_ln_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_lt_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_uc_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_un_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_ut_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );


// --- external wrapper prototypes ---------------------------------------------

FLA_Error FLA_Gemv_external( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Gemvc_external( FLA_Trans transa, FLA_Conj conjx, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Ger_external( FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Gerc_external( FLA_Conj conjx, FLA_Conj conjy, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Hemv_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Hemvc_external( FLA_Uplo uplo, FLA_Conj conja, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Her_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Herc_external( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Her2_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Her2c_external( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Symv_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Syr_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Syr2_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Trmv_external( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trmvsx_external( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Trsv_external( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trsvsx_external( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );


// --- gpu wrapper prototypes --------------------------------------------------

FLA_Error FLA_Gemv_external_gpu( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu, FLA_Obj beta, FLA_Obj y, void* y_gpu );
FLA_Error FLA_Trsv_external_gpu( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu );


// --- hip wrapper prototypes --------------------------------------------------
#ifdef FLA_ENABLE_HIP
FLA_Error FLA_Gemv_external_hip( rocblas_handle handle, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu, FLA_Obj beta, FLA_Obj y, void* y_gpu );
FLA_Error FLA_Trsv_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu );
#endif

// --- check routine prototypes ------------------------------------------------

// front-ends
FLA_Error FLA_Gemv_check( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Gemvc_check( FLA_Trans transa, FLA_Conj conjx, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Ger_check( FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Gerc_check( FLA_Conj conjx, FLA_Conj conjy, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Hemv_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Hemvc_check( FLA_Uplo uplo, FLA_Conj conja, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Her_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Herc_check( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Her2_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Her2c_check( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Symv_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Syr_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Syr2_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Trmv_check( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trmvsx_check( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Trsv_check( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trsvsx_check( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );

// internal back-ends
FLA_Error FLA_Gemv_internal_check( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );
FLA_Error FLA_Trsv_internal_check( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );

