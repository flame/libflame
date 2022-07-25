/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_common.h
 *  @brief Defines function declarations of validate functions of APIs to be
 *         used in test suite.
 *  */

#ifndef VALIDATE_COMMON_H
#define VALIDATE_COMMON_H

void validate_geqrf(integer m_A,
    integer n_A,
    void *A,
    void *A_test,
    void *T_test,
    integer datatype,
    double* residual);

void validate_gerq2(integer m_A,
    integer n_A,
    void *A,
    void *A_test,
    void *T_test,
    integer datatype,
    double* residual);

void validate_gerqf(integer m_A,
    integer n_A,
    void *A,
    void *A_test,
    void *T_test,
    integer datatype,
    double* residual);

void validate_gesdd(char *jobz,
    integer m,
    integer n,
    void* A,
    void* A_test,
    void* s,
    void* U,
    void* V,
    integer datatype,
    double* residual);

void validate_gesvd(char *jobu, char *jobvt,
    integer m,
    integer n,
    void* A,
    void* A_test,
    void* s,
    void* U,
    void* V,
    integer datatype,
    double* residual);

void validate_getrf(integer m_A,
    integer n_A,
    void* A,
    void* A_test,
    integer* IPIV,
    integer datatype,
    double* residual);

void validate_getri(integer m_A,
    integer n_A,
    void* A,
    void* A_inv,
    integer* IPIV,
    integer datatype,
    double* residual);

void validate_getrs(char *trans,
    integer m,
    integer n,
    void* A,
    void* B,
    void* X,
    integer datatype,
    double* residual);

void validate_orgqr(integer m,
    integer n,
    void *A,
    void* Q,
    void *R,
    void* work,
    integer datatype,
    double* residual);

void validate_potrf(char *uplo,
    integer m,
    void *A,
    void *A_test,
    integer datatype,
    double* residual);

void validate_potrs(char *uplo,
    integer m,
    void *A,
    void *A_test,
    integer datatype,
    void *x,
    void *b,
    double* residual);

void validate_syevd(char* jobz,
    integer n,
    void* A,
    void* A_test,
    void* w,
    integer datatype,
    double* residual);

void validate_geevx(char* jobvl,
    char* jobvr,
    char* sense,
    char* balanc,
    integer m,
    void* A,
    void* A_test,
    void* VL,
    void* VR,
    void* w,
    void* wr,
    void* wi,
    void* scale,
    void* abnrm,
    void* rconde,
    void* rcondv,
    integer datatype,
    double* residual);

void validate_geev(char* jobvl,
    char* jobvr,
    integer m,
    void* A,
    void* A_test,
    void* VL,
    void* VR,
    void* w,
    void* wr,
    void* wi,
    integer datatype,
    double* residual);

void validate_geqp3(integer m_A,
    integer n_A,
    void *A,
    void *A_test,
    integer *jpvt,
    void *T_test,
    integer datatype,
    double* residual);

void validate_gesv(integer n,
    integer nrhs,
    void* A,
    void* B,
    void* X,
    integer datatype,
    double* residual);

void validate_ggev(char* jobvl,
    char* jobvr, 
    integer n, 
    void* A,
    void* B, 
    void* alpha, 
    void* alphar, 
    void* alphai, 
    void* beta,
    void* VL, 
    integer ldvl, 
    void* VR, 
    integer ldvr,
    integer datatype, 
    double *residual);

void validate_ggevx(char* balanc,
    char* jobvl,
    char* jobvr,
    char* sense,    
    integer n, 
    void* A,
    void* B, 
    void* alpha, 
    void* alphar, 
    void* alphai, 
    void* beta,
    void* VL, 
    integer ldvl, 
    void* VR, 
    integer ldvr,
    integer datatype, 
    double *residual);

#endif // VALIDATE_COMMON_H
