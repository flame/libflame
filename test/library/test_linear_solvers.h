/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file test_linear_solvers.h
 *  @brief Defines function declarations to use in linear solver APIs
           of test suite.
 *  */

#ifndef TEST_LINEAR_SOLVERS_H
#define TEST_LINEAR_SOLVERS_H

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

void validate_potrs(char *uplo,
    integer m,
    void *A,
    void *A_test,
    integer datatype,
    void *x,
    void *b,
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

void validate_getrs(integer m,
	integer n,
	void* A,
	void* B,
	void* X,
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
#endif // TEST_LINEAR_SOLVERS_H
