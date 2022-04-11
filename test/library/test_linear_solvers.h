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

#endif // TEST_LINEAR_SOLVERS_H