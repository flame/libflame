/*
    Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/* > \brief \b validate_gelsd.c                                              */
/* =========== DOCUMENTATION ===========                                     */
/* Definition:                                                               */
/* ===========                                                               */
/* SUBROUTINE validate_gelsd(integer m, integer n, integer NRHS, void *A,....*/
/*                           integer lda, void *B, integer ldb, void *S,     */
/*                           void *X, void *rcond, integer *rank,            */
/*                           integer datatype, double* residual);            */
/* > \par Purpose:                                                           */
/* =============                                                             */
/* >                                                                         */
/* > \verbatim                                                               */
/* >                                                                         */
/* > Defines validate function of GELSD() to use in test suite               */
/* >                                                                         */
/* > This API validates GELSD() functionality by performing                  */
/* > the folllowing tests:                                                   */
/* >                                                                         */
/* > TEST 1: If m >= n and RANK = n, the residual sum-of-squares             */
/* >         for the solution in the i-th column is given by the sum of      */
/* >         squares of elements n+1:m in that column                        */
/* >                                                                         */
/* > TEST 2: Check AX = B.....................................               */
/* >                                                                         */
/* > TEST 3: checks whether X is in the row space of A or A'                 */
/* >                                                                         */
/* > Test SVD: Check if Singular values generated are in order               */
/* >                                                                         */
/* > \endverbatim                                                            */
/* ========================================================================= */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

/* compute_matrix_norm() stores its result through an untyped pointer whose width
   follows the datatype: four bytes for the single precision types and eight for
   the double precision ones. Give it a destination of the matching type, so the
   result can be widened once and the caller can work in double throughout. */
static double matrix_norm_as_double(integer datatype, char ntype, integer m, integer n, void *A,
                                    integer lda, char imatrix, void *work)
{
    if(datatype == FLOAT || datatype == COMPLEX)
    {
        float norm = 0.0f;
        compute_matrix_norm(datatype, ntype, m, n, A, lda, &norm, imatrix, work);
        return norm;
    }
    else
    {
        double norm = 0.0;
        compute_matrix_norm(datatype, ntype, m, n, A, lda, &norm, imatrix, work);
        return norm;
    }
}

void validate_gelsd(char *tst_api, integer m, integer n, integer nrhs, void *A, integer lda,
                    void *B, integer ldb, void *S, void *X, void *rcond, integer *rank,
                    integer datatype, double err_thresh, char imatrix, void *params)
{
    void *work = NULL, *B_res = NULL;
    char NORM = '1';
    integer ldx, norm_rows = *rank < n ? n : m;
    ldx = ldb;
    double residual, resid1 = 0., resid2 = 0.;
    double resid3 = 0., resid4 = 0., resid5 = 0.;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, nrhs, &B_res, n);

    /* Early return conditions */
    if(m == 0 || n == 0 || nrhs == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, n, err_thresh);
    }

    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, n, err_thresh);

    /* Test SVD
       Check the order of Singular values generated */
    resid1 = svd_check_order(datatype, S, m, n, err_thresh);

    double norm_a = 0., norm_b = 0., norm_x = 0., norm = 0.;
    if((m >= n) && (*rank == n))
    {
        residual_sum_of_squares(datatype, m, n, nrhs, X, ldb, &resid2);
    }
    else
    {
        /* Test 3 */
        /* Compute |B-AX| = 0 */
        norm_a = matrix_norm_as_double(datatype, NORM, m, n, A, lda, imatrix, work);
        norm_b = matrix_norm_as_double(datatype, NORM, m, nrhs, B, ldb, imatrix, work);
        norm_x = matrix_norm_as_double(datatype, NORM, n, nrhs, X, ldx, imatrix, work);

        /* Compute B-AX */
        fla_invoke_gemm(datatype, "N", "N", &m, &nrhs, &n, s_n_one, A, &lda, X, &ldx, s_one, B,
                        &ldb);
        /* If rank < n, use ||A^T*(Ax-B)|| to validate least square solution X */
        if(*rank < n)
        {
            fla_invoke_gemm(datatype, "T", "N", &n, &nrhs, &m, s_n_one, A, &lda, B, &ldb, s_zero,
                            B_res, &n);
            norm = matrix_norm_as_double(datatype, NORM, norm_rows, nrhs, B_res, n, imatrix, work);
            resid3 = fla_compute_norm_based_residual(datatype, norm, norm_a, params);
        }
        else
        {
            norm = matrix_norm_as_double(datatype, NORM, norm_rows, nrhs, B, ldb, imatrix, work);
            resid3 = fla_compute_residual(datatype, 'E', norm, norm_a,
                                          (norm_x * fla_max(m, fla_max(n, nrhs)) * norm_b), params);
        }
    }
    
    /* Test 4
     * checks whether X is in the row space of A or A'
     * by scaling both X and A such that their norms are in the range
     * [sqrt(eps), 1/sqrt(eps)], then computing an LQ factorization
     * of [A',X]', and returning the norm of the trailing triangle,
     * scaled by MAX(M,N,NRHS)*eps.
     * Currently disabled because of random failures: TODO.
     */
    /*
     * if(!((imatrix == 'U') || (imatrix == 'O')))
     * {
     *     check_vector_in_rowspace(datatype, "N", m, n, nrhs, A, lda, X, ldb, resid4);
     * }
     */

    free_matrix(B_res);

    /* Test 5: Check padding rows not modified */
    resid5 = check_padding(datatype, fla_max(m, n), nrhs, X, ldb);

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(resid3, residual);
    residual = fla_test_max(resid4, residual);
    residual = fla_test_max(resid5, residual);

    FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
    FLA_PRINT_SUBTEST_STATUS(resid5, err_thresh, "05");
}
