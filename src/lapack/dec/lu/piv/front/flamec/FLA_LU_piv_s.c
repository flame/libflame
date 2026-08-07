/*
    Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
#include "FLA_f2c.h"
#include "fla_lapack_x86_common.h"

#ifdef FLA_ENABLE_AMD_OPT
#ifdef FLA_OPENMP_MULTITHREADING
static real d__1 = -1;
static real c_b1 = 1;
#endif
extern int fla_thread_get_num_threads(void);

void FLA_get_optimum_params_sgetrf(fla_dim_t m, fla_dim_t n, fla_dim_t *nb, int *n_threads)
{
#ifdef FLA_OPENMP_MULTITHREADING
    int available_n_threads;

    /* Get maximum thread available*/
    available_n_threads = fla_thread_get_num_threads();

    if(m <= 512 || n <= 512)
    {
        *nb = 15;
        *n_threads = 8;
    }
    else if(m <= 1024 || n <= 1024)
    {
        *nb = 32;
        *n_threads = 8;
    }
    else if(m <= 2048 || n <= 2048)
    {
        *nb = 64;
        *n_threads = 24;
    }
    else if(m <= 6500 || n <= 6500)
    {
        *nb = 128;
        *n_threads = 24;
    }
    else if(m <= 12000 || n <= 12000)
    {
        *nb = 128;
        *n_threads = 32;
    }
    else
    {
        *nb = 128;
        *n_threads = 64;
    }

    if(*n_threads > available_n_threads)
        *n_threads = available_n_threads;

#else
    *nb = 64;
    *n_threads = 1;
#endif

    return;
}

/* LU factorization blocked varaiant */
int FLA_LU_piv_s_parallel(fla_dim_t *m, fla_dim_t *n, real *a, fla_dim_t *lda, aocl_int_t *ipiv,
                          fla_dim_t *info)
{
#ifdef FLA_OPENMP_MULTITHREADING
    /* System generated locals */
    fla_dim_t a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    fla_dim_t i__, j, s, iinfo;
    fla_dim_t jb, jb_prev, jb_offset, nb;
    fla_dim_t c__1 = 1;
    int threads_id, n_threads;

    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }

    if(*info != 0)
    {
        return *info;
    }

    /* Quick return if possible */
    if(*m == 0 || *n == 0)
        return 0;

    /* Determine optimum block and thread size for this environment */
    FLA_get_optimum_params_sgetrf(*m, *n, &nb, &n_threads);

    /* call sequencial algorithm for single thread */
    if(n_threads == 1)
    {
        aocl_lapack_sgetrf2(m, n, a, lda, ipiv, info);
        return *info;
    }

    /* Adjust dimension of the matrix */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --ipiv;

    /*----------------blocked LU algorithm-------------------------

    A00 |   A01               L00 |   0           U00 |   U01
    ----|-----------   ==>    ----|-------        ----|----------
        |                         |           *       |
    A10 |   A11               L10 |   L11          0  |   U11
        |                         |                   |

    1. Step 1 => compute L00 and U00
        A00 = L00 * U00

    2. Step 2 => Compute U01
        A01 = L00 * U01

    3. Step 3 => compute L10
        A10 = L10 * U00

    4. Compute L11 * U11
        A11 = L10 * U01 + L11 * U11
    ------------------------------------------------------------------*/

    j = 1;
    i__1 = fla_min(*m, *n);
    i__2 = nb;

    i__3 = fla_min(*m, *n) - j + 1;
    jb = fla_min(i__3, nb);

    /* Compute L00 and U00 of diagonal blocks */
    i__3 = *m - j + 1;

#if FLA_ENABLE_AOCL_BLAS
    /* OpenMP outer loop owns parallelism; pin inner BLIS to 1 thread to
    avoid oversubscription when BLIS would otherwise spawn its own 
    threads alongside OpenMP. */
    aocl_int64_t orig_blis_threads = bli_thread_get_num_threads();
    bli_thread_set_num_threads(1);
#endif

    aocl_lapack_sgetrf2(&i__3, &jb, M_PTR(a, j, j, a_dim1), lda, &ipiv[j], &iinfo);

    if(*info == 0 && iinfo > 0)
        *info = iinfo + j - 1;

    /* Computing MIN */
    i__4 = *m, i__5 = j + jb - 1;
    i__3 = fla_min(i__4, i__5);
    for(i__ = j; i__ <= i__3; ++i__)
    {
        ipiv[i__] = (aocl_int_t)(j - 1 + ipiv[i__]);
    }

#pragma omp parallel num_threads(n_threads) private(i__3, i__4, i__5, i__6, i__7, i__8, j, s, jb, \
                                                    jb_prev, jb_offset, threads_id)
    {
        threads_id = omp_get_thread_num();
        for(j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
        {
            s = j + i__2;
            /* Factorize the next block while the rest of the matrix is being updated */
            if(threads_id == 0)
            {
                /* Computing MIN */
                i__3 = fla_min(*m, *n) - j + 1;
                jb_prev = fla_min(i__3, nb);
                /* Apply interchanges to columns J+JB:N */
                i__3 = *n - j - jb_prev + 1;
                i__4 = j + jb_prev - 1;
                i__3 = fla_min(i__3, jb_prev);
                aocl_lapack_slaswp(&i__3, M_PTR(a, 1, j + jb_prev, a_dim1), lda, &j, &i__4, &ipiv[1], &c__1);
                /* compute U10 */
                i__3 = *n - j - jb_prev + 1;
                i__3 = fla_min(i__3, jb_prev);
                aocl_blas_strsm("Left", "Lower", "No transpose", "Unit", &jb_prev, &i__3, &c_b1,
                       M_PTR(a, j, j, a_dim1), lda, M_PTR(a, j, j + jb_prev, a_dim1), lda);
                /* compute L11 * U11 */
                if(j + jb_prev <= *m)
                {
                    /* Update trailing submatrix. */
                    i__3 = *m - j - jb_prev + 1;
                    i__4 = *n - j - jb_prev + 1;
                    i__4 = fla_min(i__4, jb_prev);
                    aocl_blas_sgemm("No transpose", "No transpose", &i__3, &i__4, &jb_prev, &d__1,
                           M_PTR(a, j + jb_prev, j, a_dim1), lda, M_PTR(a, j, j + jb_prev, a_dim1), lda,
                           &c_b1, M_PTR(a, j + jb_prev, j + jb_prev, a_dim1), lda);
                }

                if(s <= i__1)
                {
                    /* Computing MIN */
                    i__3 = fla_min(*m, *n) - s + 1;
                    jb = fla_min(i__3, nb);

                    /* Compute L00 and U00 of diagonal blocks */
                    i__3 = *m - s + 1;
                    aocl_lapack_sgetrf2(&i__3, &jb, M_PTR(a, s, s, a_dim1), lda, &ipiv[s], &iinfo);

                    if(*info == 0 && iinfo > 0)
                        *info = iinfo + s - 1;
                    /* Computing MIN */
                    i__4 = *m, i__5 = s + jb - 1;
                    i__3 = fla_min(i__4, i__5);
                    for(i__ = s; i__ <= i__3; ++i__)
                    {
                        ipiv[i__] = (aocl_int_t)(s - 1 + ipiv[i__]);
                    }
                }
            }
            else
            {
                /* Computing MIN */
                i__3 = fla_min(*m, *n) - j + 1;
                jb_prev = fla_min(i__3, nb);
                jb_offset = jb_prev * 2;

                if(j + jb_prev <= *n)
                {
                    /* Apply interchanges to columns J+JB:N */
                    i__3 = *n - j - jb_prev + 1 - jb_prev;
                    i__4 = j + jb_prev - 1;
                    FLA_Thread_get_subrange(threads_id - 1, n_threads - 1, i__3, &i__5, &i__6);
                    aocl_lapack_slaswp(&i__5, M_PTR(a, 1, j + jb_offset + i__6, a_dim1), lda, &j, &i__4,
                            &ipiv[1], &c__1);

                    /* compute U10 */
                    i__3 = *n - j - jb_prev + 1 - jb_prev;
                    FLA_Thread_get_subrange(threads_id - 1, n_threads - 1, i__3, &i__5, &i__6);
                    aocl_blas_strsm("Left", "Lower", "No transpose", "Unit", &jb_prev, &i__5, &c_b1,
                           M_PTR(a, j, j, a_dim1), lda, M_PTR(a, j, j + jb_offset + i__6, a_dim1), lda);

                    /* compute L11 * U11 */
                    if(j + jb_prev <= *m)
                    {
                        /* Update trailing submatrix. */
                        i__3 = *m - j - jb_prev + 1;
                        i__4 = *n - j - jb_prev + 1 - jb_prev;
                        FLA_Thread_get_subrange(threads_id - 1, n_threads - 1, i__4, &i__7, &i__8);
                        aocl_blas_sgemm("No transpose", "No transpose", &i__3, &i__7, &jb_prev, &d__1,
                               M_PTR(a, j + jb_prev, j, a_dim1), lda,
                               M_PTR(a, j, j + jb_offset + i__8, a_dim1), lda, &c_b1,
                               M_PTR(a, j + jb_prev, j + jb_offset + i__8, a_dim1), lda);
                    }
                }
            }
#pragma omp barrier
        }
    }

#pragma omp parallel num_threads(n_threads) private(j, i__3, i__4, i__5, i__6, jb, threads_id)
    {
        threads_id = omp_get_thread_num();
        for(j = 1; j <= i__1; j += i__2)
        {
            /* Computing MIN */
            i__3 = fla_min(*m, *n) - j + 1;
            jb = fla_min(i__3, nb);

            i__3 = j - 1;
            i__4 = j + jb - 1;
            FLA_Thread_get_subrange(threads_id, n_threads, i__3, &i__5, &i__6);
            aocl_lapack_slaswp(&i__5, (real *)&a[a_offset + (i__6 * a_dim1)], lda, &j, &i__4, &ipiv[1], &c__1);
#pragma omp barrier
        }
    }

#if FLA_ENABLE_AOCL_BLAS
    bli_thread_set_num_threads(orig_blis_threads);
#endif

#else
    aocl_lapack_sgetrf2(m, n, a, lda, ipiv, info);
#endif

    return *info;
}

#endif
