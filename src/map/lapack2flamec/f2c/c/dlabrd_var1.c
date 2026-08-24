/*
    Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
*/
#include "FLA_f2c.h" /* Table of constant values */
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

#ifdef FLA_OPENMP_MULTITHREADING
static doublereal neg_one = -1.;
static doublereal d_one = 1.;
static aocl_int64_t i_one = 1;
static doublereal d_zero = 0.;
void fla_dlabrd_var1(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, doublereal *a,
                     aocl_int64_t *lda, doublereal *d__, doublereal *e, doublereal *tauq,
                     doublereal *taup, doublereal *x, aocl_int64_t *ldx, doublereal *y,
                     aocl_int64_t *ldy, volatile FLA_BARRIER *barrier,
                     doublereal *gemv_a_row_buffer, int requested_num_threads)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, x_dim1, x_offset, y_dim1, y_offset, i__1, i__2, i__3;
    /* Local variables */
    aocl_int64_t i__;
    int thread_id;
#ifdef FLA_OPENMP_MULTITHREADING
    aocl_int64_t i__4, i__5;
    int actual_num_threads = 1;
    int optimal_num_threads = 1;
#endif
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */

    /* Initialize global context data */
    aocl_fla_init();
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    if(*m <= 0 || *n <= 0)
    {
        return;
    }
#ifdef FLA_OPENMP_MULTITHREADING
#if FLA_ENABLE_AOCL_BLAS
    /* Set no. of threads to BLIS as 1 to run DGEMV in ST.
     * This is to avoid isolated threading causing cache misses.
     * Done after the quick-return check so the original BLIS thread
     * count is never modified (and thus never needs restoring) on the
     * early-exit path.
     */
    aocl_int64_t orig_blis_threads = bli_thread_get_num_threads();
    bli_thread_set_num_threads(1);
#endif
#endif
    if(*m >= *n)
    {
        /* Reduce to upper bidiagonal form */
        i__1 = *nb;
#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp parallel num_threads(requested_num_threads) private( \
        i__, i__2, i__3, i__4, i__5, thread_id) shared(actual_num_threads, optimal_num_threads)
        {
            thread_id = omp_get_thread_num();
            /* Initialize barrier with actual team size inside parallel region */
            if(thread_id == 0)
            {
                actual_num_threads = omp_get_num_threads();
                FLA_BARRIER_INIT(*barrier, actual_num_threads);
                /* Chunk partitioning uses optimal_num_threads as divisor; keep >= 1 when actual
                 * < 8. */
                optimal_num_threads = fla_max(1, actual_num_threads / 8);
                if(optimal_num_threads <= 2)
                {
                    optimal_num_threads = actual_num_threads;
                }
            }
            /* Ensure all threads see the initialized barrier before proceeding */
#pragma omp barrier
#else
        {
            thread_id = 0;
#endif
            for(i__ = 1; i__ <= i__1; ++i__)
            {

                /* Update A(i:m,i) */
                i__2 = *m - i__ + 1;
                i__3 = i__ - 1;
#ifdef FLA_OPENMP_MULTITHREADING
                FLA_Thread_get_subrange(thread_id, actual_num_threads, i__2, &i__4, &i__5);
                FLA_BARRIER_WAIT(*barrier);
                {
                    aocl_blas_dgemv("No transpose", &i__4, &i__3, &neg_one, &a[i__ + a_dim1 + i__5],
                                    lda, &y[i__ + y_dim1], ldy, &d_one,
                                    &a[i__ + i__ * a_dim1 + i__5], &i_one);

                    aocl_blas_dgemv("No transpose", &i__4, &i__3, &neg_one, &x[i__ + x_dim1 + i__5],
                                    ldx, &a[i__ * a_dim1 + 1], &i_one, &d_one,
                                    &a[i__ + i__ * a_dim1 + i__5], &i_one);
                }
                FLA_BARRIER_WAIT(*barrier);
#else
                {
                    aocl_blas_dgemv("No transpose", &i__2, &i__3, &neg_one, &a[i__ + a_dim1], lda,
                                    &y[i__ + y_dim1], ldy, &d_one, &a[i__ + i__ * a_dim1], &i_one);

                    aocl_blas_dgemv("No transpose", &i__2, &i__3, &neg_one, &x[i__ + x_dim1], ldx,
                                    &a[i__ * a_dim1 + 1], &i_one, &d_one, &a[i__ + i__ * a_dim1],
                                    &i_one);
                }
#endif
                if(thread_id == 0)
                {
                    /* Generate reflection Q(i) to annihilate A(i+1:m,i) */
                    i__2 = *m - i__ + 1;
                    /* Computing MIN */
                    i__3 = i__ + 1;
                    aocl_lapack_dlarfg(&i__2, &a[i__ + i__ * a_dim1],
                                       &a[fla_min(i__3, *m) + i__ * a_dim1], &i_one, &tauq[i__]);
                    d__[i__] = a[i__ + i__ * a_dim1];
                }
                if(i__ < *n)
                {
                    if(thread_id == 0)
                    {
                        a[i__ + i__ * a_dim1] = 1.;
                    }
                    /* Compute Y(i+1:n,i) */
                    i__2 = *m - i__ + 1;
                    i__3 = *n - i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, i__3, &i__4, &i__5);
                    FLA_BARRIER_WAIT(*barrier);
                    {
                        aocl_blas_dgemv("Transpose", &i__2, &i__4, &d_one,
                                        &a[i__ + (i__5 + i__ + 1) * a_dim1], lda,
                                        &a[i__ + i__ * a_dim1], &i_one, &d_zero,
                                        &y[i__5 + i__ + 1 + i__ * y_dim1], &i_one);
                    }
#pragma omp single nowait /* nowait is safe here since there is a barrier before next gemv \
                             operation. */
                    {
                        /* This is gemv6, since all the inputs are ready, calculated and stored in
                         * extra buffer */
                        i__2 = *m - i__ + 1;
                        i__3 = i__ - 1;
                        aocl_blas_dgemv("Transpose", &i__2, &i__3, &d_one, &x[i__ + x_dim1], ldx,
                                        &a[i__ + i__ * a_dim1], &i_one, &d_zero,
                                        &gemv_a_row_buffer[(*nb)], &i_one);
                        aocl_blas_dcopy(&i__, &a[i__ + a_dim1], lda, &gemv_a_row_buffer[0], &i_one);
                    }
#else
                    aocl_blas_dgemv("Transpose", &i__2, &i__3, &d_one, &a[i__ + (i__ + 1) * a_dim1],
                                    lda, &a[i__ + i__ * a_dim1], &i_one, &d_zero,
                                    &y[i__ + 1 + i__ * y_dim1], &i_one);
#endif
                    i__2 = *m - i__ + 1;
                    i__3 = i__ - 1;
#ifdef FLA_OPENMP_MULTITHREADING
                    FLA_BARRIER_WAIT(*barrier);
                    if(thread_id < optimal_num_threads)
                    {
                        FLA_Thread_get_subrange(thread_id, optimal_num_threads, i__3, &i__4, &i__5);
                        aocl_blas_dgemv("Transpose", &i__2, &i__4, &d_one,
                                        &a[i__ + (i__5 + 1) * a_dim1], lda, &a[i__ + i__ * a_dim1],
                                        &i_one, &d_zero, &y[i__5 + i__ * y_dim1 + 1], &i_one);
                        /* Using extra buffer to store the result of gemv operation which can be
                         * reused in next iteration*/
                        memcpy(&gemv_a_row_buffer[2 * (*nb) + i__5], &y[i__5 + i__ * y_dim1 + 1],
                               i__4 * sizeof(double));
                        memcpy(&y[i__5 + i__ * y_dim1 + 1], &gemv_a_row_buffer[(*nb) + i__5],
                               i__4 * sizeof(double));
                    }
                    FLA_BARRIER_WAIT(*barrier);
#else
                    aocl_blas_dgemv("Transpose", &i__2, &i__3, &d_one, &a[i__ + a_dim1], lda,
                                    &a[i__ + i__ * a_dim1], &i_one, &d_zero, &y[i__ * y_dim1 + 1],
                                    &i_one);
#endif
                    i__2 = *n - i__;
                    i__3 = i__ - 1;
#ifdef FLA_OPENMP_MULTITHREADING
                    if(thread_id < optimal_num_threads)
                    {
                        FLA_Thread_get_subrange(thread_id, optimal_num_threads, i__2, &i__4, &i__5);
                        aocl_blas_dgemv("No transpose", &i__4, &i__3, &neg_one,
                                        &y[i__ + 1 + y_dim1 + i__5], ldy,
                                        &gemv_a_row_buffer[2 * (*nb)], &i_one, &d_one,
                                        &y[i__ + 1 + i__ * y_dim1 + i__5], &i_one);
                    }

#else
                    aocl_blas_dgemv("No transpose", &i__2, &i__3, &neg_one, &y[i__ + 1 + y_dim1],
                                    ldy, &y[i__ * y_dim1 + 1], &i_one, &d_one,
                                    &y[i__ + 1 + i__ * y_dim1], &i_one);
#endif
                    i__2 = i__ - 1;
                    i__3 = *n - i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    if(thread_id < optimal_num_threads)
                    {
                        FLA_Thread_get_subrange(thread_id, optimal_num_threads, i__3, &i__4, &i__5);
                        aocl_blas_dgemv("Transpose", &i__2, &i__4, &neg_one,
                                        &a[(i__ + i__5 + 1) * a_dim1 + 1], lda,
                                        &gemv_a_row_buffer[(*nb)], &i_one, &d_one,
                                        &y[i__5 + i__ + 1 + i__ * y_dim1], &i_one);
                        aocl_blas_dscal(&i__4, &tauq[i__], &y[i__5 + i__ + 1 + i__ * y_dim1],
                                        &i_one);
                    }
#else
                    aocl_blas_dgemv("Transpose", &i__2, &i__3, &neg_one, &a[(i__ + 1) * a_dim1 + 1],
                                    lda, &y[i__ * y_dim1 + 1], &i_one, &d_one,
                                    &y[i__ + 1 + i__ * y_dim1], &i_one);
                    aocl_blas_dscal(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &i_one);

#endif
                    /* Update A(i,i+1:n) */
                    i__2 = *n - i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    FLA_BARRIER_WAIT(*barrier);
                    if(thread_id < optimal_num_threads)
                    {
                        FLA_Thread_get_subrange(thread_id, optimal_num_threads, i__2, &i__4, &i__5);
                        aocl_blas_dgemv("No transpose", &i__4, &i__, &neg_one,
                                        &y[i__ + 1 + i__5 + y_dim1], ldy, &gemv_a_row_buffer[0],
                                        &i_one, &d_one, &a[i__ + (i__ + 1 + i__5) * a_dim1], lda);
                    }

#else
                    {
                        aocl_blas_dgemv("No transpose", &i__2, &i__, &neg_one, &y[i__ + 1 + y_dim1],
                                        ldy, &a[i__ + a_dim1], lda, &d_one,
                                        &a[i__ + (i__ + 1) * a_dim1], lda);
                    }
#endif
                    i__2 = i__ - 1;
                    i__3 = *n - i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    if(thread_id < optimal_num_threads)
                    {
                        FLA_Thread_get_subrange(thread_id, optimal_num_threads, i__3, &i__4, &i__5);
                        aocl_blas_dgemv("Transpose", &i__2, &i__4, &neg_one,
                                        &a[(i__ + 1 + i__5) * a_dim1 + 1], lda, &x[i__ + x_dim1],
                                        ldx, &d_one, &a[i__ + (i__ + 1 + i__5) * a_dim1], lda);
                    }
                    FLA_BARRIER_WAIT(*barrier);

#else
                    {
                        aocl_blas_dgemv("Transpose", &i__2, &i__3, &neg_one,
                                        &a[(i__ + 1) * a_dim1 + 1], lda, &x[i__ + x_dim1], ldx,
                                        &d_one, &a[i__ + (i__ + 1) * a_dim1], lda);
                    }
#endif
                    if(thread_id == 0)
                    {
                        /* Generate reflection P(i) to annihilate A(i,i+2:n) */
                        i__2 = *n - i__;
                        /* Computing MIN */
                        i__3 = i__ + 2;
                        aocl_lapack_dlarfg(&i__2, &a[i__ + (i__ + 1) * a_dim1],
                                           &a[i__ + fla_min(i__3, *n) * a_dim1], lda, &taup[i__]);
                        e[i__] = a[i__ + (i__ + 1) * a_dim1];
                        a[i__ + (i__ + 1) * a_dim1] = 1.;

                        i__3 = *n - i__;
                        if(gemv_a_row_buffer != NULL)
                        {
                            aocl_blas_dcopy(&i__3, &a[i__ + (i__ + 1) * a_dim1], lda,
                                            gemv_a_row_buffer, &i_one);
                        }
                    }
                    /* Compute X(i+1:m,i) */
                    i__2 = *m - i__;
                    i__3 = *n - i__;
#ifdef FLA_OPENMP_MULTITHREADING

                    FLA_Thread_get_subrange(thread_id, actual_num_threads, i__2, &i__4, &i__5);
                    FLA_BARRIER_WAIT(*barrier);
                    {
                        if(gemv_a_row_buffer != NULL)
                        {
                            aocl_blas_dgemv("No transpose", &i__4, &i__3, &d_one,
                                            &a[i__5 + i__ + 1 + (i__ + 1) * a_dim1], lda,
                                            &gemv_a_row_buffer[0], &i_one, &d_zero,
                                            &x[i__5 + i__ + 1 + i__ * x_dim1], &i_one);
                        }
                        else
                        {
                            aocl_blas_dgemv("No transpose", &i__4, &i__3, &d_one,
                                            &a[i__5 + i__ + 1 + (i__ + 1) * a_dim1], lda,
                                            &a[i__ + (i__ + 1) * a_dim1], lda, &d_zero,
                                            &x[i__5 + i__ + 1 + i__ * x_dim1], &i_one);
                        }
                    }

#else
                    aocl_blas_dgemv("No transpose", &i__2, &i__3, &d_one,
                                    &a[i__ + 1 + (i__ + 1) * a_dim1], lda,
                                    &a[i__ + (i__ + 1) * a_dim1], lda, &d_zero,
                                    &x[i__ + 1 + i__ * x_dim1], &i_one);
#endif

#ifdef FLA_OPENMP_MULTITHREADING

                    i__2 = *n - i__;
                    if(thread_id < optimal_num_threads)
                    {
                        FLA_Thread_get_subrange(thread_id, optimal_num_threads, i__, &i__4, &i__5);
                        aocl_blas_dgemv("Transpose", &i__2, &i__4, &d_one,
                                        &y[i__ + 1 + (i__5 + 1) * y_dim1], ldy,
                                        &a[i__ + (i__ + 1) * a_dim1], lda, &d_zero,
                                        &x[i__ * x_dim1 + 1 + i__5], &i_one);
                        memcpy(&gemv_a_row_buffer[(*n) + i__5], &x[i__ * x_dim1 + 1 + i__5],
                               i__4 * sizeof(double));
                    }
                    FLA_BARRIER_WAIT(*barrier);

#else
                    {
                        i__2 = *n - i__;
                        aocl_blas_dgemv("Transpose", &i__2, &i__, &d_one, &y[i__ + 1 + y_dim1], ldy,
                                        &a[i__ + (i__ + 1) * a_dim1], lda, &d_zero,
                                        &x[i__ * x_dim1 + 1], &i_one);
                    }
#endif
                    i__2 = *m - i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    if(thread_id < optimal_num_threads)
                    {
                        FLA_Thread_get_subrange(thread_id, optimal_num_threads, i__2, &i__4, &i__5);
                        aocl_blas_dgemv("No transpose", &i__4, &i__, &neg_one,
                                        &a[i__ + 1 + a_dim1 + i__5], lda, &gemv_a_row_buffer[(*n)],
                                        &i_one, &d_one, &x[i__ + 1 + i__ * x_dim1 + i__5], &i_one);
                    }

#else
                    {
                        aocl_blas_dgemv("No transpose", &i__2, &i__, &neg_one, &a[i__ + 1 + a_dim1],
                                        lda, &x[i__ * x_dim1 + 1], &i_one, &d_one,
                                        &x[i__ + 1 + i__ * x_dim1], &i_one);
                    }
#endif
                    i__2 = i__ - 1;
                    i__3 = *n - i__;

#ifdef FLA_OPENMP_MULTITHREADING
                    if(thread_id < optimal_num_threads)
                    {
                        FLA_Thread_get_subrange(thread_id, optimal_num_threads, i__2, &i__4, &i__5);
                        aocl_blas_dgemv("No transpose", &i__4, &i__3, &d_one,
                                        &a[(i__ + 1) * a_dim1 + 1 + i__5], lda,
                                        &a[i__ + (i__ + 1) * a_dim1], lda, &d_zero,
                                        &x[i__ * x_dim1 + 1 + i__5], &i_one);
                    }

#else
                    {
                        aocl_blas_dgemv("No transpose", &i__2, &i__3, &d_one,
                                        &a[(i__ + 1) * a_dim1 + 1], lda,
                                        &a[i__ + (i__ + 1) * a_dim1], lda, &d_zero,
                                        &x[i__ * x_dim1 + 1], &i_one);
                    }
#endif

                    i__2 = *m - i__;
                    i__3 = i__ - 1;
#ifdef FLA_OPENMP_MULTITHREADING
                    FLA_BARRIER_WAIT(*barrier);
                    if(thread_id < optimal_num_threads)
                    {
                        FLA_Thread_get_subrange(thread_id, optimal_num_threads, i__2, &i__4, &i__5);
                        aocl_blas_dgemv("No transpose", &i__4, &i__3, &neg_one,
                                        &x[i__ + 1 + x_dim1 + i__5], ldx, &x[i__ * x_dim1 + 1],
                                        &i_one, &d_one, &x[i__ + 1 + i__ * x_dim1 + i__5], &i_one);
                        aocl_blas_dscal(&i__4, &taup[i__], &x[i__ + 1 + i__ * x_dim1 + i__5],
                                        &i_one);
                    }
#else
                    {
                        aocl_blas_dgemv("No transpose", &i__2, &i__3, &neg_one,
                                        &x[i__ + 1 + x_dim1], ldx, &x[i__ * x_dim1 + 1], &i_one,
                                        &d_one, &x[i__ + 1 + i__ * x_dim1], &i_one);
                        aocl_blas_dscal(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &i_one);
                    }
#endif
                }
                /* L10: */
            }
        }
    }
    else
    {
        /* Reduce to lower bidiagonal form */
        i__1 = *nb;
#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp parallel num_threads(requested_num_threads) private( \
        i__, i__2, i__3, i__4, i__5, thread_id) shared(actual_num_threads)
        {
            thread_id = omp_get_thread_num();
            /* Initialize barrier with actual team size inside parallel region.
             * OpenMP may grant fewer than requested_num_threads, so the barrier
             * must be sized to the team that will actually use it. */
            if(thread_id == 0)
            {
                actual_num_threads = omp_get_num_threads();
                FLA_BARRIER_INIT(*barrier, actual_num_threads);
            }
            /* Ensure all threads see the initialized barrier before proceeding */
#pragma omp barrier
#else
        {
            thread_id = 0;
#endif
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                if(thread_id == 0)
                {
                    /* Update A(i,i:n) */
                    i__2 = *n - i__ + 1;
                    i__3 = i__ - 1;
                    aocl_blas_dgemv("No transpose", &i__2, &i__3, &neg_one, &y[i__ + y_dim1], ldy,
                                    &a[i__ + a_dim1], lda, &d_one, &a[i__ + i__ * a_dim1], lda);
                    i__2 = i__ - 1;
                    i__3 = *n - i__ + 1;
                    aocl_blas_dgemv("Transpose", &i__2, &i__3, &neg_one, &a[i__ * a_dim1 + 1], lda,
                                    &x[i__ + x_dim1], ldx, &d_one, &a[i__ + i__ * a_dim1], lda);
                    /* Generate reflection P(i) to annihilate A(i,i+1:n) */
                    i__2 = *n - i__ + 1;
                    /* Computing MIN */
                    i__3 = i__ + 1;
                    aocl_lapack_dlarfg(&i__2, &a[i__ + i__ * a_dim1],
                                       &a[i__ + fla_min(i__3, *n) * a_dim1], lda, &taup[i__]);
                    d__[i__] = a[i__ + i__ * a_dim1];
                }
                if(i__ < *m)
                {
                    if(thread_id == 0)
                    {
                        a[i__ + i__ * a_dim1] = 1.;
                    }
                    /* Compute X(i+1:m,i) */
                    i__2 = *m - i__;
                    i__3 = *n - i__ + 1;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, i__2, &i__4, &i__5);
                    FLA_BARRIER_WAIT(*barrier);
                    aocl_blas_dgemv("No transpose", &i__4, &i__3, &d_one,
                                    &a[i__5 + i__ + 1 + i__ * a_dim1], lda, &a[i__ + i__ * a_dim1],
                                    lda, &d_zero, &x[i__5 + i__ + 1 + i__ * x_dim1], &i_one);
                    FLA_BARRIER_WAIT(*barrier);
#else
                    aocl_blas_dgemv("No transpose", &i__2, &i__3, &d_one,
                                    &a[i__ + 1 + i__ * a_dim1], lda, &a[i__ + i__ * a_dim1], lda,
                                    &d_zero, &x[i__ + 1 + i__ * x_dim1], &i_one);
#endif
                    if(thread_id == 0)
                    {
                        i__2 = *n - i__ + 1;
                        i__3 = i__ - 1;
                        aocl_blas_dgemv("Transpose", &i__2, &i__3, &d_one, &y[i__ + y_dim1], ldy,
                                        &a[i__ + i__ * a_dim1], lda, &d_zero, &x[i__ * x_dim1 + 1],
                                        &i_one);
                        i__2 = *m - i__;
                        i__3 = i__ - 1;
                        aocl_blas_dgemv("No transpose", &i__2, &i__3, &neg_one,
                                        &a[i__ + 1 + a_dim1], lda, &x[i__ * x_dim1 + 1], &i_one,
                                        &d_one, &x[i__ + 1 + i__ * x_dim1], &i_one);
                        i__2 = i__ - 1;
                        i__3 = *n - i__ + 1;
                        aocl_blas_dgemv("No transpose", &i__2, &i__3, &d_one, &a[i__ * a_dim1 + 1],
                                        lda, &a[i__ + i__ * a_dim1], lda, &d_zero,
                                        &x[i__ * x_dim1 + 1], &i_one);
                        i__2 = *m - i__;
                        i__3 = i__ - 1;
                        aocl_blas_dgemv("No transpose", &i__2, &i__3, &neg_one,
                                        &x[i__ + 1 + x_dim1], ldx, &x[i__ * x_dim1 + 1], &i_one,
                                        &d_one, &x[i__ + 1 + i__ * x_dim1], &i_one);
                        i__2 = *m - i__;
                        aocl_blas_dscal(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &i_one);
                        /* Update A(i+1:m,i) */
                        i__2 = *m - i__;
                        i__3 = i__ - 1;
                        aocl_blas_dgemv("No transpose", &i__2, &i__3, &neg_one,
                                        &a[i__ + 1 + a_dim1], lda, &y[i__ + y_dim1], ldy, &d_one,
                                        &a[i__ + 1 + i__ * a_dim1], &i_one);
                        i__2 = *m - i__;
                        aocl_blas_dgemv("No transpose", &i__2, &i__, &neg_one, &x[i__ + 1 + x_dim1],
                                        ldx, &a[i__ * a_dim1 + 1], &i_one, &d_one,
                                        &a[i__ + 1 + i__ * a_dim1], &i_one);
                        /* Generate reflection Q(i) to annihilate A(i+2:m,i) */
                        i__2 = *m - i__;
                        /* Computing MIN */
                        i__3 = i__ + 2;
                        aocl_lapack_dlarfg(&i__2, &a[i__ + 1 + i__ * a_dim1],
                                           &a[fla_min(i__3, *m) + i__ * a_dim1], &i_one,
                                           &tauq[i__]);
                        e[i__] = a[i__ + 1 + i__ * a_dim1];
                        a[i__ + 1 + i__ * a_dim1] = 1.;
                    }
                    /* Compute Y(i+1:n,i) */
                    i__2 = *m - i__;
                    i__3 = *n - i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, i__3, &i__4, &i__5);
                    FLA_BARRIER_WAIT(*barrier);
                    aocl_blas_dgemv("Transpose", &i__2, &i__4, &d_one,
                                    &a[i__ + 1 + (i__5 + i__ + 1) * a_dim1], lda,
                                    &a[i__ + 1 + i__ * a_dim1], &i_one, &d_zero,
                                    &y[i__5 + i__ + 1 + i__ * y_dim1], &i_one);
                    FLA_BARRIER_WAIT(*barrier);
#else
                    aocl_blas_dgemv("Transpose", &i__2, &i__3, &d_one,
                                    &a[i__ + 1 + (i__ + 1) * a_dim1], lda,
                                    &a[i__ + 1 + i__ * a_dim1], &i_one, &d_zero,
                                    &y[i__ + 1 + i__ * y_dim1], &i_one);
#endif
                    if(thread_id == 0)
                    {
                        i__2 = *m - i__;
                        i__3 = i__ - 1;
                        aocl_blas_dgemv("Transpose", &i__2, &i__3, &d_one, &a[i__ + 1 + a_dim1],
                                        lda, &a[i__ + 1 + i__ * a_dim1], &i_one, &d_zero,
                                        &y[i__ * y_dim1 + 1], &i_one);
                        i__2 = *n - i__;
                        i__3 = i__ - 1;
                        aocl_blas_dgemv("No transpose", &i__2, &i__3, &neg_one,
                                        &y[i__ + 1 + y_dim1], ldy, &y[i__ * y_dim1 + 1], &i_one,
                                        &d_one, &y[i__ + 1 + i__ * y_dim1], &i_one);
                        i__2 = *m - i__;
                        aocl_blas_dgemv("Transpose", &i__2, &i__, &d_one, &x[i__ + 1 + x_dim1], ldx,
                                        &a[i__ + 1 + i__ * a_dim1], &i_one, &d_zero,
                                        &y[i__ * y_dim1 + 1], &i_one);
                        i__2 = *n - i__;
                        aocl_blas_dgemv("Transpose", &i__, &i__2, &neg_one,
                                        &a[(i__ + 1) * a_dim1 + 1], lda, &y[i__ * y_dim1 + 1],
                                        &i_one, &d_one, &y[i__ + 1 + i__ * y_dim1], &i_one);
                        i__2 = *n - i__;
                        aocl_blas_dscal(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &i_one);
                    }
                }
                /* L20: */
            }
        }
    }
#ifdef FLA_OPENMP_MULTITHREADING
#if FLA_ENABLE_AOCL_BLAS
    /* reset no. of threads back to original for BLIS */
    bli_thread_set_num_threads(orig_blis_threads);
#endif
#endif
    return;
    /* End of DLABRD */
}
#endif