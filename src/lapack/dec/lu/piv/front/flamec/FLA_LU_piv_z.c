/*
    Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"
#include "FLA_f2c.h"


#define FLA_LU_SMALL_BLOCK_SIZE 16


void FLA_get_optimum_params_zgetrf(integer m, integer n, integer *nb, int *n_threads)
{
    int available_n_threads;

    /* Get maximum thread available*/
    available_n_threads = fla_thread_get_num_threads();

#ifdef FLA_OPENMP_MULTITHREADING
    if(m >= 2048 || n >= 2048)
        *nb = 128;
    else if(m >= 512 || n >= 512)
        *nb = 64;
    else
        *nb = 32;

    if(m >= 2048 || n >= 2048)
        *n_threads = 64;
    else if(m >= 512 || n >= 512)
        *n_threads = 24;
    else if(m >= 256 || n >= 256)
        *n_threads = 8;
    else
        *n_threads = 4;
#else
    *nb = 64;
    *n_threads = 1;
#endif

    if(*n_threads > available_n_threads)
        *n_threads = available_n_threads;

return;
}


/* LU factorization recursive variant*/
integer FLA_LU_piv_z_var0(integer *m, integer *n, doublecomplex *a, integer *lda, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, i__1, i__2;
    doublecomplex z__1;
    /* Local variables */
    integer i__, n1, n2;
    doublecomplex temp;
    integer iinfo;
    doublereal sfmin;
    doublecomplex c_b1 = {1.,0.};
    integer c__1 = 1;

    /* Parameter adjustments */
    a_dim1 = *lda;

    /* Function Body */
    *info = 0;

    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        return 0;
    }

    if (*m <= FLA_LU_SMALL_BLOCK_SIZE || *n <= FLA_LU_SMALL_BLOCK_SIZE)
    {
        FLA_LU_piv_small_z_avx2(m, n, a, lda, ipiv, &iinfo);
    }
    else
    {
        /* Use recursive code */
        n1 = fla_min(*m,*n) / 2;
        n2 = *n - n1;
        /* [ A11 ] */
        /* Factor [ --- ] */
        /* [ A21 ] */

        FLA_LU_piv_z_var0(m, &n1, a, lda, ipiv, &iinfo);

        if (*info == 0 && iinfo > 0)
        {
            *info = iinfo;
        }

        /* [ A12 ] */
        /* Apply interchanges to [ --- ] */
        /* [ A22 ] */
        zlaswp_(&n2, &a[(n1 ) * a_dim1 ], lda, &c__1, &n1, ipiv, & c__1);

        /* Solve A12 */
        ztrsm_("L", "L", "N", "U", &n1, &n2, &c_b1, a, lda, &a[(n1 ) * a_dim1], lda);

        /* Update A22 */
        i__1 = *m - n1;
        z__1.r = -1.;
        z__1.i = -0.; // , expr subst
        zgemm_("N", "N", &i__1, &n2, &n1, &z__1, &a[n1], lda, &a[ (n1) * a_dim1 ], lda, &c_b1, &a[n1 +  (n1) * a_dim1], lda);

        /* Factor A22 */
        i__1 = *m - n1;
        FLA_LU_piv_z_var0(&i__1, &n2, &a[n1 + (n1) * a_dim1], lda, &ipiv[n1], &iinfo);

        /* Adjust INFO and the pivot indices */
        if (*info == 0 && iinfo > 0)
        {
            *info = iinfo + n1;
        }

        i__1 = fla_min(*m,*n);
        for (i__ = n1; i__ < i__1; ++i__)
        {
            ipiv[i__] += n1;
            /* L20: */
        }

        /* Apply interchanges to A21 */
        i__1 = n1 + 1;
        i__2 = fla_min(*m,*n);
        zlaswp_(&n1, a, lda, &i__1, &i__2, ipiv, &c__1);
    }
    return 0;
}

/* LU factorization blocked varaiant */
integer FLA_LU_piv_z_var1( integer *m, integer *n,
                                   doublecomplex *a, integer *lda,
                                   integer *ipiv,
                                   integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10;
    doublecomplex z__1 = {-1, 0};
    integer i__, j, iinfo;
    integer jb, nb;
    doublecomplex c_b1 = {1.,0.};
    integer c__1 = 1;
    integer c_n1 = -1;
    integer x;
    #define a_subscr(a_1,a_2) (a_2)*a_dim1 + a_1
    #define a_ref(a_1,a_2) a[a_subscr(a_1,a_2)]
    int n_threads;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --ipiv;

    // Quick return if possible
    if (*m == 0 || *n == 0)
        return 0;

    // Determine optimum block and thread size for this environment
    FLA_get_optimum_params_zgetrf(*m, *n, &nb, &n_threads);

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

    i__1 = fla_min(*m,*n);
    i__2 = nb;

    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
    {
        // Computing MIN
        i__3 = fla_min(*m,*n) - j + 1;
        jb = fla_min(i__3,nb);

        // Compute L00 and U00 of diagonal blocks
        i__3 = *m - j + 1;
        FLA_LU_piv_z_var0(&jb, &jb, &a_ref(j, j), lda, &ipiv[j], &iinfo);

        if (*info == 0 && iinfo > 0)
            *info = iinfo + j - 1;

        // Computing MIN
        i__4 = *m, i__5 = j + jb - 1;
        i__3 = fla_min(i__4,i__5);
        for (i__ = j; i__ <= i__3; ++i__)
        {
            ipiv[i__] = j - 1 + ipiv[i__];
        }

        // Apply interchanges to rows J:J+JB
        i__3 = j - 1;
        i__4 = j + jb - 1;
        zlaswp_(&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);

        // compute L10
        if (j + jb <= *m)
        {
            i__3 = *m - j - jb + 1;
            ztrsm_("Right", "Upper", "No transpose", "N", &i__3, &jb, &c_b1, &a_ref(j, j), lda, &a_ref(j + jb, j), lda);
        }

        if (j + jb <= *n)
        {
            // Apply interchanges to columns J+JB:N
            i__3 = *n - j - jb + 1;
            i__4 = j + jb - 1;
            zlaswp_(&i__3, &a_ref(1, j + jb), lda, &j, &i__4, &ipiv[1], &c__1);

            // compute U10
            i__3 = *n - j - jb + 1;
            ztrsm_("Left", "Lower", "No transpose", "Unit", &jb, &i__3, &c_b1, &a_ref(j, j), lda, &a_ref(j, j + jb), lda);


            // compute L11 * U11
            if (j + jb <= *m)
            {
                /* Update trailing submatrix. */
                i__3 = *m - j - jb + 1;
                i__4 = *n - j - jb + 1;
                zgemm_("No transpose", "No transpose", &i__3, &i__4, &jb, &z__1, &a_ref(j + jb, j), lda, &a_ref(j, j + jb), lda, &c_b1, &a_ref(j + jb, j + jb), lda);
            }
        }
    }

    return *info;
}


#ifdef FLA_OPENMP_MULTITHREADING
/* LU factorization varaiant for mutli-threading environment */
integer FLA_LU_piv_z_var1_parallel( integer *m, integer *n,
                                   doublecomplex *a, integer *lda,
                                   integer *ipiv,
                                   integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10;
    doublecomplex z__1 = {-1, 0};
    integer i__, j, iinfo;
    integer jb, nb;
    doublecomplex c_b1 = {1.,0.};
    integer c__1 = 1;
    integer c_n1 = -1;
    integer x;
    #define a_subscr(a_1,a_2) (a_2)*a_dim1 + a_1
    #define a_ref(a_1,a_2) a[a_subscr(a_1,a_2)]
    int threads_id, n_threads;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --ipiv;

    // Quick return if possible
    if (*m == 0 || *n == 0)
        return 0;

    // Determine optimum block and thread size for this environment
    FLA_get_optimum_params_zgetrf(*m, *n, &nb, &n_threads);

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

    i__1 = fla_min(*m,*n);
    i__2 = nb;

    #pragma omp parallel num_threads(n_threads) private(i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, j, threads_id)
    {
        threads_id = omp_get_thread_num();
        for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
        {
            #pragma omp single
            {
                // Computing MIN
                i__3 = fla_min(*m,*n) - j + 1;
                jb = fla_min(i__3,nb);

                // Compute L00 and U00 of diagonal blocks
                i__3 = *m - j + 1;
                FLA_LU_piv_z_var0(&jb, &jb, &a_ref(j, j), lda, &ipiv[j], &iinfo);

                if (*info == 0 && iinfo > 0)
                    *info = iinfo + j - 1;

                // Computing MIN
                i__4 = *m, i__5 = j + jb - 1;
                i__3 = fla_min(i__4,i__5);
                for (i__ = j; i__ <= i__3; ++i__)
                {
                    ipiv[i__] = j - 1 + ipiv[i__];
                }
            }


            // Apply interchanges to rows J:J+JB
            i__3 = j - 1;
            i__4 = j + jb - 1;
            FLA_Thread_get_subrange(threads_id, n_threads, i__3, &i__5, &i__6);
            zlaswp_(&i__5, &a[a_offset + (i__6 * a_dim1)], lda, &j, &i__4, &ipiv[1], &c__1);

            // compute L10
            if (j + jb <= *m)
            {
                i__3 = *m - j - jb + 1;
                FLA_Thread_get_subrange(threads_id, n_threads, i__3, &i__5, &i__6);
                ztrsm_("Right", "Upper", "No transpose", "N", &i__5, &jb, &c_b1, &a_ref(j, j), lda, &a_ref(j + jb + i__6, j), lda);
            }

            if (j + jb <= *n)
            {
                // Apply interchanges to columns J+JB:N
                i__3 = *n - j - jb + 1;
                i__4 = j + jb - 1;
                FLA_Thread_get_subrange(threads_id, n_threads, i__3, &i__5, &i__6);
                zlaswp_(&i__5, &a_ref(1, j + jb + i__6), lda, &j, &i__4, &ipiv[1], &c__1);

                // compute U10
                i__3 = *n - j - jb + 1;
                FLA_Thread_get_subrange(threads_id, n_threads, i__3, &i__5, &i__6);
                ztrsm_("Left", "Lower", "No transpose", "Unit", &jb, &i__5, &c_b1, &a_ref(j, j), lda, &a_ref(j, j + jb + i__6), lda);

                #pragma omp barrier

                // compute L11 * U11
                if (j + jb <= *m)
                {
                    /* Update trailing submatrix. */
                    i__3 = *m - j - jb + 1;
                    i__4 = *n - j - jb + 1;
                    FLA_Thread_get_subrange(threads_id, n_threads, i__4, &i__7, &i__8);
                    zgemm_("No transpose", "No transpose", &i__3, &i__7, &jb, &z__1, &a_ref(j + jb, j), lda, &a_ref(j, j + jb + i__8), lda, &c_b1, &a_ref(j + jb, j + jb + i__8), lda);
                }
            }

            #pragma omp barrier
        }
    }

    return *info;
}
#endif