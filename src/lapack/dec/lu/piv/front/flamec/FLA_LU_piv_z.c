/*
    Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
#include "FLA_f2c.h"
#include "fla_lapack_x86_common.h"

#define FLA_LU_SMALL_BLOCK_SIZE 4096
#define FLA_LU_SMALL_DIM 32

static dcomplex z__1 = { -1, 0};
static dcomplex c_b1 = {1.,0.};
static integer c__1 = 1;

#ifdef FLA_ENABLE_AMD_OPT

void FLA_get_optimum_params_zgetrf(integer m, integer n, integer *nb, int *n_threads)
{
    int available_n_threads;
    extern int fla_thread_get_num_threads(void);

    /* Get maximum thread available*/
    available_n_threads = fla_thread_get_num_threads();

#ifdef FLA_OPENMP_MULTITHREADING

    if(m <= 100 || n <= 100)
    {
        *nb = 15;
        *n_threads = 4;
    }
    else if(m <= 512 || n <= 512)
    {
        *nb = 15;
        *n_threads = 8;
    }
    else if(m <= 920 || n <= 920)
    {
        *nb = 15;
        *n_threads = 16;
    }
    else if(m <= 2048 || n <= 2048)
    {
        *nb = 15;
        *n_threads = 32;
    }
    else if(m <= 6144 || n <= 6144)
    {
        *nb = 15;
        *n_threads = 96;
    }
    else if(m <= 12288 || n <= 12288)
    {
        *nb = 32;
        *n_threads = 96;
    }
    else
    {
        *nb = 64;
        *n_threads = 192;
    }

    if(*n_threads > available_n_threads)
        *n_threads = available_n_threads;

#else
    *nb = 64;
    *n_threads = 1;
#endif

    return;
}

/*
 * LU with partial pivoting for tiny matrices
 *
 * All the computations are done inline without using
 * corresponding BLAS APIs to reduce function overheads.
 */
int FLA_LU_piv_small_z_var0( integer *m, integer *n,
                                   dcomplex *a, integer *lda,
                                   integer *ipiv,
                                   integer *info)
{
    integer mi, ni;
    integer i, j, i_1, i_2, i_3;
    double max_val, t_val, z_val;
    dcomplex *acur, *apiv, *asrc;
    integer p_idx;
    integer min_m_n = fla_min(*m, *n);
#ifndef _WIN32
    dcomplex z__1;
    double _Complex pinv;
#else
    double piv_r, piv_i;
    double pinv;
#endif

    *info = 0;

    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -4;
    }

    if (*info != 0)
    {
        return 0;
    }

    for( i = 0; i < min_m_n; i++ )
    {
        mi = *m - i;
        ni = *n - i;

        acur = &a[i + *lda * i];

        /* Find the pivot element */
        max_val = 0;
        p_idx = i;
        for( i_1 = 0; i_1 < mi; i_1++ )
        {
            t_val = f2c_abs(acur[i_1].real) + f2c_abs(acur[i_1].imag);
            if( t_val > max_val )
            {
                max_val = t_val;
                p_idx = i + i_1;
            }
        }

        apiv = a + p_idx;
        asrc = a + i;
        ipiv[i] = p_idx + 1;

        /* Swap rows and calculate a column of L */
        if( apiv[*lda * i].real != 0. || apiv[*lda * i].imag != 0. )
        {
            /* Swap entire rows */
            if( p_idx != i )
            {
                for( i_1 = 0; i_1 < *n ; i_1++ )
                {
                    i_2 = i_1 * *lda;
                    t_val = apiv[i_2].real;
                    z_val = apiv[i_2].imag;
                    apiv[i_2].real = asrc[i_2].real;
                    apiv[i_2].imag = asrc[i_2].imag;
                    asrc[i_2].real = t_val;
                    asrc[i_2].imag = z_val;
                 }
            }

            /* Calculate scalefactors (L) & update trailing matrix */

#ifndef _WIN32
            pinv = 1.0 / ((*acur).real + (I * (*acur).imag));
            z__1.real = creal(pinv);
            z__1.imag = cimag(pinv);
#else
            piv_r = (*acur).real;
            piv_i = (*acur).imag;
            pinv = piv_r * piv_r + piv_i * piv_i;
#endif

            for( i_1 = 1; i_1 < mi; i_1++ )
            {
                t_val = acur[i_1].real;
#ifndef _WIN32
                acur[i_1].real = (t_val * z__1.real - acur[i_1].imag * z__1.imag);
                acur[i_1].imag = (t_val * z__1.imag + acur[i_1].imag * z__1.real);
#else
                acur[i_1].real = (acur[i_1].imag * piv_i + t_val * piv_r) / pinv;
                acur[i_1].imag = (acur[i_1].imag * piv_r - t_val * piv_i) / pinv;
#endif
                t_val = acur[i_1].real;
                z_val = acur[i_1].imag;

                for( j = 1; j < ni; j++ )
                {
                    i_2 = i_1 + j * *lda;
                    i_3 = j * *lda;

                    acur[i_2].real = acur[i_2].real - t_val * acur[i_3].real + z_val * acur[i_3].imag;
                    acur[i_2].imag = acur[i_2].imag - t_val * acur[i_3].imag - z_val * acur[i_3].real;
                }
            }
        }
        else
        {
            *info = ( *info == 0 ) ? p_idx + 1 : *info;
        }
    }

    return *info;
}


/* LU factorization recursive variant*/
int FLA_LU_piv_z_var0(integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info)
{
    integer a_dim1, i__1, i__2, i__, n1, n2, block_size;
    integer iinfo;

    /* Adjust dimension of the matrix */
    a_dim1 = *lda;

    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -4;
    }

    if (*info != 0)
    {
        return 0;
    }

    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        return 0;
    }

    /* compute matrix size*/
    block_size = *m * *n;

    if ( (block_size <= FLA_LU_SMALL_BLOCK_SIZE) && ( *m <= FLA_LU_SMALL_DIM || *n <= FLA_LU_SMALL_DIM ) )
    {
        fla_zgetrf_small_simd(m, n, a, lda, ipiv, &iinfo);

        if (*info == 0 && iinfo > 0)
        {
            *info = iinfo;
        }
    }
    else if (*m == 1 || *n == 1)
    {
        lapack_zgetf2(m, n, a, lda, ipiv, &iinfo);

        if (*info == 0 && iinfo > 0)
        {
            *info = iinfo;
        }
    }
    else
    {
        /* Use recursive code */

        /* calculate n1 and n2 for recursive call*/
        n1 = fla_min(*m,*n) / 2;
        n2 = *n - n1;

        /*        [ A11 ] */
        /* Factor [ --- ] */
        /*        [ A21 ] */
        FLA_LU_piv_z_var0(m, &n1, a, lda, ipiv, &iinfo);

        /* Update info if necessary*/
        if (*info == 0 && iinfo > 0)
        {
            *info = iinfo;
        }

        /*                       [ A12 ] */
        /* Apply interchanges to [ --- ] */
        /*                       [ A22 ] */
        zlaswp_(&n2, &a[(n1 ) * a_dim1 ], lda, &c__1, &n1, ipiv, & c__1);

        /* Solve A12 */
        ztrsm_("L", "L", "N", "U", &n1, &n2, &c_b1, a, lda, &a[(n1 ) * a_dim1], lda);

        /* Update A22 */
        i__1 = *m - n1;
        zgemm_("N", "N", &i__1, &n2, &n1, &z__1, &a[n1], lda, &a[ (n1) * a_dim1 ], lda, &c_b1, &a[n1 +  (n1) * a_dim1], lda);

        /* Factor A22 */
        i__1 = *m - n1;
        FLA_LU_piv_z_var0(&i__1, &n2, &a[n1 + (n1) * a_dim1], lda, &ipiv[n1], &iinfo);

        /* Adjust INFO and the pivot indices */
        if (*info == 0 && iinfo > 0)
        {
            *info = iinfo + n1;
        }

        /* Update pivot indices*/
        i__1 = fla_min(*m,*n);
        for (i__ = n1; i__ < i__1; ++i__)
        {
            ipiv[i__] += n1;
        }

        /* Apply interchanges to A21 */
        i__1 = n1 + 1;
        i__2 = fla_min(*m,*n);
        zlaswp_(&n1, a, lda, &i__1, &i__2, ipiv, &c__1);
    }

    return 0;
}

#ifdef FLA_OPENMP_MULTITHREADING

/* LU factorization blocked varaiant */
int FLA_LU_piv_z_var1_parallel( integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    integer i__, j, s, iinfo;
    integer jb, jb_prev, jb_offset, nb;
    integer c__1 = 1;
    #define a_subscr(a_1,a_2) (a_2)*a_dim1 + a_1
    #define a_ref(a_1,a_2) a[a_subscr(a_1,a_2)]
    int threads_id, n_threads;

    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -4;
    }

    if (*info != 0)
    {
        return *info;
    }

    // Quick return if possible
    if (*m == 0 || *n == 0)
        return 0;

    // Determine optimum block and thread size for this environment
    FLA_get_optimum_params_zgetrf(*m, *n, &nb, &n_threads);

    /* call sequencial algorithm for single thread*/
    if(n_threads == 1)
    {
        FLA_LU_piv_z_var0( m, n, a, lda, ipiv, info);
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
    i__1 = fla_min(*m,*n);
    i__2 = nb;

    i__3 = fla_min(*m,*n) - j + 1;
    jb = fla_min(i__3,nb);

    // Compute L00 and U00 of diagonal blocks
    i__3 = *m - j + 1;
    FLA_LU_piv_z_var0(&i__3, &jb, (dcomplex *) &a_ref(j, j), lda, &ipiv[j], &iinfo);

    if (*info == 0 && iinfo > 0)
        *info = iinfo + j - 1;

    // Computing MIN
    i__4 = *m, i__5 = j + jb - 1;
    i__3 = fla_min(i__4,i__5);
    for (i__ = j; i__ <= i__3; ++i__)
    {
        ipiv[i__] = j - 1 + ipiv[i__];
    }

    #pragma omp parallel num_threads(n_threads) private(i__3, i__4, i__5, i__6, i__7, i__8, j, s, jb, jb_prev, jb_offset, threads_id)
    {
        threads_id = omp_get_thread_num();
        for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
        {
            s = j + i__2;

            /* Factorize the next block while the rest of the matrix is being updated */
            if(threads_id == 0)
            {
                // Computing MIN
                i__3 = fla_min(*m,*n) - j + 1;
                jb_prev = fla_min(i__3,nb);

                // Apply interchanges to columns J+JB:N
                i__3 = *n - j - jb_prev + 1;
                i__4 = j + jb_prev - 1;
                i__3 = fla_min(i__3, jb_prev);
                zlaswp_(&i__3, (dcomplex*) &a_ref(1, j + jb_prev), lda, &j, &i__4, &ipiv[1], &c__1);

                // compute U10
                i__3 = *n - j - jb_prev + 1;
                i__3 = fla_min(i__3, jb_prev);
                ztrsm_("Left", "Lower", "No transpose", "Unit", &jb_prev, &i__3, &c_b1, (dcomplex*) &a_ref(j, j), lda, (dcomplex*) &a_ref(j, j + jb_prev), lda);

                // compute L11 * U11
                if (j + jb_prev <= *m)
                {
                    /* Update trailing submatrix. */
                    i__3 = *m - j - jb_prev + 1;
                    i__4 = *n - j - jb_prev + 1;
                    i__4 = fla_min(i__4, jb_prev);
                    zgemm_("No transpose", "No transpose", &i__3, &i__4, &jb_prev, &z__1,(dcomplex*)  &a_ref(j + jb_prev, j), lda,(dcomplex*) &a_ref(j, j + jb_prev), lda, &c_b1,(dcomplex*) &a_ref(j + jb_prev, j + jb_prev), lda);
                }

                if(s <= i__1)
                {

                    // Computing MIN
                    i__3 = fla_min(*m,*n) - s + 1;
                    jb = fla_min(i__3,nb);

                    // Compute L00 and U00 of diagonal blocks
                    i__3 = *m - s + 1;
                    FLA_LU_piv_z_var0(&i__3, &jb,(dcomplex*) &a_ref(s, s), lda, &ipiv[s], &iinfo);

                    if (*info == 0 && iinfo > 0)
                        *info = iinfo + s - 1;

                    // Computing MIN
                    i__4 = *m, i__5 = s + jb - 1;
                    i__3 = fla_min(i__4,i__5);
                    for (i__ = s; i__ <= i__3; ++i__)
                    {
                        ipiv[i__] = s - 1 + ipiv[i__];
                    }
                }
            }
            else
            {
                // Computing MIN
                i__3 = fla_min(*m,*n) - j + 1;
                jb_prev = fla_min(i__3,nb);
                jb_offset = jb_prev * 2;

                if (j + jb_prev <= *n)
                {
                    // Apply interchanges to columns J+JB:N
                    i__3 = *n - j - jb_prev + 1 - jb_prev;
                    i__4 = j + jb_prev - 1;
                    FLA_Thread_get_subrange(threads_id - 1, n_threads - 1, i__3, &i__5, &i__6);
                    zlaswp_(&i__5,(dcomplex*) &a_ref(1, j + jb_offset + i__6), lda, &j, &i__4, &ipiv[1], &c__1);

                    // compute U10
                    i__3 = *n - j - jb_prev + 1 - jb_prev;
                    FLA_Thread_get_subrange(threads_id - 1, n_threads - 1, i__3, &i__5, &i__6);
                    ztrsm_("Left", "Lower", "No transpose", "Unit", &jb_prev, &i__5, &c_b1,(dcomplex*) &a_ref(j, j), lda,(dcomplex*) &a_ref(j, j + jb_offset + i__6), lda);

                    // compute L11 * U11
                    if (j + jb_prev <= *m)
                    {
                        /* Update trailing submatrix. */
                        i__3 = *m - j - jb_prev + 1;
                        i__4 = *n - j - jb_prev + 1 - jb_prev;
                        FLA_Thread_get_subrange(threads_id - 1, n_threads - 1, i__4, &i__7, &i__8);
                        zgemm_("No transpose", "No transpose", &i__3, &i__7, &jb_prev, &z__1,(dcomplex*) &a_ref(j + jb_prev, j), lda,(dcomplex*) &a_ref(j, j + jb_offset + i__8), lda, &c_b1,(dcomplex*) &a_ref(j + jb_prev, j + jb_offset + i__8), lda);
                    }
                }
            }

            #pragma omp barrier
        }
    }

    #pragma omp parallel num_threads(n_threads) private(j, i__3, i__4, i__5, i__6, jb, threads_id)
    {
        threads_id = omp_get_thread_num();
        for (j = 1; j <= i__1; j += i__2)
        {
            // Computing MIN
            i__3 = fla_min(*m,*n) - j + 1;
            jb = fla_min(i__3,nb);

            i__3 = j - 1;
            i__4 = j + jb - 1;
            FLA_Thread_get_subrange(threads_id, n_threads, i__3, &i__5, &i__6);
            zlaswp_(&i__5,(dcomplex*) &a[a_offset + (i__6 * a_dim1)], lda, &j, &i__4, &ipiv[1], &c__1);
            #pragma omp barrier
        }
    }

    return *info;
}

#endif
#endif
