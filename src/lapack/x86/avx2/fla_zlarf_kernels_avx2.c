/******************************************************************************
 * * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT
/* Apply a left Householder reflector when the full vector v is stored with
   unit stride. This fuses the column dot products and rank-1 update for small
   complex-double panels using AVX2/FMA operations. */
void fla_zlarf_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff,
                                     aocl_int64_t ldr, dcomplex *v, dcomplex *ntau, dcomplex *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp1, vd2_dtmp2, vd2_vj, vd2_vjr, vd2_vji;
    __m256d vd4_inp, vd4_dtmp1, vd4_dtmp2, vd4_vj, vd4_vjr, vd4_vji;
    __m128d vd2_ltmp, vd2_htmp;
    __m128d vd2_one_neg_one = _mm_set_pd(-1.0, 1.0);

    /* Apply the Householder rotation                      */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */
    arows = m;
    acols = n;
    vd2_ntau = _mm_loadu_pd((const doublereal *)ntau);

    /* Compute A**T * v */
    for(j = 1; j <= acols; j++) /* for every column c_A of A */
    {
        vd2_dtmp1 = _mm_setzero_pd();
        vd2_dtmp2 = _mm_setzero_pd();
        vd4_dtmp1 = _mm256_setzero_pd();
        vd4_dtmp2 = _mm256_setzero_pd();
        /* Compute tmp = c_A**T . v */
        for(k = 1; k <= (arows - 1); k += 2)
        {
            /* load column elements of A and v */
            /* Column loaded in the following format */
            /* [Ar0, Ai0, Ar1, Ai1] */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            /* Vector v is loaded in following format */
            /* [Vr0, Vi0, Vr1, Vi1] */
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* Rearrange vd4 as follows */
            /* vd4_vjr = [Vr0, Vr0, Vr1, Vr1] */
            vd4_vjr = _mm256_permute_pd(vd4_vj, 0b0000);
            /* vd4_vji = [Vi0, Vi0, Vi1, Vi1] */
            vd4_vji = _mm256_permute_pd(vd4_vj, 0b1111);

            /* take dot product */
            vd4_dtmp1 = _mm256_fmadd_pd(vd4_inp, vd4_vjr, vd4_dtmp1);
            vd4_dtmp2 = _mm256_fmadd_pd(vd4_inp, vd4_vji, vd4_dtmp2);
        }
        if(k == arows)
        {
            /* load column elements of A and v */
            /* Column is loaded in following format */
            /* [Ar0, Ai0] */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            /* Vector c is loaded in following format */
            /* [Vr0, Vi0] */
            vd2_vj = _mm_loadu_pd((const doublereal *)&v[k]);

            /* Rearrange vd2 as follows */
            /* vd2vjr = [Vr0, Vr0] */
            vd2_vjr = _mm_permute_pd(vd2_vj, 0b00);
            /* vd2vji = [Vi0, Vi0] */
            vd2_vji = _mm_permute_pd(vd2_vj, 0b11);

            /* take dot product */
            vd2_dtmp1 = _mm_fmadd_pd(vd2_inp, vd2_vjr, vd2_dtmp1);
            vd2_dtmp2 = _mm_fmadd_pd(vd2_inp, vd2_vji, vd2_dtmp2);
            k += 1;
        }
        /* Reduce add the values in vd8_dtmp, vd4_dtmp and vd2_dtmp */

        /* Etract Upper and lower 128 bits of vd4_dtmp1 */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp1);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp1, 0x1);

        /* Add the lower and upper 128 bits and store in vd2_dtmp1 */
        vd2_dtmp1 = _mm_add_pd(vd2_dtmp1, vd2_ltmp);
        vd2_dtmp1 = _mm_add_pd(vd2_dtmp1, vd2_htmp);

        /* Etract Upper and lower 128 bits of vd4_dtmp2 */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp2);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp2, 0x1);

        /* Add the lower and upper 128 bits and store in vd2_dtmp2 */
        vd2_dtmp2 = _mm_add_pd(vd2_dtmp2, vd2_ltmp);
        vd2_dtmp2 = _mm_add_pd(vd2_dtmp2, vd2_htmp);

        /*
            Register vd2_dtmp1 = [ Ar * Vr, Ai * Vr ]
            Regsiter vd2_dtmp2 = [ Ar * Vi, Ai * Vi ]
            Since taking conjugate of A, the multiplcation result would
            bas as follows
            [ Ar * Vr + Ai * Vi, Ar * Vi - Ai * Vr ]
        */
        /* permuting vd2_dtmp2 = [ Ai * Vi, Ar * Vi ] */
        vd2_dtmp2 = _mm_permute_pd(vd2_dtmp2, 0b01);
        /* multiplying vd2_dtmp2 with [1.0, -1.0] */
        /* vd2_dtmp1 = [ Ar * Vr, - Ai * Vr ] */
        vd2_dtmp1 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);
        /* adding vd2_dtmp1 and vd2_dtmp2 */
        vd2_dtmp1 = _mm_add_pd(vd2_dtmp1, vd2_dtmp2);

        /* Store the result in work */
        _mm_storeu_pd((doublereal *)&work[j], vd2_dtmp1);

        /* Take conjugate of tmp */
        vd2_dtmp1 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);

        /* vd2_vjr = [ Tr, Tr ] */
        vd2_vjr = _mm_permute_pd(vd2_ntau, 0b00);
        /* vd2_vji = [ Ti. Ti ] */
        vd2_vji = _mm_permute_pd(vd2_ntau, 0b11);
        /*
            tmp = [ Kr, Ki ]
            tau = [ Tr, Ti ]
            Muliplication will be as follows
            [Kr * Tr - Ki * Ti, Kr * Ti + Ki * Tr]
        */
        /* vd2_temp = [ Kr, -Ki ] */
        vd2_dtmp2 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);
        /* Permute vd2_dtmp2 = [ -Ki, Kr ] */
        vd2_dtmp2 = _mm_permute_pd(vd2_dtmp2, 0b01);
        /* vd2_dtmp1 = [ Kr * Tr, Ki * Tr ] */
        vd2_dtmp1 = _mm_mul_pd(vd2_dtmp1, vd2_vjr);
        /* vd2_dtmp1 = [Kr * Tr - Ki * Ti, Kr * Ti + Ki * Tr] = [ Pr Pi ] */
        vd2_dtmp1 = _mm_fmadd_pd(vd2_vji, vd2_dtmp2, vd2_dtmp1);

        /* Set the first element as negative [ Pr -Pi ] */
        vd2_dtmp2 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);
        /* Shuffle vd2_dtmp2 = [ -Pi Pr ] */
        vd2_dtmp2 = _mm_permute_pd(vd2_dtmp2, 0b01);

        vd4_dtmp1 = _mm256_castpd128_pd256(vd2_dtmp1);
        vd4_dtmp1 = _mm256_insertf128_pd(vd4_dtmp1, vd2_dtmp1, 0x1);

        vd4_dtmp2 = _mm256_castpd128_pd256(vd2_dtmp2);
        vd4_dtmp2 = _mm256_insertf128_pd(vd4_dtmp2, vd2_dtmp2, 0x1);

        /* alternate for above 2 instructions which do not  */
        /* compile for older gcc versions (7 and below).    */
        /* Both will be same in terms of latency though     */
        /* vd4_dtmp = _mm256_set_m128d(vd2_dtmp, vd2_dtmp); */

        /* Compute c_A + tmp * v */
        for(k = 1; k <= (arows - 1); k += 2)
        {
            /* load column elements of c_A and v */
            /* vd4_inp = [ Ar0, Ai0, Ar1, Ai1 ] */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            /* vd4_vj = [ Vr0, Vi0, Vr1, Vi1 ] */
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* vd4_vjr = [ Vr0, Vr0, Vr1, Vr1 ] */
            vd4_vjr = _mm256_permute_pd(vd4_vj, 0b0000);
            /* vd4_vji = [ Vi0, Vi0, Vi1, Vi1 ] */
            vd4_vji = _mm256_permute_pd(vd4_vj, 0b1111);

            /* inp + [  Pr * Vr, Pi * Vr ] */
            vd4_inp = _mm256_fmadd_pd(vd4_dtmp1, vd4_vjr, vd4_inp);
            /* inp + [ Pr * Vr - Pi * Vi, Pi * Vr + Pr * Vi ] */
            vd4_inp = _mm256_fmadd_pd(vd4_dtmp2, vd4_vji, vd4_inp);

            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
        }
        if(k == arows)
        {
            /* load column elements of c_A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj = _mm_loadu_pd((const doublereal *)&v[k]);

            vd2_vjr = _mm_permute_pd(vd2_vj, 0b00);
            vd2_vji = _mm_permute_pd(vd2_vj, 0b11);

            /* mul by dtmp, add and store */
            vd2_inp = _mm_fmadd_pd(vd2_dtmp1, vd2_vjr, vd2_inp);
            vd2_inp = _mm_fmadd_pd(vd2_dtmp2, vd2_vji, vd2_inp);

            _mm_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
        }
    }
}

/* Apply a left ZLARF1F reflector where v(1) is implicit and equal to 1+0i.
   v uses f2c 1-based indexing (no pointer adjustment). The head block is
   synthesized in vector registers as [1+0i, v(2)]; ZLARF1F does not
   reference v(1), so its real/imag lanes are masked out of the load. */
void fla_zlarf1f_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff,
                                       aocl_int64_t ldr, dcomplex *v, dcomplex *ntau,
                                       dcomplex *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp1, vd2_dtmp2, vd2_vj, vd2_vjr, vd2_vji;
    __m256d vd4_inp, vd4_dtmp1, vd4_dtmp2, vd4_vj, vd4_vjr, vd4_vji, vd4_vhead;
    __m128d vd2_ltmp, vd2_htmp;
    __m128d vd2_one_neg_one = _mm_set_pd(-1.0, 1.0);

    arows = m;
    acols = n;
    vd2_ntau = _mm_loadu_pd((const doublereal *)ntau);

    /* Load v(2) from &v[1] with lanes 2-3 enabled; v(1) lanes are masked off.
       Blend in 1+0i so vd4_vhead = [Vr1, Vi1, Vr2, Vi2] = [1, 0, v(2).r, v(2).i]. */
    vd4_vhead = _mm256_blend_pd(
        _mm256_maskload_pd((const doublereal *)&v[1], _mm256_set_epi64x(-1, -1, 0, 0)),
        _mm256_set_pd(0.0, 0.0, 0.0, 1.0), 0x3);

    for(j = 1; j <= acols; j++)
    {
        vd2_dtmp1 = _mm_setzero_pd();
        vd2_dtmp2 = _mm_setzero_pd();
        vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[1 + j * ldr]);
        vd4_vjr = _mm256_permute_pd(vd4_vhead, 0b0000);
        vd4_vji = _mm256_permute_pd(vd4_vhead, 0b1111);
        vd4_dtmp1 = _mm256_mul_pd(vd4_inp, vd4_vjr);
        vd4_dtmp2 = _mm256_mul_pd(vd4_inp, vd4_vji);

        for(k = 3; k <= (arows - 1); k += 2)
        {
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            vd4_vjr = _mm256_permute_pd(vd4_vj, 0b0000);
            vd4_vji = _mm256_permute_pd(vd4_vj, 0b1111);

            vd4_dtmp1 = _mm256_fmadd_pd(vd4_inp, vd4_vjr, vd4_dtmp1);
            vd4_dtmp2 = _mm256_fmadd_pd(vd4_inp, vd4_vji, vd4_dtmp2);
        }
        if(k == arows)
        {
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj = _mm_loadu_pd((const doublereal *)&v[k]);

            vd2_vjr = _mm_permute_pd(vd2_vj, 0b00);
            vd2_vji = _mm_permute_pd(vd2_vj, 0b11);

            vd2_dtmp1 = _mm_fmadd_pd(vd2_inp, vd2_vjr, vd2_dtmp1);
            vd2_dtmp2 = _mm_fmadd_pd(vd2_inp, vd2_vji, vd2_dtmp2);
        }

        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp1);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp1, 0x1);
        vd2_dtmp1 = _mm_add_pd(vd2_dtmp1, vd2_ltmp);
        vd2_dtmp1 = _mm_add_pd(vd2_dtmp1, vd2_htmp);

        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp2);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp2, 0x1);
        vd2_dtmp2 = _mm_add_pd(vd2_dtmp2, vd2_ltmp);
        vd2_dtmp2 = _mm_add_pd(vd2_dtmp2, vd2_htmp);

        vd2_dtmp2 = _mm_permute_pd(vd2_dtmp2, 0b01);
        vd2_dtmp1 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);
        vd2_dtmp1 = _mm_add_pd(vd2_dtmp1, vd2_dtmp2);

        /* The blended head block already includes C(1,j) with v(1)=1, so do
           not add the first row again here. */
        _mm_storeu_pd((doublereal *)&work[j], vd2_dtmp1);

        vd2_dtmp1 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);

        vd2_vjr = _mm_permute_pd(vd2_ntau, 0b00);
        vd2_vji = _mm_permute_pd(vd2_ntau, 0b11);

        vd2_dtmp2 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);
        vd2_dtmp2 = _mm_permute_pd(vd2_dtmp2, 0b01);
        vd2_dtmp1 = _mm_mul_pd(vd2_dtmp1, vd2_vjr);
        vd2_dtmp1 = _mm_fmadd_pd(vd2_vji, vd2_dtmp2, vd2_dtmp1);

        vd2_dtmp2 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);
        vd2_dtmp2 = _mm_permute_pd(vd2_dtmp2, 0b01);

        vd4_dtmp1 = _mm256_castpd128_pd256(vd2_dtmp1);
        vd4_dtmp1 = _mm256_insertf128_pd(vd4_dtmp1, vd2_dtmp1, 0x1);

        vd4_dtmp2 = _mm256_castpd128_pd256(vd2_dtmp2);
        vd4_dtmp2 = _mm256_insertf128_pd(vd4_dtmp2, vd2_dtmp2, 0x1);

        vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[1 + j * ldr]);
        vd4_vjr = _mm256_permute_pd(vd4_vhead, 0b0000);
        vd4_vji = _mm256_permute_pd(vd4_vhead, 0b1111);
        vd4_inp = _mm256_fmadd_pd(vd4_dtmp1, vd4_vjr, vd4_inp);
        vd4_inp = _mm256_fmadd_pd(vd4_dtmp2, vd4_vji, vd4_inp);
        _mm256_storeu_pd((doublereal *)&a_buff[1 + j * ldr], vd4_inp);

        for(k = 3; k <= (arows - 1); k += 2)
        {
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            vd4_vjr = _mm256_permute_pd(vd4_vj, 0b0000);
            vd4_vji = _mm256_permute_pd(vd4_vj, 0b1111);

            vd4_inp = _mm256_fmadd_pd(vd4_dtmp1, vd4_vjr, vd4_inp);
            vd4_inp = _mm256_fmadd_pd(vd4_dtmp2, vd4_vji, vd4_inp);

            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
        }
        if(k == arows)
        {
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj = _mm_loadu_pd((const doublereal *)&v[k]);

            vd2_vjr = _mm_permute_pd(vd2_vj, 0b00);
            vd2_vji = _mm_permute_pd(vd2_vj, 0b11);

            vd2_inp = _mm_fmadd_pd(vd2_dtmp1, vd2_vjr, vd2_inp);
            vd2_inp = _mm_fmadd_pd(vd2_dtmp2, vd2_vji, vd2_inp);

            _mm_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
        }
    }
}

/* Scalar complex multiply used for the cleanup row in the right-apply path. */
static inline dcomplex fla_zlarf1f_mul_avx2(dcomplex a, dcomplex b)
{
    dcomplex c;
    c.real = a.real * b.real - a.imag * b.imag;
    c.imag = a.real * b.imag + a.imag * b.real;
    return c;
}

/* Scalar conjugate helper for the cleanup row in the right-apply path. */
static inline dcomplex fla_zlarf1f_conj_avx2(dcomplex a)
{
    a.imag = -a.imag;
    return a;
}

/* Multiply two packed complex values by a scalar complex value. The input
   vector layout is [r0, i0, r1, i1]. */
static inline __m256d fla_zlarf1f_mul2_avx2(__m256d a, doublereal b_real, doublereal b_imag)
{
    const __m256d sign = _mm256_set_pd(1.0, -1.0, 1.0, -1.0);
    __m256d real_part = _mm256_mul_pd(a, _mm256_set1_pd(b_real));
    __m256d imag_part
        = _mm256_mul_pd(_mm256_permute_pd(a, 0b0101), _mm256_mul_pd(_mm256_set1_pd(b_imag), sign));
    return _mm256_add_pd(real_part, imag_part);
}

/*
 * Apply a right ZLARF1F reflector where v(1) is implicit.
 * The column of A that pairs with the implicit v element 1 is not part of
 * a_buff, it is passed separately in c1. The remaining n columns start at
 * a_buff and pair with the n elements of v.
 *
 * @note: n is the number of columns of A excluding c1 (i.e. the number of columns of A minus 1)
 */
void fla_zlarf1_right_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff,
                                       aocl_int64_t ldr, dcomplex *v, dcomplex *c1, dcomplex *ntau,
                                       dcomplex *work)
{
    aocl_int64_t i, j;
    __m256d vd4_w, vd4_a, vd4_update;

    /* Process two complex rows at a time. c1 holds the matrix column paired
       with the implicit v(1)=1 and is folded into w before the rank-1 update. */
    for(i = 0; i <= m - 2; i += 2)
    {
        vd4_w = _mm256_loadu_pd((const doublereal *)&c1[i]);
        for(j = 0; j < n; ++j)
        {
            vd4_a = _mm256_loadu_pd((const doublereal *)&a_buff[i + j * ldr]);
            vd4_w = _mm256_add_pd(vd4_w, fla_zlarf1f_mul2_avx2(vd4_a, v[j].real, v[j].imag));
        }

        _mm256_storeu_pd((doublereal *)&work[i], vd4_w);
        vd4_update = fla_zlarf1f_mul2_avx2(vd4_w, ntau->real, ntau->imag);

        vd4_a = _mm256_loadu_pd((const doublereal *)&c1[i]);
        vd4_a = _mm256_add_pd(vd4_a, vd4_update);
        _mm256_storeu_pd((doublereal *)&c1[i], vd4_a);

        for(j = 0; j < n; ++j)
        {
            vd4_a = _mm256_loadu_pd((const doublereal *)&a_buff[i + j * ldr]);
            vd4_a = _mm256_add_pd(vd4_a, fla_zlarf1f_mul2_avx2(vd4_update, v[j].real, -v[j].imag));
            _mm256_storeu_pd((doublereal *)&a_buff[i + j * ldr], vd4_a);
        }
    }

    for(; i < m; ++i)
    {
        dcomplex w = c1[i];
        dcomplex update;

        for(j = 0; j < n; ++j)
        {
            dcomplex prod = fla_zlarf1f_mul_avx2(a_buff[i + j * ldr], v[j]);
            w.real += prod.real;
            w.imag += prod.imag;
        }

        work[i] = w;
        update = fla_zlarf1f_mul_avx2(*ntau, w);
        c1[i].real += update.real;
        c1[i].imag += update.imag;

        for(j = 0; j < n; ++j)
        {
            dcomplex v_conj = fla_zlarf1f_conj_avx2(v[j]);
            dcomplex scaled = fla_zlarf1f_mul_avx2(update, v_conj);
            a_buff[i + j * ldr].real += scaled.real;
            a_buff[i + j * ldr].imag += scaled.imag;
        }
    }
}

#endif