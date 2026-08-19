/******************************************************************************
 * * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* Apply a left Householder reflector when the full vector v is stored with
   unit stride. This fuses the column dot products and rank-1 update for small
   complex-double panels using AVX512/FMA operations. */
void fla_zlarf_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff,
                                       aocl_int64_t ldr, dcomplex *v, dcomplex *ntau,
                                       dcomplex *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp1, vd2_dtmp2, vd2_vj, vd2_vjr, vd2_vji;
    __m128d vd2_ltmp, vd2_htmp;
    __m256d vd4_inp, vd4_dtmp1, vd4_dtmp2, vd4_vj, vd4_vjr, vd4_vji;
    __m512d vd8_dtmp1, vd8_dtmp2, vd8_inp, vd8_vj, vd8_vjr, vd8_vji;
    __m256d vd4_ltmp, vd4_htmp;
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
        vd8_dtmp1 = _mm512_setzero_pd();
        vd8_dtmp2 = _mm512_setzero_pd();

        /* Compute tmp = c_A**T . v */
        for(k = 1; k <= (arows - 3); k += 4)
        {
            /* load column elements of A and v */
            /* Column is loaded in the following format */
            /* [Ar0, Ai0, Ar1, Ai1, Ar2, Ai2, Ar3, Ai3] */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            /* Vector v is loaded in following format */
            /* [Vr0, Vi0, Vr1, Vi1, Vr2, Vi2, Vr3, Vi3] */
            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* Rearrange vd8 as follows */
            /* vd8_vjr = [Vr0, Vr0, Vr1, Vr1, Vr2, Vr2, Vr3, Vr3] */
            vd8_vjr = _mm512_permute_pd(vd8_vj, 0b00000000);
            /* vd8_vji = [Vi0, Vi0, Vi1, Vi1, Vi2, Vi2, Vi3, Vi3] */
            vd8_vji = _mm512_permute_pd(vd8_vj, 0b11111111);

            /* take dot product */
            vd8_dtmp1 = _mm512_fmadd_pd(vd8_inp, vd8_vjr, vd8_dtmp1);
            vd8_dtmp2 = _mm512_fmadd_pd(vd8_inp, vd8_vji, vd8_dtmp2);
        }
        if(k <= (arows - 1))
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

            k += 2;
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

        /* Etract Upper and lower 256 bits of vd8_dtmp1 */
        vd4_ltmp = _mm512_castpd512_pd256(vd8_dtmp1);
        vd4_htmp = _mm512_extractf64x4_pd(vd8_dtmp1, 0x1);

        /* Add the lower and upper 256 bits with vd4_dtmp1 */
        vd4_dtmp1 = _mm256_add_pd(vd4_dtmp1, vd4_ltmp);
        vd4_dtmp1 = _mm256_add_pd(vd4_dtmp1, vd4_htmp);

        /* Etract Upper and lower 256 bits of vd8_dtmp2 */
        vd4_ltmp = _mm512_castpd512_pd256(vd8_dtmp2);
        vd4_htmp = _mm512_extractf64x4_pd(vd8_dtmp2, 0x1);

        /* Add the lower and upper 256 bits with vd4_dtmp2 */
        vd4_dtmp2 = _mm256_add_pd(vd4_dtmp2, vd4_ltmp);
        vd4_dtmp2 = _mm256_add_pd(vd4_dtmp2, vd4_htmp);

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

        /* Compute tmp = ntau * tmp */
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
        /* vd2_dtmp1 = [Kr * Tr - Ki * Ti, Kr * Ti + Ki * Tr] */
        vd2_dtmp1 = _mm_fmadd_pd(vd2_vji, vd2_dtmp2, vd2_dtmp1);

        /* Set the first element as negative [ Pr -Pi ] */
        vd2_dtmp2 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);
        /* Shuffle vd2_dtmp2 = [ -Pi Pr ] */
        vd2_dtmp2 = _mm_permute_pd(vd2_dtmp2, 0b01);

        /* Broadcast value to __m512d and _m256d */
        vd4_dtmp1 = _mm256_castpd128_pd256(vd2_dtmp1);
        vd4_dtmp1 = _mm256_insertf128_pd(vd4_dtmp1, vd2_dtmp1, 0x1);

        vd4_dtmp2 = _mm256_castpd128_pd256(vd2_dtmp2);
        vd4_dtmp2 = _mm256_insertf128_pd(vd4_dtmp2, vd2_dtmp2, 0x1);

        vd8_dtmp1 = _mm512_castpd256_pd512(vd4_dtmp1);
        vd8_dtmp1 = _mm512_insertf64x4(vd8_dtmp1, vd4_dtmp1, 0x1);

        vd8_dtmp2 = _mm512_castpd256_pd512(vd4_dtmp2);
        vd8_dtmp2 = _mm512_insertf64x4(vd8_dtmp2, vd4_dtmp2, 0x1);

        /* Compute c_A + (tmp * v')' = c_A + (v * tmp')*/
        for(k = 1; k <= (arows - 3); k += 4)
        {
            /* load column elements of c_A and v */
            /* vd8_inp = [ Ar0, Ai0, Ar1, Ai1, Ar2, Ai2, Ar3, Ai3  ] */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            /* vd8_vj =  [ Vr0, Vi0, Vr1, Vi1, Vr2, Vi2, Vr3, Vi3 ] */
            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* vd8_vjr = [ Vr0, Vr0, Vr1, Vr1, Vr2, Vr2, Vr3, Vr3 ] */
            vd8_vjr = _mm512_permute_pd(vd8_vj, 0b00000000);
            /* vd8_vji = [ Vi0, Vi0, Vi1, Vi1, Vi2, Vi2, Vi3, Vi3 ] */
            vd8_vji = _mm512_permute_pd(vd8_vj, 0b11111111);

            /* mul by dtmp, add and store */
            /* inp + [  Pr * Vr, Pi * Vr ] */
            vd8_inp = _mm512_fmadd_pd(vd8_dtmp1, vd8_vjr, vd8_inp);
            /* inp + [ Pr * Vr - Pi * Vi, Pi * Vr + Pr * Vi ] */
            vd8_inp = _mm512_fmadd_pd(vd8_dtmp2, vd8_vji, vd8_inp);

            _mm512_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd8_inp);
        }
        if(k <= (arows - 1))
        {
            /* Same steps followed as in above loop */
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            vd4_vjr = _mm256_permute_pd(vd4_vj, 0b0000);
            vd4_vji = _mm256_permute_pd(vd4_vj, 0b1111);

            vd4_inp = _mm256_fmadd_pd(vd4_dtmp1, vd4_vjr, vd4_inp);
            vd4_inp = _mm256_fmadd_pd(vd4_dtmp2, vd4_vji, vd4_inp);

            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
            k += 2;
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
   synthesized in vector registers as [1+0i, v(2), v(3), v(4)]; ZLARF1F does
   not reference v(1), so its real/imag lanes are masked out of the load. */
void fla_zlarf1f_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff,
                                         aocl_int64_t ldr, dcomplex *v, dcomplex *ntau,
                                         dcomplex *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp1, vd2_dtmp2, vd2_vj, vd2_vjr, vd2_vji;
    __m128d vd2_ltmp, vd2_htmp;
    __m256d vd4_inp, vd4_dtmp1, vd4_dtmp2, vd4_vj, vd4_vjr, vd4_vji;
    __m512d vd8_dtmp1, vd8_dtmp2, vd8_inp, vd8_vj, vd8_vjr, vd8_vji, vd8_vhead;
    __m256d vd4_ltmp, vd4_htmp;
    __m128d vd2_one_neg_one = _mm_set_pd(-1.0, 1.0);

    arows = m;
    acols = n;

    vd2_ntau = _mm_loadu_pd((const doublereal *)ntau);

    /* Load v(2:4) from &v[1] with mask 0xfc; v(1) lanes are masked off.
       Blend in 1+0i so vd8_vhead = [1, 0, v(2), v(3), v(4)] as packed doubles. */
    vd8_vhead = _mm512_mask_loadu_pd(_mm512_set_pd(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0), 0xfc,
                                     (const doublereal *)&v[1]);

    for(j = 1; j <= acols; j++)
    {
        vd2_dtmp1 = _mm_setzero_pd();
        vd2_dtmp2 = _mm_setzero_pd();
        vd4_dtmp1 = _mm256_setzero_pd();
        vd4_dtmp2 = _mm256_setzero_pd();
        vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[1 + j * ldr]);
        vd8_vjr = _mm512_permute_pd(vd8_vhead, 0b00000000);
        vd8_vji = _mm512_permute_pd(vd8_vhead, 0b11111111);
        vd8_dtmp1 = _mm512_mul_pd(vd8_inp, vd8_vjr);
        vd8_dtmp2 = _mm512_mul_pd(vd8_inp, vd8_vji);

        for(k = 5; k <= (arows - 3); k += 4)
        {
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            vd8_vjr = _mm512_permute_pd(vd8_vj, 0b00000000);
            vd8_vji = _mm512_permute_pd(vd8_vj, 0b11111111);

            vd8_dtmp1 = _mm512_fmadd_pd(vd8_inp, vd8_vjr, vd8_dtmp1);
            vd8_dtmp2 = _mm512_fmadd_pd(vd8_inp, vd8_vji, vd8_dtmp2);
        }
        if(k <= (arows - 1))
        {
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            vd4_vjr = _mm256_permute_pd(vd4_vj, 0b0000);
            vd4_vji = _mm256_permute_pd(vd4_vj, 0b1111);

            vd4_dtmp1 = _mm256_fmadd_pd(vd4_inp, vd4_vjr, vd4_dtmp1);
            vd4_dtmp2 = _mm256_fmadd_pd(vd4_inp, vd4_vji, vd4_dtmp2);

            k += 2;
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

        vd4_ltmp = _mm512_castpd512_pd256(vd8_dtmp1);
        vd4_htmp = _mm512_extractf64x4_pd(vd8_dtmp1, 0x1);
        vd4_dtmp1 = _mm256_add_pd(vd4_dtmp1, vd4_ltmp);
        vd4_dtmp1 = _mm256_add_pd(vd4_dtmp1, vd4_htmp);

        vd4_ltmp = _mm512_castpd512_pd256(vd8_dtmp2);
        vd4_htmp = _mm512_extractf64x4_pd(vd8_dtmp2, 0x1);
        vd4_dtmp2 = _mm256_add_pd(vd4_dtmp2, vd4_ltmp);
        vd4_dtmp2 = _mm256_add_pd(vd4_dtmp2, vd4_htmp);

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

        /* The masked head block already includes C(1,j) with v(1)=1, so do
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

        vd8_dtmp1 = _mm512_castpd256_pd512(vd4_dtmp1);
        vd8_dtmp1 = _mm512_insertf64x4(vd8_dtmp1, vd4_dtmp1, 0x1);

        vd8_dtmp2 = _mm512_castpd256_pd512(vd4_dtmp2);
        vd8_dtmp2 = _mm512_insertf64x4(vd8_dtmp2, vd4_dtmp2, 0x1);

        vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[1 + j * ldr]);
        vd8_vjr = _mm512_permute_pd(vd8_vhead, 0b00000000);
        vd8_vji = _mm512_permute_pd(vd8_vhead, 0b11111111);
        vd8_inp = _mm512_fmadd_pd(vd8_dtmp1, vd8_vjr, vd8_inp);
        vd8_inp = _mm512_fmadd_pd(vd8_dtmp2, vd8_vji, vd8_inp);
        _mm512_storeu_pd((doublereal *)&a_buff[1 + j * ldr], vd8_inp);

        for(k = 5; k <= (arows - 3); k += 4)
        {
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            vd8_vjr = _mm512_permute_pd(vd8_vj, 0b00000000);
            vd8_vji = _mm512_permute_pd(vd8_vj, 0b11111111);

            vd8_inp = _mm512_fmadd_pd(vd8_dtmp1, vd8_vjr, vd8_inp);
            vd8_inp = _mm512_fmadd_pd(vd8_dtmp2, vd8_vji, vd8_inp);

            _mm512_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd8_inp);
        }
        if(k <= (arows - 1))
        {
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            vd4_vjr = _mm256_permute_pd(vd4_vj, 0b0000);
            vd4_vji = _mm256_permute_pd(vd4_vj, 0b1111);

            vd4_inp = _mm256_fmadd_pd(vd4_dtmp1, vd4_vjr, vd4_inp);
            vd4_inp = _mm256_fmadd_pd(vd4_dtmp2, vd4_vji, vd4_inp);

            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
            k += 2;
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
static inline dcomplex fla_zlarf1f_mul_avx512(dcomplex a, dcomplex b)
{
    dcomplex c;
    c.real = a.real * b.real - a.imag * b.imag;
    c.imag = a.real * b.imag + a.imag * b.real;
    return c;
}

/* Scalar conjugate helper for the cleanup row in the right-apply path. */
static inline dcomplex fla_zlarf1f_conj_avx512(dcomplex a)
{
    a.imag = -a.imag;
    return a;
}

/* Multiply four packed complex values by a scalar complex value. The input
   vector layout is [r0, i0, r1, i1, r2, i2, r3, i3]. */
static inline __m512d fla_zlarf1f_mul4_avx512(__m512d a, doublereal b_real, doublereal b_imag)
{
    const __m512d sign = _mm512_set_pd(1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0);
    __m512d real_part = _mm512_mul_pd(a, _mm512_set1_pd(b_real));
    __m512d imag_part = _mm512_mul_pd(_mm512_permute_pd(a, 0b01010101),
                                      _mm512_mul_pd(_mm512_set1_pd(b_imag), sign));
    return _mm512_add_pd(real_part, imag_part);
}

/*
 * Apply a right ZLARF1F reflector where v(1) is implicit.
 * The column of A that pairs with the implicit v element 1 is not part of
 * a_buff, it is passed separately in c1. The remaining n columns start at
 * a_buff and pair with the n elements of v.
 *
 * @note: n is the number of columns of A excluding c1 (i.e. the number of columns of A minus 1)
 */
void fla_zlarf1_right_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff,
                                         aocl_int64_t ldr, dcomplex *v, dcomplex *c1,
                                         dcomplex *ntau, dcomplex *work)
{
    aocl_int64_t i, j;
    __m512d vd8_w, vd8_a, vd8_update;

    /* Process four complex rows at a time. c1 is the column paired with the
       implicit v(1)=1 and is updated together with the explicit columns. */
    for(i = 0; i <= m - 4; i += 4)
    {
        vd8_w = _mm512_loadu_pd((const doublereal *)&c1[i]);
        for(j = 0; j < n; ++j)
        {
            vd8_a = _mm512_loadu_pd((const doublereal *)&a_buff[i + j * ldr]);
            vd8_w = _mm512_add_pd(vd8_w, fla_zlarf1f_mul4_avx512(vd8_a, v[j].real, v[j].imag));
        }

        _mm512_storeu_pd((doublereal *)&work[i], vd8_w);
        vd8_update = fla_zlarf1f_mul4_avx512(vd8_w, ntau->real, ntau->imag);

        vd8_a = _mm512_loadu_pd((const doublereal *)&c1[i]);
        vd8_a = _mm512_add_pd(vd8_a, vd8_update);
        _mm512_storeu_pd((doublereal *)&c1[i], vd8_a);

        for(j = 0; j < n; ++j)
        {
            vd8_a = _mm512_loadu_pd((const doublereal *)&a_buff[i + j * ldr]);
            vd8_a
                = _mm512_add_pd(vd8_a, fla_zlarf1f_mul4_avx512(vd8_update, v[j].real, -v[j].imag));
            _mm512_storeu_pd((doublereal *)&a_buff[i + j * ldr], vd8_a);
        }
    }

    for(; i < m; ++i)
    {
        dcomplex w = c1[i];
        dcomplex update;

        for(j = 0; j < n; ++j)
        {
            dcomplex prod = fla_zlarf1f_mul_avx512(a_buff[i + j * ldr], v[j]);
            w.real += prod.real;
            w.imag += prod.imag;
        }

        work[i] = w;
        update = fla_zlarf1f_mul_avx512(*ntau, w);
        c1[i].real += update.real;
        c1[i].imag += update.imag;

        for(j = 0; j < n; ++j)
        {
            dcomplex v_conj = fla_zlarf1f_conj_avx512(v[j]);
            dcomplex scaled = fla_zlarf1f_mul_avx512(update, v_conj);
            a_buff[i + j * ldr].real += scaled.real;
            a_buff[i + j * ldr].imag += scaled.imag;
        }
    }
}

#endif