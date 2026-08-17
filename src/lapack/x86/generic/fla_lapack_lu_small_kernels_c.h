/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 ******************************************************************************/

/*! @file fla_lapack_lu_small_kernels_c.h
 *  @brief Fixed-size LU micro-kernels for single-complex GETRF (macro-only).
 */

#include "FLAME.h"

#if FLA_ENABLE_AMD_OPT

#define FLA_LU_C_ABS1_SC(z) (f2c_abs((z).real) + f2c_abs((z).imag))

/* Reciprocal 1/denom, used to scale a column by multiplication as
 * bl1_cinvscalv does. |denom|^2 is formed directly whenever the result is a
 * normal float, which is the case for any denominator between roughly 1e-19 and
 * 1e19 and costs a single division. Outside that band the square would overflow
 * to infinity or underflow to zero, so the guard falls back to Smith's scaled
 * form, matching bl1_cinvert2s exactly. Callers must ensure denom != 0. */
#define FLA_LU_C_RECIP_SC(rptr, denom)                     \
    do                                                     \
    {                                                      \
        real _dr = (denom).real;                           \
        real _di = (denom).imag;                           \
        real _den = _dr * _dr + _di * _di;                 \
        if(_den >= FLT_MIN && _den <= FLT_MAX)             \
        {                                                  \
            real _inv = 1.0f / _den;                       \
            (rptr)->real = _dr * _inv;                     \
            (rptr)->imag = -_di * _inv;                    \
        }                                                  \
        else                                               \
        {                                                  \
            real _s = fla_max(f2c_abs(_dr), f2c_abs(_di)); \
            real _dr_s = _dr / _s;                         \
            real _di_s = _di / _s;                         \
            real _d2 = _dr_s * _dr + _di_s * _di;          \
            (rptr)->real = _dr_s / _d2;                    \
            (rptr)->imag = -_di_s / _d2;                   \
        }                                                  \
    } while(0)

#define FLA_LU_C_MUL_SC(yptr, alpha)          \
    do                                        \
    {                                         \
        real _ar = (alpha).real;              \
        real _ai = (alpha).imag;              \
        real _yr = (yptr)->real;              \
        real _yi = (yptr)->imag;              \
        (yptr)->real = _yr * _ar - _yi * _ai; \
        (yptr)->imag = _yr * _ai + _yi * _ar; \
    } while(0)

#define FLA_LU_C_SUBMUL_SC(y, alpha, x)    \
    do                                     \
    {                                      \
        real _ar = (alpha).real;           \
        real _ai = (alpha).imag;           \
        real _xr = (x).real;               \
        real _xi = (x).imag;               \
        (y).real -= _ar * _xr - _ai * _xi; \
        (y).imag -= _ar * _xi + _ai * _xr; \
    } while(0)

#define FLA_LU_C_SWAP_SC(a, b, t) \
    do                            \
    {                             \
        (t) = (a);                \
        (a) = (b);                \
        (b) = (t);                \
    } while(0)

/*
 * LU 1x1 with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_C_1x1(i, n, buff_A, ldim_A, buff_p, info)                           \
    buff_p[i] = (aocl_int_t)((i) + 1);                                                       \
    if((buff_A[i + (*ldim_A) * i].real == 0.0f) && (buff_A[i + (*ldim_A) * i].imag == 0.0f)) \
        info = (info == 0) ? ((i) + 1) : info;

/*
 * LU 2x2 with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_GEN_C_2x2(i, n, buff_A, ldim_A, buff_p, info)            \
    real max_val_2 = 0.0f;                                                        \
    scomplex t_sc_2, recip_2;                                                     \
    scomplex *acur_2, *apiv_2, *asrc_2;                                           \
    aocl_int64_t i_2, p_idx_2 = i, lda2 = *ldim_A;                                \
    acur_2 = &buff_A[i + lda2 * i];                                               \
    for(i_2 = 0; i_2 < 2; i_2++)                                                  \
    {                                                                             \
        real t_abs_2 = FLA_LU_C_ABS1_SC(acur_2[i_2]);                             \
        if(t_abs_2 > max_val_2)                                                   \
        {                                                                         \
            max_val_2 = t_abs_2;                                                  \
            p_idx_2 = i + i_2;                                                    \
        }                                                                         \
    }                                                                             \
    apiv_2 = buff_A + p_idx_2;                                                    \
    asrc_2 = buff_A + i;                                                          \
    buff_p[i] = (aocl_int_t)(p_idx_2 + 1);                                        \
    if(max_val_2 != 0.0f)                                                         \
    {                                                                             \
        if(p_idx_2 != i)                                                          \
        {                                                                         \
            FLA_LU_C_SWAP_SC(apiv_2[0], asrc_2[0], t_sc_2);                       \
            FLA_LU_C_SWAP_SC(apiv_2[lda2], asrc_2[lda2], t_sc_2);                 \
            if(n >= 3)                                                            \
            {                                                                     \
                FLA_LU_C_SWAP_SC(apiv_2[2 * lda2], asrc_2[2 * lda2], t_sc_2);     \
                if(n == 4)                                                        \
                    FLA_LU_C_SWAP_SC(apiv_2[3 * lda2], asrc_2[3 * lda2], t_sc_2); \
            }                                                                     \
        }                                                                         \
        FLA_LU_C_RECIP_SC(&recip_2, *acur_2);                                     \
        FLA_LU_C_MUL_SC(&acur_2[1], recip_2);                                     \
        FLA_LU_C_SUBMUL_SC(acur_2[1 + lda2], acur_2[1], acur_2[lda2]);            \
    }                                                                             \
    else                                                                          \
    {                                                                             \
        info = (info == 0) ? p_idx_2 + 1 : info;                                  \
    }                                                                             \
    i = i + 1;                                                                    \
    acur_2 = &buff_A[i + lda2 * i];                                               \
    p_idx_2 = i;                                                                  \
    max_val_2 = FLA_LU_C_ABS1_SC(acur_2[0]);                                      \
    apiv_2 = buff_A + p_idx_2;                                                    \
    asrc_2 = buff_A + i;                                                          \
    buff_p[i] = (aocl_int_t)(p_idx_2 + 1);                                        \
    if(max_val_2 != 0.0f)                                                         \
    {                                                                             \
        if(p_idx_2 != i)                                                          \
        {                                                                         \
            FLA_LU_C_SWAP_SC(apiv_2[0], asrc_2[0], t_sc_2);                       \
            FLA_LU_C_SWAP_SC(apiv_2[lda2], asrc_2[lda2], t_sc_2);                 \
            if(n >= 3)                                                            \
            {                                                                     \
                FLA_LU_C_SWAP_SC(apiv_2[2 * lda2], asrc_2[2 * lda2], t_sc_2);     \
                if(n == 4)                                                        \
                    FLA_LU_C_SWAP_SC(apiv_2[3 * lda2], asrc_2[3 * lda2], t_sc_2); \
            }                                                                     \
        }                                                                         \
    }                                                                             \
    else                                                                          \
    {                                                                             \
        info = (info == 0) ? p_idx_2 + 1 : info;                                  \
    }

/*
 * LU 3x3 with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_GEN_C_3x3(i, n, buff_A, ldim_A, buff_p, info)         \
    aocl_int64_t i_3;                                                          \
    real max_val_3 = 0.0f;                                                     \
    scomplex t_sc_3, recip_3;                                                  \
    scomplex *acur_3, *apiv_3, *asrc_3;                                        \
    aocl_int64_t p_idx_3 = i, lda3 = *ldim_A;                                  \
    acur_3 = &buff_A[i + lda3 * i];                                            \
    for(i_3 = 0; i_3 < 3; i_3++)                                               \
    {                                                                          \
        real t_abs_3 = FLA_LU_C_ABS1_SC(acur_3[i_3]);                          \
        if(t_abs_3 > max_val_3)                                                \
        {                                                                      \
            max_val_3 = t_abs_3;                                               \
            p_idx_3 = i + i_3;                                                 \
        }                                                                      \
    }                                                                          \
    apiv_3 = buff_A + p_idx_3;                                                 \
    asrc_3 = buff_A + i;                                                       \
    buff_p[i] = (aocl_int_t)(p_idx_3 + 1);                                     \
    if(max_val_3 != 0.0f)                                                      \
    {                                                                          \
        if(p_idx_3 != i)                                                       \
        {                                                                      \
            FLA_LU_C_SWAP_SC(apiv_3[0], asrc_3[0], t_sc_3);                    \
            FLA_LU_C_SWAP_SC(apiv_3[lda3], asrc_3[lda3], t_sc_3);              \
            FLA_LU_C_SWAP_SC(apiv_3[2 * lda3], asrc_3[2 * lda3], t_sc_3);      \
            if(n == 4)                                                         \
                FLA_LU_C_SWAP_SC(apiv_3[3 * lda3], asrc_3[3 * lda3], t_sc_3);  \
        }                                                                      \
        FLA_LU_C_RECIP_SC(&recip_3, *acur_3);                                  \
        FLA_LU_C_MUL_SC(&acur_3[1], recip_3);                                  \
        FLA_LU_C_SUBMUL_SC(acur_3[1 + lda3], acur_3[1], acur_3[lda3]);         \
        FLA_LU_C_SUBMUL_SC(acur_3[1 + 2 * lda3], acur_3[1], acur_3[2 * lda3]); \
        FLA_LU_C_MUL_SC(&acur_3[2], recip_3);                                  \
        FLA_LU_C_SUBMUL_SC(acur_3[2 + lda3], acur_3[2], acur_3[lda3]);         \
        FLA_LU_C_SUBMUL_SC(acur_3[2 + 2 * lda3], acur_3[2], acur_3[2 * lda3]); \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        info = (info == 0) ? p_idx_3 + 1 : info;                               \
    }                                                                          \
    i = i + 1;                                                                 \
    FLA_LU_PIV_SMALL_GEN_C_2x2(i, n, buff_A, ldim_A, buff_p, info);

/*
 * LU 4x4 with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_GEN_C_4x4(i, n, buff_A, ldim_A, buff_p, info) \
    aocl_int64_t i_1;                                                  \
    real max_val = 0.0f;                                               \
    scomplex t_sc, recip;                                              \
    scomplex *acur, *apiv, *asrc;                                      \
    aocl_int64_t p_idx = i, lda = *ldim_A;                             \
    acur = &buff_A[i + lda * i];                                       \
    for(i_1 = 0; i_1 < 4; i_1++)                                       \
    {                                                                  \
        real t_abs = FLA_LU_C_ABS1_SC(acur[i_1]);                      \
        if(t_abs > max_val)                                            \
        {                                                              \
            max_val = t_abs;                                           \
            p_idx = i + i_1;                                           \
        }                                                              \
    }                                                                  \
    apiv = buff_A + p_idx;                                             \
    asrc = buff_A + i;                                                 \
    buff_p[i] = (aocl_int_t)(p_idx + 1);                               \
    if(max_val != 0.0f)                                                \
    {                                                                  \
        if(p_idx != i)                                                 \
        {                                                              \
            FLA_LU_C_SWAP_SC(apiv[0], asrc[0], t_sc);                  \
            FLA_LU_C_SWAP_SC(apiv[lda], asrc[lda], t_sc);              \
            FLA_LU_C_SWAP_SC(apiv[2 * lda], asrc[2 * lda], t_sc);      \
            FLA_LU_C_SWAP_SC(apiv[3 * lda], asrc[3 * lda], t_sc);      \
        }                                                              \
        FLA_LU_C_RECIP_SC(&recip, *acur);                              \
        FLA_LU_C_MUL_SC(&acur[1], recip);                              \
        FLA_LU_C_SUBMUL_SC(acur[1 + lda], acur[1], acur[lda]);         \
        FLA_LU_C_SUBMUL_SC(acur[1 + 2 * lda], acur[1], acur[2 * lda]); \
        FLA_LU_C_SUBMUL_SC(acur[1 + 3 * lda], acur[1], acur[3 * lda]); \
        FLA_LU_C_MUL_SC(&acur[2], recip);                              \
        FLA_LU_C_SUBMUL_SC(acur[2 + lda], acur[2], acur[lda]);         \
        FLA_LU_C_SUBMUL_SC(acur[2 + 2 * lda], acur[2], acur[2 * lda]); \
        FLA_LU_C_SUBMUL_SC(acur[2 + 3 * lda], acur[2], acur[3 * lda]); \
        FLA_LU_C_MUL_SC(&acur[3], recip);                              \
        FLA_LU_C_SUBMUL_SC(acur[3 + lda], acur[3], acur[lda]);         \
        FLA_LU_C_SUBMUL_SC(acur[3 + 2 * lda], acur[3], acur[2 * lda]); \
        FLA_LU_C_SUBMUL_SC(acur[3 + 3 * lda], acur[3], acur[3 * lda]); \
    }                                                                  \
    else                                                               \
    {                                                                  \
        info = (info == 0) ? p_idx + 1 : info;                         \
    }                                                                  \
    i = i + 1;                                                         \
    FLA_LU_PIV_SMALL_GEN_C_3x3(i, n, buff_A, ldim_A, buff_p, info);

#define FLA_LU_PIV_SMALL_C_2x2(i, n, buff_A, ldim_A, buff_p, info) \
    FLA_LU_PIV_SMALL_GEN_C_2x2(i, n, buff_A, ldim_A, buff_p, info)

#define FLA_LU_PIV_SMALL_C_3x3(i, n, buff_A, ldim_A, buff_p, info) \
    FLA_LU_PIV_SMALL_GEN_C_3x3(i, n, buff_A, ldim_A, buff_p, info)

#define FLA_LU_PIV_SMALL_C_4x4(i, n, buff_A, ldim_A, buff_p, info) \
    FLA_LU_PIV_SMALL_GEN_C_4x4(i, n, buff_A, ldim_A, buff_p, info)

#endif
