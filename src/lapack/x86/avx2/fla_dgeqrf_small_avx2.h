/******************************************************************************
 * Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_DGEQRF_SMALL_AVX2_DEFS_H
#define FLA_DGEQRF_SMALL_AVX2_DEFS_H

/*! @file fla_dgeqrf_small_avx2.h
 *  @brief QR Kernels for small sizes.
 *  */

#if FLA_ENABLE_AMD_OPT

#include "fla_geqrf_small_avx2_kernel.h"

/* Application of Givens Rotation ** T)
 * over rows row & row + 1
 * from the left */
#define FLA_APPLY_GIVENS_LVX(idim, imat, ldi, row, cs, sn)   \
    {                                                        \
        aocl_int64_t im;                                     \
        doublereal tv0, tv1;                                 \
        for(im = 1; im <= *idim; im++)                       \
        {                                                    \
            tv0 = imat[row + 0 + im * *ldi];                 \
            tv1 = imat[row + 1 + im * *ldi];                 \
                                                             \
            imat[row + 0 + im * *ldi] = cs * tv0 + sn * tv1; \
            imat[row + 1 + im * *ldi] = cs * tv1 - sn * tv0; \
        }                                                    \
    }
/* Application of Givens Rotation ** T)
 * over columns col & col + 1
 * from the right */
#define FLA_APPLY_GIVENS_RVX(idim, imat, ldi, col, cs, sn)     \
    {                                                          \
        aocl_int64_t im;                                       \
        doublereal tv0, tv1;                                   \
        for(im = 1; im <= *idim; im++)                         \
        {                                                      \
            tv0 = imat[im + (col + 0) * *ldi];                 \
            tv1 = imat[im + (col + 1) * *ldi];                 \
                                                               \
            imat[im + (col + 0) * *ldi] = cs * tv0 + sn * tv1; \
            imat[im + (col + 1) * *ldi] = cs * tv1 - sn * tv0; \
        }                                                      \
    }

#endif /* FLA_ENABLE_AMD_OPT */
#endif /* FLA_DGEQRF_SMALL_AVX2_DEFS_H */
