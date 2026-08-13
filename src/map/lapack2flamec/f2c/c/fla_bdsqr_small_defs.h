/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_BDSQR_SMALL_DEFS_H
#define FLA_BDSQR_SMALL_DEFS_H

/*! @file fla_bdsqr_small_defs.h
 *  @brief Precision bindings for the small bidiagonal QR iteration.
 *  */

#if FLA_ENABLE_AMD_OPT

#define FLA_BDSQR_SMALL_CAT_(a, b) a##b
#define FLA_BDSQR_SMALL_CAT(a, b) FLA_BDSQR_SMALL_CAT_(a, b)

#define d_BDSQR_FNAME lapack_dbdsqr_small
#define d_BDSQR_ELT doublereal
#define d_BDSQR_LAMCH dlamch_
#define d_BDSQR_LARTG dlartg_
#define d_BDSQR_LAS2 dlas2_
#define d_BDSQR_LASV2 dlasv2_
#define d_BDSQR_SIGN d_sign
#define d_BDSQR_ROT aocl_blas_drot
#define d_BDSQR_SCAL aocl_blas_dscal
#define d_BDSQR_SWAP aocl_blas_dswap

#define s_BDSQR_FNAME lapack_sbdsqr_small
#define s_BDSQR_ELT real
#define s_BDSQR_LAMCH slamch_
#define s_BDSQR_LARTG slartg_
#define s_BDSQR_LAS2 slas2_
#define s_BDSQR_LASV2 slasv2_
#define s_BDSQR_SIGN r_sign
#define s_BDSQR_ROT aocl_blas_srot
#define s_BDSQR_SCAL aocl_blas_sscal
#define s_BDSQR_SWAP aocl_blas_sswap

#define d_BDSQR_DECL_EXTRA
#define s_BDSQR_DECL_EXTRA doublereal dpow;

#define d_BDSQR_TOLMUL           \
    d__3 = 100.;                 \
    d__4 = pow_dd(&eps, &c_b15); \
    d__1 = 10.;                  \
    d__2 = fla_min(d__3, d__4);  \
    tolmul = fla_max(d__1, d__2);

#define s_BDSQR_TOLMUL            \
    dpow = (doublereal)eps;       \
    d__3 = 100.f;                 \
    d__4 = pow_dd(&dpow, &c_b15); \
    d__1 = 10.f;                  \
    d__2 = fla_min(d__3, d__4);   \
    tolmul = fla_max(d__1, d__2);

/* Application of Givens Rotation ** T)
 * over rows row & row + 1
 * from the left */
#define FLA_APPLY_GIVENS_GLVX(P, idim, imat, ldi, row, cs, sn) \
    {                                                          \
        aocl_int64_t im;                                       \
        FLA_BDSQR_SMALL_CAT(P, _BDSQR_ELT) tv0, tv1;           \
        for(im = 1; im <= *idim; im++)                         \
        {                                                      \
            tv0 = imat[row + 0 + im * *ldi];                   \
            tv1 = imat[row + 1 + im * *ldi];                   \
                                                               \
            imat[row + 0 + im * *ldi] = cs * tv0 + sn * tv1;   \
            imat[row + 1 + im * *ldi] = cs * tv1 - sn * tv0;   \
        }                                                      \
    }
/* Application of Givens Rotation ** T)
 * over columns col & col + 1
 * from the right */
#define FLA_APPLY_GIVENS_GRVX(P, idim, imat, ldi, col, cs, sn) \
    {                                                          \
        aocl_int64_t im;                                       \
        FLA_BDSQR_SMALL_CAT(P, _BDSQR_ELT) tv0, tv1;           \
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
#endif /* FLA_BDSQR_SMALL_DEFS_H */
