/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_GESVD_DISPATCH_H
#define FLA_GESVD_DISPATCH_H

#include "FLAME.h"

/*! @file fla_gesvd_dispatch.h
 *  @brief Predicates selecting the small-size GESVD paths.
 *  */

#if FLA_ENABLE_AMD_OPT

/* Path 1: M much larger than N, JOBU='N' */
#define FLA_GESVD_SMALL_PATH1(WNTVO, M) \
    ((!(WNTVO)) && (M) <= FLA_GESVD_SMALL_SIZE_THRESH1 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))

/* Path 6: M much larger than N, JOBU='S', JOBVT='S' or 'A' */
#define FLA_GESVD_SMALL_PATH6(M) \
    ((M) <= FLA_GESVD_SMALL_SIZE_THRESH1 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))

/* Path 10: M at least N, but not much larger */
#define FLA_GESVD_SMALL_PATH10(WNTUN, WNTUS, WNTVN, WNTVS, M)                           \
    (((WNTUN) || (WNTUS)) && ((WNTVN) || (WNTVS)) && (M) < FLA_GESVD_SMALL_SIZE_THRESH1 \
     && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))

/* Path 1t: N much larger than M, JOBVT='N' */
#define FLA_GESVD_SMALL_PATH1T(WNTUN, WNTVN, M, N)             \
    ((WNTUN) && (WNTVN) && (N) <= FLA_GESVD_SMALL_SIZE_THRESH2 \
     && (M) < FLA_GESVD_SMALL_SIZE_THRESH0 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))

/* Path 6t: N much larger than M, JOBU='S' or 'A', JOBVT='S' */
#define FLA_GESVD_SMALL_PATH6T(N) \
    ((N) <= FLA_GESVD_SMALL_SIZE_THRESH1 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))

/* Path 10t: N greater than M, but not much larger */
#define FLA_GESVD_SMALL_PATH10T(WNTUAS, WNTVS, N)               \
    (((WNTUAS) & (WNTVS)) && (N) < FLA_GESVD_SMALL_SIZE_THRESH0 \
     && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))

#else

#define FLA_GESVD_SMALL_PATH1(WNTVO, M) (0)
#define FLA_GESVD_SMALL_PATH6(M) (0)
#define FLA_GESVD_SMALL_PATH10(WNTUN, WNTUS, WNTVN, WNTVS, M) (0)
#define FLA_GESVD_SMALL_PATH1T(WNTUN, WNTVN, M, N) (0)
#define FLA_GESVD_SMALL_PATH6T(N) (0)
#define FLA_GESVD_SMALL_PATH10T(WNTUAS, WNTVS, N) (0)

#endif /* FLA_ENABLE_AMD_OPT */
#endif /* FLA_GESVD_DISPATCH_H */
