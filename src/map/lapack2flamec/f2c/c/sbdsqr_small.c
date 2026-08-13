/*
 *     Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
#include "FLA_f2c.h" /* Table of constant values */
#include "fla_bdsqr_small_defs.h"

#if FLA_ENABLE_AMD_OPT
static doublereal c_b15 = -.125;
static aocl_int64_t c__1 = 1;
static real c_b49 = 1.f;
static real c_b72 = -1.f;

#define BDSQR_PRE s
#include "fla_bdsqr_small_kernel.h"

#endif
