/**
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */

#include <stdint.h>

#ifdef INT_64BIT
typedef int64_t aocl_int_t;
#else
typedef int32_t aocl_int_t;
#endif

int64_t aocl_lapack_ilaenv(int64_t *ispec, char *name__, char *opts, int64_t *n1, int64_t *n2,
                           int64_t *n3, int64_t *n4)
{
    /* System generated locals */
    extern aocl_int_t ilaenv_(aocl_int_t *, char *, char *, aocl_int_t *, aocl_int_t *,
                              aocl_int_t *, aocl_int_t *);

    aocl_int_t ispec_lp = (aocl_int_t)(*ispec);
    aocl_int_t n1_lp = (aocl_int_t)(*n1);
    aocl_int_t n2_lp = (aocl_int_t)(*n2);
    aocl_int_t n3_lp = (aocl_int_t)(*n3);
    aocl_int_t n4_lp = (aocl_int_t)(*n4);

    return ilaenv_(&ispec_lp, name__, opts, &n1_lp, &n2_lp, &n3_lp, &n4_lp);

} /* ilaenv_ */

int64_t aocl_lapack_ilaenv2stage(int64_t *ispec, char *name__, char *opts, int64_t *n1, int64_t *n2,
                                 int64_t *n3, int64_t *n4)
{
    extern aocl_int_t ilaenv2stage_(aocl_int_t *, char *, char *, aocl_int_t *, aocl_int_t *,
                                    aocl_int_t *, aocl_int_t *);

    aocl_int_t ispec_lp = (aocl_int_t)(*ispec);
    aocl_int_t n1_lp = (aocl_int_t)(*n1);
    aocl_int_t n2_lp = (aocl_int_t)(*n2);
    aocl_int_t n3_lp = (aocl_int_t)(*n3);
    aocl_int_t n4_lp = (aocl_int_t)(*n4);

    return ilaenv2stage_(&ispec_lp, name__, opts, &n1_lp, &n2_lp, &n3_lp, &n4_lp);

} /* aocl_lapack_ilaenv2stage */

int64_t aocl_lapack_iparmq(int64_t *ispec, char *name__, char *opts, int64_t *n, int64_t *ilo,
                           int64_t *ihi, int64_t *lwork)
{
    extern aocl_int_t iparmq_(aocl_int_t *, char *, char *, aocl_int_t *, aocl_int_t *,
                              aocl_int_t *, aocl_int_t *);

    aocl_int_t ispec_lp = (aocl_int_t)(*ispec);
    aocl_int_t n_lp = (aocl_int_t)(*n);
    aocl_int_t ilo_lp = (aocl_int_t)(*ilo);
    aocl_int_t ihi_lp = (aocl_int_t)(*ihi);
    aocl_int_t lwork_lp = (aocl_int_t)(*lwork);

    return iparmq_(&ispec_lp, name__, opts, &n_lp, &ilo_lp, &ihi_lp, &lwork_lp);
} /* aocl_lapack_iparmq */
