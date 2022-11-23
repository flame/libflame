/* ************************************************************************
 * Copyright (c) 2022 Advanced Micro Devices, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * ************************************************************************ */

#include "FLAME.h"

#define ALC_CPU_FEATURE_REG(ftr, idx, reg) ({   \
            uint32_t val;                       \
            struct alc_cpuid_regs *r;               \
            r = &(ftr)->available[0];             \
            val = r[(idx)].reg;                 \
            val;                                \
        })

#define ALC_CPU_FEATURE(ptr, idx, reg, bit) ({      \
            uint32_t __reg =                        \
                ALC_CPU_FEATURE_REG(ptr, idx, reg); \
            (__reg & bit);                          \
        })

/* For AVX512 instructions */
#define ALC_CPU_HAS_AVX512F(f)     ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ebx, ALC_CPUID_BIT_AVX512F) /* For AVX512 foundation flag */
#define ALC_CPU_HAS_AVX512DQ(f)    ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ebx, ALC_CPUID_BIT_AVX512DQ)
#define ALC_CPU_HAS_AVX512BW(f)    ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ebx, ALC_CPUID_BIT_AVX512BW)
#define ALC_CPU_HAS_AVX512ER(f)    ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ebx, ALC_CPUID_BIT_AVX512ER)
#define ALC_CPU_HAS_AVX512CD(f)    ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ebx, ALC_CPUID_BIT_AVX512CD)
#define ALC_CPU_HAS_AVX512VL(f)    ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ebx, ALC_CPUID_BIT_AVX512VL)
#define ALC_CPU_HAS_AVX512PF(f)    ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ebx, ALC_CPUID_BIT_AVX512PF)
#define ALC_CPU_HAS_AVX512_IFMA(f) ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ebx, ALC_CPUID_BIT_AVX512_IFMA)

#define ALC_CPU_HAS_AVX512_VNNI(f)      ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ecx, ALC_CPUID_BIT_AVX512_VNNI)
#define ALC_CPU_HAS_AVX512_BITALG(f)    ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ecx, ALC_CPUID_BIT_AVX512_BITALG)
#define ALC_CPU_HAS_AVX512_VBMI(f)      ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ecx, ALC_CPUID_BIT_AVX512_VBMI)
#define ALC_CPU_HAS_AVX512_VBMI2(f)     ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ecx, ALC_CPUID_BIT_AVX512_VBMI2)
#define ALC_CPU_HAS_AVX512_VPOPCNTDQ(f) ALC_CPU_FEATURE(f, ALC_CPUID_EAX_7, ecx, ALC_CPUID_BIT_AVX512_VPOPCNTDQ)

struct alc_cpu_features cpu_features;

struct
{
    uint32_t eax;
    uint32_t ecx;
} __cpuid_values[ALC_CPUID_MAX] = {
    [ALC_CPUID_EAX_1]    = { 0x1, 0x0 },        /* eax = 0, ecx=0 */
    [ALC_CPUID_EAX_7]    = { 0x7, 0x0 },        /* eax = 7,  -"- */
    [ALC_CPUID_EAX_8_01] = { 0x80000001, 0x0 }, /* eax = 0x80000001 */
    [ALC_CPUID_EAX_8_07] = { 0x80000007, 0x0 }, /* eax = 0x80000007 */
    [ALC_CPUID_EAX_8_08] = { 0x80000008, 0x0 }, /* eax = 0x80000008 */
};

static inline uint32_t
__extract32(uint32_t value, int start, int length)
{
    assert(start >= 0 && length > 0 && length <= 32 - start);
    return (value >> start) & (~0U >> (32 - length));
}


static inline uint16_t
alc_cpuid_get_family(uint32_t var)
{
    return (uint16_t)(__extract32(var, 20, 8) +
                      __extract32(var, 8, 4));
}

static inline uint16_t
alc_cpuid_get_model(uint32_t var)
{
    return (uint16_t)(__extract32(var, 16, 4) << 4 |
                      __extract32(var, 4, 4));
}

static inline uint16_t
alc_cpuid_get_stepping(uint32_t var)
{
    return (uint16_t)(__extract32(var, 20, 8) +
                      __extract32(var, 8, 4));
}

static inline void __cpuid(struct alc_cpuid_regs *out)
{
    __asm__ volatile
        (
         "cpuid"
         :"=a"(out->eax), "=b"(out->ebx), "=c"(out->ecx), "=d"(out->edx)
         );
}

static inline void __cpuid_1(uint32_t eax, struct alc_cpuid_regs *out)
{
    __asm__ volatile
        (
         "cpuid"
         :"=a"(out->eax), "=b"(out->ebx), "=c"(out->ecx), "=d"(out->edx)
         :"0"(eax)
         );
}

static inline void __cpuid_2(uint32_t eax, uint32_t ecx, struct alc_cpuid_regs *out)
{
    __asm__ volatile
        (
         "cpuid"
         :"=a"(out->eax), "=b"(out->ebx), "=c"(out->ecx), "=d"(out->edx)
         :"0"(eax), "2"(ecx)
         );
}

static void
__get_mfg_info(struct alc_cpuid_regs*   cpuid_regs,
               struct alc_cpu_mfg_info* mfg_info)
{
    uint16_t model;
    uint16_t family;

    if (mfg_info) {
        struct alc_cpuid_regs regs;

        __cpuid_1(1, &regs);

        family = alc_cpuid_get_family(regs.eax);
        model  = alc_cpuid_get_model(regs.eax);

        if (family >= ALC_CPU_FAMILY_ZEN) {
            mfg_info->family = (uint16_t)family;
            mfg_info->model  = (uint16_t)model;
        }

        mfg_info->stepping = alc_cpuid_get_stepping(regs.eax);
    }
}

static void
__init_cpu_features(void)
{
    static unsigned initialized = 0;

    struct alc_cpu_mfg_info* mfg_info = &cpu_features.cpu_mfg_info;
    int                      arr_size = ARRAY_SIZE(__cpuid_values);

    if (initialized == INITIALIZED_MAGIC)
        return;

    struct alc_cpuid_regs regs;
    __cpuid_1(0, &regs);

    /* "AuthenticAMD" */
    if (regs.ebx == 0x68747541 && regs.ecx == 0x444d4163
        && regs.edx == 0x69746e65) {
        cpu_features.cpu_mfg_info.mfg_type = ALC_CPU_MFG_AMD;
    }

    for (int i = 0; i < arr_size; i++) {
        struct alc_cpuid_regs ft;

        __cpuid_2(__cpuid_values[i].eax, __cpuid_values[i].ecx, &ft);

        cpu_features.available[i].eax = ft.eax;
        cpu_features.available[i].ebx = ft.ebx;
        cpu_features.available[i].ecx = ft.ecx;
        cpu_features.available[i].edx = ft.edx;
    }

    __get_mfg_info(&cpu_features.available[ALC_CPUID_EAX_1], mfg_info);

    /*
     * Globally disable some *_USEABLE flags, so that all ifunc's
     * sees them
     */
    if (mfg_info->mfg_type == ALC_CPU_MFG_AMD
        && mfg_info->family >= ALC_CPU_FAMILY_ZEN) {
        memcpy(&cpu_features.usable[0],
               &cpu_features.available[0],
               sizeof(cpu_features.usable));

        }
    initialized = INITIALIZED_MAGIC;
}



uint32_t
alc_cpu_has_avx512f(void)
{
    __init_cpu_features();

    return ALC_CPU_HAS_AVX512F(&cpu_features);
}

uint32_t
alc_cpu_has_avx512dq(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512DQ(&cpu_features);
}

uint32_t
alc_cpu_has_avx512bw(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512BW(&cpu_features);
}

uint32_t
alc_cpu_has_avx512er(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512ER(&cpu_features);
}

uint32_t
alc_cpu_has_avx512cd(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512CD(&cpu_features);
}

uint32_t
alc_cpu_has_avx512vl(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512VL(&cpu_features);
}

uint32_t
alc_cpu_has_avx512pf(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512PF(&cpu_features);
}

uint32_t
alc_cpu_has_avx512_ifma(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512_IFMA(&cpu_features);
}

uint32_t
alc_cpu_has_avx512_vnni(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512_VNNI(&cpu_features);
}

uint32_t
alc_cpu_has_avx512_bitalg(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512_BITALG(&cpu_features);
}

uint32_t
alc_cpu_has_avx512_vbmi(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512_VBMI(&cpu_features);
}

uint32_t
alc_cpu_has_avx512_vbmi2(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512_VBMI2(&cpu_features);
}

uint32_t
alc_cpu_has_avx512_vpopcntdq(void)
{
    if (alc_cpu_has_avx512f() == 0)
        return 0;

    __init_cpu_features();

    return ALC_CPU_HAS_AVX512_VPOPCNTDQ(&cpu_features);
}

