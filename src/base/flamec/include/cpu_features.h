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
#ifndef FLA_CPUFEATURES_H
#define FLA_CPUFEATURES_H

#include <assert.h>
#include <stdint.h>
#include <string.h>

enum {
    ALC_CPUID_EAX_1 = 0,
    ALC_CPUID_EAX_7,
    ALC_CPUID_EAX_8_01,           /* 8000.0001 */
    ALC_CPUID_EAX_8_07,           /* 8000.0007 */
    ALC_CPUID_EAX_8_08,           /* 8000.0008 */

    /* Last entry */
    ALC_CPUID_MAX,
};

enum {
    /*EBX Values*/
    ALC_CPUID_BIT_FSGSBASE        = (1u << 0),
    ALC_CPUID_BIT_TSC_ADJUST      = (1u << 1),
    ALC_CPUID_BIT_SGX             = (1u << 2),
    ALC_CPUID_BIT_BMI1            = (1u << 3),
    ALC_CPUID_BIT_HLE             = (1u << 4),
    ALC_CPUID_BIT_AVX2            = (1u << 5),
    ALC_CPUID_BIT_SMEP            = (1u << 7),
    ALC_CPUID_BIT_BMI2            = (1u << 8),
    ALC_CPUID_BIT_ERMS            = (1u << 9),
    ALC_CPUID_BIT_INVPCID         = (1u << 10),
    ALC_CPUID_BIT_RTM             = (1u << 11),
    ALC_CPUID_BIT_TSX             = ALC_CPUID_BIT_RTM,
    ALC_CPUID_BIT_PQM             = (1u << 12),
    ALC_CPUID_BIT_MPX             = (1u << 14),
    ALC_CPUID_BIT_PQE             = (1u << 15),
    ALC_CPUID_BIT_AVX512F         = (1u << 16),
    ALC_CPUID_BIT_AVX512DQ        = (1u << 17),
    ALC_CPUID_BIT_RDSEED          = (1u << 18),
    ALC_CPUID_BIT_ADX             = (1u << 19),
    ALC_CPUID_BIT_SMAP            = (1u << 20),
    ALC_CPUID_BIT_AVX512_IFMA     = (1u << 21),
    ALC_CPUID_BIT_CLFLUSHOPT      = (1u << 22),
    ALC_CPUID_BIT_CLWB            = (1u << 24),
    ALC_CPUID_BIT_TRACE           = (1u << 25),
    ALC_CPUID_BIT_AVX512PF        = (1u << 26),
    ALC_CPUID_BIT_AVX512ER        = (1u << 27),
    ALC_CPUID_BIT_AVX512CD        = (1u << 28),
    ALC_CPUID_BIT_SHA             = (1u << 29),
    ALC_CPUID_BIT_AVX512BW        = (1u << 30),
    ALC_CPUID_BIT_AVX512VL        = (1u << 31),

    /* ECX Values*/
    ALC_CPUID_BIT_PREFETCHWT1     = (1u << 0),
    ALC_CPUID_BIT_AVX512_VBMI     = (1u << 1),
    ALC_CPUID_BIT_UMIP            = (1u << 2),
    ALC_CPUID_BIT_PKU             = (1u << 3),
    ALC_CPUID_BIT_OSPKE           = (1u << 4),
    ALC_CPUID_BIT_WAITPKG         = (1u << 5),
    ALC_CPUID_BIT_AVX512_VBMI2    = (1u << 6),
    ALC_CPUID_BIT_SHSTK           = (1u << 7),
    ALC_CPUID_BIT_GFNI            = (1u << 8),
    ALC_CPUID_BIT_VAES            = (1u << 9),
    ALC_CPUID_BIT_VPCLMULQDQ      = (1u << 10),
    ALC_CPUID_BIT_AVX512_VNNI     = (1u << 11),
    ALC_CPUID_BIT_AVX512_BITALG   = (1u << 12),
    ALC_CPUID_BIT_AVX512_VPOPCNTDQ = (1u << 14),
    ALC_CPUID_BIT_RDPID           = (1u << 22),
    ALC_CPUID_BIT_CLDEMOTE        = (1u << 25),
    ALC_CPUID_BIT_MOVDIRI         = (1u << 27),
    ALC_CPUID_BIT_MOVDIR64B       = (1u << 28),
    ALC_CPUID_BIT_SGX_LC          = (1u << 30),

    /* EDX Values */
    ALC_CPUID_BIT_AVX512_4VNNIW   = (1u << 2),
    ALC_CPUID_BIT_AVX512_4FMAPS   = (1u << 3),
    ALC_CPUID_BIT_FSRM            = (1u << 4),
    ALC_CPUID_BIT_PCONFIG         = (1u << 18),
    ALC_CPUID_BIT_IBT             = (1u << 20),
    ALC_CPUID_BIT_IBRS_IBPB       = (1u << 26),
    ALC_CPUID_BIT_STIBP           = (1u << 27),
    ALC_CPUID_BIT_CAPABILITIES    = (1u << 29),
    ALC_CPUID_BIT_SSBD            = (1u << 31),
};

#define ALC_CPU_FAMILY_ZEN              0x17
#define ALC_CPU_FAMILY_ZEN_PLUS         0x17
#define ALC_CPU_FAMILY_ZEN2             0x17
#define ALC_CPU_FAMILY_ZEN3             0x19
#define ALC_CPU_FAMILY_ZEN4             0x19

static inline uint32_t
__extract32(uint32_t value, int start, int length);

static inline uint16_t
alc_cpuid_get_family(uint32_t var);

static inline uint16_t
alc_cpuid_get_model(uint32_t var);

static inline uint16_t
alc_cpuid_get_stepping(uint32_t var);

/* ID return values */
struct alc_cpuid_regs {
    uint32_t eax;
    uint32_t ebx;
    uint32_t ecx;
    uint32_t edx;
};

typedef enum {
    ALC_CPU_MFG_INTEL,
    ALC_CPU_MFG_AMD,
    ALC_CPU_MFG_OTHER,
} alc_cpu_mfg_t;

struct alc_cpu_mfg_info {
    alc_cpu_mfg_t     mfg_type;
    uint16_t          family;
    uint16_t          model;
    uint16_t          stepping;
};

struct alc_cpu_features {
    struct alc_cpu_mfg_info cpu_mfg_info;
    struct alc_cpuid_regs   available[ALC_CPUID_MAX];
    struct alc_cpuid_regs   usable[ALC_CPUID_MAX];
};

static inline void __cpuid(struct alc_cpuid_regs *out);

static inline void __cpuid_1(uint32_t eax, struct alc_cpuid_regs *out);

static inline void __cpuid_2(uint32_t eax, uint32_t ecx, struct alc_cpuid_regs *out);


#ifndef ARRAY_SIZE
#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))
#endif

#define INITIALIZED_MAGIC 0xdeadbeaf


static void
__get_mfg_info(struct alc_cpuid_regs*   cpuid_regs,
               struct alc_cpu_mfg_info* mfg_info);

static void
__init_cpu_features(void);


uint32_t
alc_cpu_has_avx512f(void);

uint32_t
alc_cpu_has_avx512dq(void);

uint32_t
alc_cpu_has_avx512bw(void);

uint32_t
alc_cpu_has_avx512er(void);

uint32_t
alc_cpu_has_avx512cd(void);

uint32_t
alc_cpu_has_avx512vl(void);

uint32_t
alc_cpu_has_avx512pf(void);

uint32_t
alc_cpu_has_avx512_ifma(void);

uint32_t
alc_cpu_has_avx512_vnni(void);

uint32_t
alc_cpu_has_avx512_bitalg(void);

uint32_t
alc_cpu_has_avx512_vbmi(void);

uint32_t
alc_cpu_has_avx512_vbmi2(void);

uint32_t
alc_cpu_has_avx512_vpopcntdq(void);

#endif //FLA_CPUFEATURES_H
