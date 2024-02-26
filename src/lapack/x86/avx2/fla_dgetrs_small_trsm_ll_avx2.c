/******************************************************************************
 * Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dgetrs_small_trsm_ll_avx2.c
 *  @brief solves a system of linear equations: A * X = B  or  A**T * X = B in AVX2.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

#define TRANSPOSE_4x4(r_1, r_2, r_3, r_4)                              \
    t_reg[0] = _mm256_unpacklo_pd(b_reg[r_1], b_reg[r_2]);             \
    t_reg[1] = _mm256_unpacklo_pd(b_reg[r_3], b_reg[r_4]);             \
    b_reg[r_1] = _mm256_unpackhi_pd(b_reg[r_1], b_reg[r_2]);           \
    b_reg[r_2] = _mm256_unpackhi_pd(b_reg[r_3], b_reg[r_4]);           \
    b_reg[r_4] = _mm256_permute2f128_pd(b_reg[r_1], b_reg[r_2], 0x31); \
    b_reg[r_2] = _mm256_permute2f128_pd(b_reg[r_1], b_reg[r_2], 0x20); \
    b_reg[r_1] = _mm256_permute2f128_pd(t_reg[0], t_reg[1], 0x20);     \
    b_reg[r_3] = _mm256_permute2f128_pd(t_reg[0], t_reg[1], 0x31);

#define TRANSPOSE_8x8()                 \
    TRANSPOSE_4x4(0, 1, 2, 3)           \
        TRANSPOSE_4x4(4, 5, 6, 7)       \
            TRANSPOSE_4x4(8, 9, 10, 11) \
                TRANSPOSE_4x4(12, 13, 14, 15)

static const __m256i mask_reg[6] = {{0, 0, 0, 0},
                                    {-1, 0, 0, 0},
                                    {-1, -1, 0, 0},
                                    {-1, -1, -1, 0},
                                    {-1, -1, -1, -1},
                                    {-1, -1, -1, -1}};

static void n_8(integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *ipiv)
{
    __m256d b_reg[16];
    __m256i mask1, mask2;
    __m256d t_reg[3];
    mask1 = mask_reg[4];
    mask2 = mask_reg[4];
    int i, j;
    /*****************************************/
    /* B matrix is loaded into the folliwing */
    /* b_reg registers.                      */
    /*example: col 3, row 4-7 are loaded     */
    /*         into b_reg[7]                 */
    /*   n->  0  1  2  3  4  5  6  7         */
    /*       -_  _  _  _  _  _  _  _ -       */
    /*    0  |0  1  2  3  8  9  10 11|       */
    /*    1  |0  1  2  3  8  9  10 11|       */
    /*    2  |0  1  2  3  8  9  10 11|       */
    /*    3  |0  1  2  3  8  9  10 11|       */
    /*    4  |4  5  6  7  12 13 14 15|       */
    /*    5  |4  5  6  7  12 13 14 15|       */
    /*    6  |4  5  6  7  12 13 14 15|       */
    /*    7  |4  5  6  7  12 13 14 15|       */
    /*       -_  _  _  _  _  _  _  _ -       */
    /*                                       */
    /*****************************************/
    for (i = 0; i < (*nrhs); ++i)
    {
        j = (int)(i / 4);
        j *= 4;
        b_reg[i + 0 + j] = _mm256_maskload_pd((void const *)(b + (i * (*ldb)) + 0), mask1);
        b_reg[i + 4 + j] = _mm256_maskload_pd((void const *)(b + (i * (*ldb)) + 4), mask2);
    }
    /* To vectorize 'left' variants of TRSM, B matrix is required to be */
    /* stored in row major format, to convert it to row major format,   */
    /* B matrix is transposed.                                          */
 
 
    /*****************************************/
    /* After Transpose the registers are     */
    /* changed to following order.           */
    /*                                       */
    /*   n->  0  1  2  3  4  5  6  7         */
    /*       -_  _  _  _  _  _  _  _ -       */
    /*    0  |0  0  0  0  8  8  8  8 |       */
    /*    1  |1  1  1  1  9  9  9  9 |       */
    /*    2  |2  2  2  2  10 10 10 10|       */
    /*    3  |3  3  3  3  11 11 11 11|       */
    /*    4  |4  4  4  4  12 12 12 12|       */
    /*    5  |5  5  5  5  13 13 13 13|       */
    /*    6  |6  6  6  6  14 14 14 14|       */
    /*    7  |7  7  7  7  15 15 15 15|       */
    /*       -_  _  _  _  _  _  _  _ -       */
    /*                                       */
    /*****************************************/
    TRANSPOSE_8x8()

    // # REGION - ROW Swap
 
    /* After tranpose, B matrix is stored in row major format in the b_reg registers*/
    /* so in order to swap row, we only need to swap registers                      */
    /* Swap Row [n] with Row [Ipiv[n]]                                              */
    t_reg[0] = b_reg[0], t_reg[1] = b_reg[8];                        // store row 0(b_reg[0] and b_reg[8]) into temporary registers
    b_reg[0] = b_reg[ipiv[0] - 1], b_reg[8] = b_reg[ipiv[0] + 7];    // copy row [ipiv[0]] into row 0 registers
    b_reg[ipiv[0] - 1] = t_reg[0], b_reg[ipiv[0] + 7] = t_reg[1];    // copy row 0(from temp registers) to row [ipiv[0]]

    t_reg[0] = b_reg[1], t_reg[1] = b_reg[9];
    b_reg[1] = b_reg[ipiv[1] - 1], b_reg[9] = b_reg[ipiv[1] + 7];
    b_reg[ipiv[1] - 1] = t_reg[0], b_reg[ipiv[1] + 7] = t_reg[1];

    t_reg[0] = b_reg[2], t_reg[1] = b_reg[10];
    b_reg[2] = b_reg[ipiv[2] - 1], b_reg[10] = b_reg[ipiv[2] + 7];
    b_reg[ipiv[2] - 1] = t_reg[0], b_reg[ipiv[2] + 7] = t_reg[1];

    t_reg[0] = b_reg[3], t_reg[1] = b_reg[11];
    b_reg[3] = b_reg[ipiv[3] - 1], b_reg[11] = b_reg[ipiv[3] + 7];
    b_reg[ipiv[3] - 1] = t_reg[0], b_reg[ipiv[3] + 7] = t_reg[1];

    t_reg[0] = b_reg[4], t_reg[1] = b_reg[12];
    b_reg[4] = b_reg[ipiv[4] - 1], b_reg[12] = b_reg[ipiv[4] + 7];
    b_reg[ipiv[4] - 1] = t_reg[0], b_reg[ipiv[4] + 7] = t_reg[1];

    t_reg[0] = b_reg[5], t_reg[1] = b_reg[13];
    b_reg[5] = b_reg[ipiv[5] - 1], b_reg[13] = b_reg[ipiv[5] + 7];
    b_reg[ipiv[5] - 1] = t_reg[0], b_reg[ipiv[5] + 7] = t_reg[1];

    t_reg[0] = b_reg[6], t_reg[1] = b_reg[14];
    b_reg[6] = b_reg[ipiv[6] - 1], b_reg[14] = b_reg[ipiv[6] + 7];
    b_reg[ipiv[6] - 1] = t_reg[0], b_reg[ipiv[6] + 7] = t_reg[1];

    t_reg[0] = b_reg[7], t_reg[1] = b_reg[15];
    b_reg[7] = b_reg[ipiv[7] - 1], b_reg[15] = b_reg[ipiv[7] + 7];
    b_reg[ipiv[7] - 1] = t_reg[0], b_reg[ipiv[7] + 7] = t_reg[1];
    // # ENDREGION - ROW Swap
 
    // # REGION - TRSM LLNU
 
    // ROW 0 compute is not needed because diagonal is unit.
    // REGION - TRSM row 1 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1)); // t_reg[2] = [ a[1][0], a[1][0], a[1][0], a[1][0] ]
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]); // t_reg[0] = b_reg[0] * a[1][0]
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]); // t_reg[1] = b_reg[8] * a[1][0]  // row0 * a[1][0]
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]); // b_reg[1] -= t_teg[0]
    b_reg[9] = _mm256_sub_pd(b_reg[9], t_reg[1]); // b_reg[9] -= t_teg[1]   //row1 -= row0 * a[1][0]
    // ENDREGION - TRSM row 1 computation
 
    // REGION - TRSM row 2 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]); // t_reg[0 to 1] = row0 * a[2][0]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]); // t_reg[0 to 1] = row0 * a[2][0] + row1 * a[2][1]
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);
    b_reg[10] = _mm256_sub_pd(b_reg[10], t_reg[1]); // row2 -= row0 * a[2][0] + row1 * a[2][1]
    // ENDREGION - TRSM row 2 computation
 
    // REGION - TRSM row 3 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]); // t_reg[0 to 1] = row0 * a[3][0]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]); // t_reg[0 to 1] = row0 * a[3][0] + row1 * a[3][1]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]); // t_reg[0 to 1] = row0 * a[3][0] + row1 * a[3][1] + row2 * a[3][2]
    b_reg[3] = _mm256_sub_pd(b_reg[3], t_reg[0]);
    b_reg[11] = _mm256_sub_pd(b_reg[11], t_reg[1]); // row3 -= row0 * a[3][0] + row1 * a[3][1] + row2 * a[3][2]
    // ENDREGION - TRSM row 3 computation
 
    // REGION - TRSM row 4 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]); // t_reg[0 to 1] = row0 * a[4][0]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]); // t_reg[0 to 1] += row1 * a[4][1]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]); // t_reg[0 to 1] += row2 * a[4][2]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]); // t_reg[0 to 1] += row3 * a[4][3]
    b_reg[4] = _mm256_sub_pd(b_reg[4], t_reg[0]);
    b_reg[12] = _mm256_sub_pd(b_reg[12], t_reg[1]); //row4 -= t_reg [0 to 1]
    // ENDREGION - TRSM row 4 computation
 
    // REGION - TRSM row 5 computation

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    b_reg[5] = _mm256_sub_pd(b_reg[5], t_reg[0]);
    b_reg[13] = _mm256_sub_pd(b_reg[13], t_reg[1]);
    // ENDREGION - TRSM row 5 computation
 
    // REGION - TRSM row 6 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    b_reg[6] = _mm256_sub_pd(b_reg[6], t_reg[0]);
    b_reg[14] = _mm256_sub_pd(b_reg[14], t_reg[1]);
    // ENDREGION - TRSM row 6 computation
 
    // REGION - TRSM row 7 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 7));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 7 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 7 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 7 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 7 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 7 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 7 + (6 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[6], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[14], t_reg[2], t_reg[1]);
    b_reg[7] = _mm256_sub_pd(b_reg[7], t_reg[0]);
    b_reg[15] = _mm256_sub_pd(b_reg[15], t_reg[1]);

    // ENDREGION - TRSM row 7 computation
 
    // # ENDREGION - TRSM LLNU
 
    // # REGION - TRSM LUNN
 
    // REGION - TRSM row 7 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 7 + (7 * (*lda))));
    b_reg[7] = _mm256_div_pd(b_reg[7], t_reg[2]);
    b_reg[15] = _mm256_div_pd(b_reg[15], t_reg[2]); // row7 /= a[7][7]
    // ENDREGION - TRSM row 7 computation
 
    // REGION - TRSM row 6 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (7 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[7], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[15], t_reg[2]); // t_reg[0 to 1] = row7 * a[6][7]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (6 * (*lda))));
    b_reg[6] = _mm256_sub_pd(b_reg[6], t_reg[0]);
    b_reg[6] = _mm256_div_pd(b_reg[6], t_reg[2]);
    b_reg[14] = _mm256_sub_pd(b_reg[14], t_reg[1]); // row6 -= t_reg[0 to 1]
    b_reg[14] = _mm256_div_pd(b_reg[14], t_reg[2]); // row6 /= a[6][6]
    // ENDREGION - TRSM row 6 computation
 
    // REGION - TRSM row 5 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (7 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[7], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[15], t_reg[2]); // t_reg[0 to 1] = row7 * a[5][7]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (6 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[6], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[14], t_reg[2], t_reg[1]); // t_reg[0 to 1] += row6 * a[5][6]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (5 * (*lda))));
    b_reg[5] = _mm256_sub_pd(b_reg[5], t_reg[0]);
    b_reg[5] = _mm256_div_pd(b_reg[5], t_reg[2]);
    b_reg[13] = _mm256_sub_pd(b_reg[13], t_reg[1]); // row5 -= t_reg[0 to 1]
    b_reg[13] = _mm256_div_pd(b_reg[13], t_reg[2]); // row5 /= a[5][5]
    // ENDREGION - TRSM row 5 computation
 
    // REGION - TRSM row 4 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (7 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[7], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[15], t_reg[2]); // t_reg[0 to 1] = row7 * a[4][7]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (6 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[6], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[14], t_reg[2], t_reg[1]); // t_reg[0 to 1] += row6 * a[4][6]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]); // t_reg[0 to 1] += row5 * a[4][5]
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (4 * (*lda))));
    b_reg[4] = _mm256_sub_pd(b_reg[4], t_reg[0]);
    b_reg[4] = _mm256_div_pd(b_reg[4], t_reg[2]);
    b_reg[12] = _mm256_sub_pd(b_reg[12], t_reg[1]); // row4 -= t_reg[0 to 1]
    b_reg[12] = _mm256_div_pd(b_reg[12], t_reg[2]); // row4 /= a[4][4]
    // ENDREGION - TRSM row 4 computation
 
    // REGION - TRSM row 3 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (7 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[7], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[15], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (6 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[6], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[14], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (3 * (*lda))));
    b_reg[3] = _mm256_sub_pd(b_reg[3], t_reg[0]);
    b_reg[3] = _mm256_div_pd(b_reg[3], t_reg[2]);
    b_reg[11] = _mm256_sub_pd(b_reg[11], t_reg[1]);
    b_reg[11] = _mm256_div_pd(b_reg[11], t_reg[2]);
    // ENDREGION - TRSM row 3 computation
 
    // REGION - TRSM row 2 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (7 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[7], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[15], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (6 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[6], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[14], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (2 * (*lda))));
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);
    b_reg[2] = _mm256_div_pd(b_reg[2], t_reg[2]);
    b_reg[10] = _mm256_sub_pd(b_reg[10], t_reg[1]);
    b_reg[10] = _mm256_div_pd(b_reg[10], t_reg[2]);
    // ENDREGION - TRSM row 2 computation
 
    // REGION - TRSM row 1 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (7 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[7], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[15], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (6 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[6], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[14], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (1 * (*lda))));
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    b_reg[1] = _mm256_div_pd(b_reg[1], t_reg[2]);
    b_reg[9] = _mm256_sub_pd(b_reg[9], t_reg[1]);
    b_reg[9] = _mm256_div_pd(b_reg[9], t_reg[2]);
    // ENDREGION - TRSM row 1 computation
 
    // REGION - TRSM row 0 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (7 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[7], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[15], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (6 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[6], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[14], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a));
    b_reg[0] = _mm256_sub_pd(b_reg[0], t_reg[0]);
    b_reg[0] = _mm256_div_pd(b_reg[0], t_reg[2]);
    b_reg[8] = _mm256_sub_pd(b_reg[8], t_reg[1]);
    b_reg[8] = _mm256_div_pd(b_reg[8], t_reg[2]);
    // ENDREGION - TRSM row 0 computation
    // # ENDREGION - TRSM LUNN
 
    // transpose B matrix again to convert it to column major format
    /*****************************************/
    /*Order of registers after transpose     */
    /*                                       */
    /*   n->  0  1  2  3  4  5  6  7         */
    /*       -_  _  _  _  _  _  _  _ -       */
    /*    0  |0  1  2  3  8  9  10 11|       */
    /*    1  |0  1  2  3  8  9  10 11|       */
    /*    2  |0  1  2  3  8  9  10 11|       */
    /*    3  |0  1  2  3  8  9  10 11|       */
    /*    4  |4  5  6  7  12 13 14 15|       */
    /*    5  |4  5  6  7  12 13 14 15|       */
    /*    6  |4  5  6  7  12 13 14 15|       */
    /*    7  |4  5  6  7  12 13 14 15|       */
    /*       -_  _  _  _  _  _  _  _ -       */
    /*                                       */
    /*****************************************/
    TRANSPOSE_8x8()
    for (i = 0; i < (*nrhs); ++i)
    {
        j = (int)(i / 4);
        j *= 4;
        _mm256_maskstore_pd((b + (i * (*ldb)) + 0), mask1, b_reg[i + 0 + j]);
        _mm256_maskstore_pd((b + (i * (*ldb)) + 4), mask2, b_reg[i + 4 + j]);
    }
}

static void n_7(integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *ipiv)
{
    __m256d b_reg[16];
    __m256i mask1, mask2;
    __m256d t_reg[3];
    mask1 = mask_reg[4];
    mask2 = mask_reg[3];
    int i, j;
    for (i = 0; i < (*nrhs); ++i)
    {
        j = (int)(i / 4);
        j *= 4;
        b_reg[i + 0 + j] = _mm256_maskload_pd((void const *)(b + (i * (*ldb)) + 0), mask1);
        b_reg[i + 4 + j] = _mm256_maskload_pd((void const *)(b + (i * (*ldb)) + 4), mask2);
    }
    TRANSPOSE_8x8()
    /* After tranpose, B matrix is stored in row major format in the b_reg registers*/
    /* so in order to swap row, we only need to swap registers                      */
    /* Swap Row [n] with Row [Ipiv[n]]                                              */
    t_reg[0] = b_reg[0],
    t_reg[1] = b_reg[8];
    b_reg[0] = b_reg[ipiv[0] - 1], b_reg[8] = b_reg[ipiv[0] + 7];
    b_reg[ipiv[0] - 1] = t_reg[0], b_reg[ipiv[0] + 7] = t_reg[1];

    t_reg[0] = b_reg[1], t_reg[1] = b_reg[9];
    b_reg[1] = b_reg[ipiv[1] - 1], b_reg[9] = b_reg[ipiv[1] + 7];
    b_reg[ipiv[1] - 1] = t_reg[0], b_reg[ipiv[1] + 7] = t_reg[1];

    t_reg[0] = b_reg[2], t_reg[1] = b_reg[10];
    b_reg[2] = b_reg[ipiv[2] - 1], b_reg[10] = b_reg[ipiv[2] + 7];
    b_reg[ipiv[2] - 1] = t_reg[0], b_reg[ipiv[2] + 7] = t_reg[1];

    t_reg[0] = b_reg[3], t_reg[1] = b_reg[11];
    b_reg[3] = b_reg[ipiv[3] - 1], b_reg[11] = b_reg[ipiv[3] + 7];
    b_reg[ipiv[3] - 1] = t_reg[0], b_reg[ipiv[3] + 7] = t_reg[1];

    t_reg[0] = b_reg[4], t_reg[1] = b_reg[12];
    b_reg[4] = b_reg[ipiv[4] - 1], b_reg[12] = b_reg[ipiv[4] + 7];
    b_reg[ipiv[4] - 1] = t_reg[0], b_reg[ipiv[4] + 7] = t_reg[1];

    t_reg[0] = b_reg[5], t_reg[1] = b_reg[13];
    b_reg[5] = b_reg[ipiv[5] - 1], b_reg[13] = b_reg[ipiv[5] + 7];
    b_reg[ipiv[5] - 1] = t_reg[0], b_reg[ipiv[5] + 7] = t_reg[1];

    t_reg[0] = b_reg[6], t_reg[1] = b_reg[14];
    b_reg[6] = b_reg[ipiv[6] - 1], b_reg[14] = b_reg[ipiv[6] + 7];
    b_reg[ipiv[6] - 1] = t_reg[0], b_reg[ipiv[6] + 7] = t_reg[1];
    // REGION - TRSM row 1 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    b_reg[9] = _mm256_sub_pd(b_reg[9], t_reg[1]);
    // REGION - TRSM row 2 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);
    b_reg[10] = _mm256_sub_pd(b_reg[10], t_reg[1]);
    // REGION - TRSM row 3 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    b_reg[3] = _mm256_sub_pd(b_reg[3], t_reg[0]);
    b_reg[11] = _mm256_sub_pd(b_reg[11], t_reg[1]);
    // REGION - TRSM row 4 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    b_reg[4] = _mm256_sub_pd(b_reg[4], t_reg[0]);
    b_reg[12] = _mm256_sub_pd(b_reg[12], t_reg[1]);
    // REGION - TRSM row 5 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    b_reg[5] = _mm256_sub_pd(b_reg[5], t_reg[0]);
    b_reg[13] = _mm256_sub_pd(b_reg[13], t_reg[1]);
    // REGION - TRSM row 6 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    b_reg[6] = _mm256_sub_pd(b_reg[6], t_reg[0]);
    b_reg[14] = _mm256_sub_pd(b_reg[14], t_reg[1]);
    // REGION - TRSM row 6 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 6 + (6 * (*lda))));
    b_reg[6] = _mm256_div_pd(b_reg[6], t_reg[2]);
    b_reg[14] = _mm256_div_pd(b_reg[14], t_reg[2]);
    // REGION - TRSM row 5 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (6 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[6], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[14], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (5 * (*lda))));
    b_reg[5] = _mm256_sub_pd(b_reg[5], t_reg[0]);
    b_reg[5] = _mm256_div_pd(b_reg[5], t_reg[2]);
    b_reg[13] = _mm256_sub_pd(b_reg[13], t_reg[1]);
    b_reg[13] = _mm256_div_pd(b_reg[13], t_reg[2]);
    // REGION - TRSM row 4 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (6 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[6], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[14], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (4 * (*lda))));
    b_reg[4] = _mm256_sub_pd(b_reg[4], t_reg[0]);
    b_reg[4] = _mm256_div_pd(b_reg[4], t_reg[2]);
    b_reg[12] = _mm256_sub_pd(b_reg[12], t_reg[1]);
    b_reg[12] = _mm256_div_pd(b_reg[12], t_reg[2]);
    // REGION - TRSM row 3 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (6 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[6], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[14], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (3 * (*lda))));
    b_reg[3] = _mm256_sub_pd(b_reg[3], t_reg[0]);
    b_reg[3] = _mm256_div_pd(b_reg[3], t_reg[2]);
    b_reg[11] = _mm256_sub_pd(b_reg[11], t_reg[1]);
    b_reg[11] = _mm256_div_pd(b_reg[11], t_reg[2]);
    // REGION - TRSM row 2 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (6 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[6], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[14], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (2 * (*lda))));
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);
    b_reg[2] = _mm256_div_pd(b_reg[2], t_reg[2]);
    b_reg[10] = _mm256_sub_pd(b_reg[10], t_reg[1]);
    b_reg[10] = _mm256_div_pd(b_reg[10], t_reg[2]);
    // REGION - TRSM row 1 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (6 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[6], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[14], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (1 * (*lda))));
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    b_reg[1] = _mm256_div_pd(b_reg[1], t_reg[2]);
    b_reg[9] = _mm256_sub_pd(b_reg[9], t_reg[1]);
    b_reg[9] = _mm256_div_pd(b_reg[9], t_reg[2]);
    // REGION - TRSM row 0 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (6 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[6], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[14], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (5 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[5], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[13], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a));
    b_reg[0] = _mm256_sub_pd(b_reg[0], t_reg[0]);
    b_reg[0] = _mm256_div_pd(b_reg[0], t_reg[2]);
    b_reg[8] = _mm256_sub_pd(b_reg[8], t_reg[1]);
    b_reg[8] = _mm256_div_pd(b_reg[8], t_reg[2]);

    TRANSPOSE_8x8() for (i = 0; i < (*nrhs); ++i)
    {
        j = (int)(i / 4);
        j *= 4;
        _mm256_maskstore_pd((b + (i * (*ldb)) + 0), mask1, b_reg[i + 0 + j]);
        _mm256_maskstore_pd((b + (i * (*ldb)) + 4), mask2, b_reg[i + 4 + j]);
    }
}

static void n_6(integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *ipiv)
{
    __m256d b_reg[16];
    __m256i mask1, mask2;
    __m256d t_reg[3];
    mask1 = mask_reg[4];
    mask2 = mask_reg[2];
    int i, j;
    for (i = 0; i < (*nrhs); ++i)
    {
        j = (int)(i / 4);
        j *= 4;
        b_reg[i + 0 + j] = _mm256_maskload_pd((void const *)(b + (i * (*ldb)) + 0), mask1);
        b_reg[i + 4 + j] = _mm256_maskload_pd((void const *)(b + (i * (*ldb)) + 4), mask2);
    }
    TRANSPOSE_8x8()
    /* After tranpose, B matrix is stored in row major format in the b_reg registers*/
    /* so in order to swap row, we only need to swap registers                      */
    /* Swap Row [n] with Row [Ipiv[n]]                                              */
    t_reg[0] = b_reg[0],
    t_reg[1] = b_reg[8];
    b_reg[0] = b_reg[ipiv[0] - 1], b_reg[8] = b_reg[ipiv[0] + 7];
    b_reg[ipiv[0] - 1] = t_reg[0], b_reg[ipiv[0] + 7] = t_reg[1];

    t_reg[0] = b_reg[1], t_reg[1] = b_reg[9];
    b_reg[1] = b_reg[ipiv[1] - 1], b_reg[9] = b_reg[ipiv[1] + 7];
    b_reg[ipiv[1] - 1] = t_reg[0], b_reg[ipiv[1] + 7] = t_reg[1];

    t_reg[0] = b_reg[2], t_reg[1] = b_reg[10];
    b_reg[2] = b_reg[ipiv[2] - 1], b_reg[10] = b_reg[ipiv[2] + 7];
    b_reg[ipiv[2] - 1] = t_reg[0], b_reg[ipiv[2] + 7] = t_reg[1];

    t_reg[0] = b_reg[3], t_reg[1] = b_reg[11];
    b_reg[3] = b_reg[ipiv[3] - 1], b_reg[11] = b_reg[ipiv[3] + 7];
    b_reg[ipiv[3] - 1] = t_reg[0], b_reg[ipiv[3] + 7] = t_reg[1];

    t_reg[0] = b_reg[4], t_reg[1] = b_reg[12];
    b_reg[4] = b_reg[ipiv[4] - 1], b_reg[12] = b_reg[ipiv[4] + 7];
    b_reg[ipiv[4] - 1] = t_reg[0], b_reg[ipiv[4] + 7] = t_reg[1];

    t_reg[0] = b_reg[5], t_reg[1] = b_reg[13];
    b_reg[5] = b_reg[ipiv[5] - 1], b_reg[13] = b_reg[ipiv[5] + 7];
    b_reg[ipiv[5] - 1] = t_reg[0], b_reg[ipiv[5] + 7] = t_reg[1];
    // REGION - TRSM row 1 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    b_reg[9] = _mm256_sub_pd(b_reg[9], t_reg[1]);
    // REGION - TRSM row 2 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);
    b_reg[10] = _mm256_sub_pd(b_reg[10], t_reg[1]);
    // REGION - TRSM row 3 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    b_reg[3] = _mm256_sub_pd(b_reg[3], t_reg[0]);
    b_reg[11] = _mm256_sub_pd(b_reg[11], t_reg[1]);
    // REGION - TRSM row 4 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    b_reg[4] = _mm256_sub_pd(b_reg[4], t_reg[0]);
    b_reg[12] = _mm256_sub_pd(b_reg[12], t_reg[1]);
    // REGION - TRSM row 5 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    b_reg[5] = _mm256_sub_pd(b_reg[5], t_reg[0]);
    b_reg[13] = _mm256_sub_pd(b_reg[13], t_reg[1]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 5 + (5 * (*lda))));
    b_reg[5] = _mm256_div_pd(b_reg[5], t_reg[2]);
    b_reg[13] = _mm256_div_pd(b_reg[13], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (5 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[5], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[13], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (4 * (*lda))));
    b_reg[4] = _mm256_sub_pd(b_reg[4], t_reg[0]);
    b_reg[4] = _mm256_div_pd(b_reg[4], t_reg[2]);
    b_reg[12] = _mm256_sub_pd(b_reg[12], t_reg[1]);
    b_reg[12] = _mm256_div_pd(b_reg[12], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (5 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[5], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[13], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (3 * (*lda))));
    b_reg[3] = _mm256_sub_pd(b_reg[3], t_reg[0]);
    b_reg[3] = _mm256_div_pd(b_reg[3], t_reg[2]);
    b_reg[11] = _mm256_sub_pd(b_reg[11], t_reg[1]);
    b_reg[11] = _mm256_div_pd(b_reg[11], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (5 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[5], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[13], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (2 * (*lda))));
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);
    b_reg[2] = _mm256_div_pd(b_reg[2], t_reg[2]);
    b_reg[10] = _mm256_sub_pd(b_reg[10], t_reg[1]);
    b_reg[10] = _mm256_div_pd(b_reg[10], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (5 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[5], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[13], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (1 * (*lda))));
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    b_reg[1] = _mm256_div_pd(b_reg[1], t_reg[2]);
    b_reg[9] = _mm256_sub_pd(b_reg[9], t_reg[1]);
    b_reg[9] = _mm256_div_pd(b_reg[9], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (5 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[5], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[13], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (4 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[4], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[12], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a));
    b_reg[0] = _mm256_sub_pd(b_reg[0], t_reg[0]);
    b_reg[0] = _mm256_div_pd(b_reg[0], t_reg[2]);
    b_reg[8] = _mm256_sub_pd(b_reg[8], t_reg[1]);
    b_reg[8] = _mm256_div_pd(b_reg[8], t_reg[2]);

    TRANSPOSE_8x8() for (i = 0; i < (*nrhs); ++i)
    {
        j = (int)(i / 4);
        j *= 4;
        _mm256_maskstore_pd((b + (i * (*ldb)) + 0), mask1, b_reg[i + 0 + j]);
        _mm256_maskstore_pd((b + (i * (*ldb)) + 4), mask2, b_reg[i + 4 + j]);
    }
}

static void n_5(integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *ipiv)
{
    __m256d b_reg[16];
    __m256i mask1, mask2;
    __m256d t_reg[3];
    mask1 = mask_reg[4];
    mask2 = mask_reg[1];
    int i, j;
    for (i = 0; i < (*nrhs); ++i)
    {
        j = (int)(i / 4);
        j *= 4;
        b_reg[i + 0 + j] = _mm256_maskload_pd((void const *)(b + (i * (*ldb)) + 0), mask1);
        b_reg[i + 4 + j] = _mm256_maskload_pd((void const *)(b + (i * (*ldb)) + 4), mask2);
    }

    b_reg[0] = _mm256_maskload_pd((void const *)(b + (0 * (*ldb)) + 0), mask1);
    b_reg[4] = _mm256_maskload_pd((void const *)(b + (0 * (*ldb)) + 4), mask2);
    b_reg[1] = _mm256_maskload_pd((void const *)(b + (1 * (*ldb)) + 0), mask1);
    b_reg[5] = _mm256_maskload_pd((void const *)(b + (1 * (*ldb)) + 4), mask2);
    b_reg[2] = _mm256_maskload_pd((void const *)(b + (2 * (*ldb)) + 0), mask1);
    b_reg[6] = _mm256_maskload_pd((void const *)(b + (2 * (*ldb)) + 4), mask2);
    b_reg[3] = _mm256_maskload_pd((void const *)(b + (3 * (*ldb)) + 0), mask1);
    b_reg[7] = _mm256_maskload_pd((void const *)(b + (3 * (*ldb)) + 4), mask2);
    b_reg[8] = _mm256_maskload_pd((void const *)(b + (4 * (*ldb)) + 0), mask1);
    b_reg[12] = _mm256_maskload_pd((void const *)(b + (4 * (*ldb)) + 4), mask2);
    TRANSPOSE_8x8()
    /* After tranpose, B matrix is stored in row major format in the b_reg registers*/
    /* so in order to swap row, we only need to swap registers                      */
    /* Swap Row [n] with Row [Ipiv[n]]                                              */
    t_reg[0] = b_reg[0],
    t_reg[1] = b_reg[8];
    b_reg[0] = b_reg[ipiv[0] - 1], b_reg[8] = b_reg[ipiv[0] + 7];
    b_reg[ipiv[0] - 1] = t_reg[0], b_reg[ipiv[0] + 7] = t_reg[1];

    t_reg[0] = b_reg[1], t_reg[1] = b_reg[9];
    b_reg[1] = b_reg[ipiv[1] - 1], b_reg[9] = b_reg[ipiv[1] + 7];
    b_reg[ipiv[1] - 1] = t_reg[0], b_reg[ipiv[1] + 7] = t_reg[1];

    t_reg[0] = b_reg[2], t_reg[1] = b_reg[10];
    b_reg[2] = b_reg[ipiv[2] - 1], b_reg[10] = b_reg[ipiv[2] + 7];
    b_reg[ipiv[2] - 1] = t_reg[0], b_reg[ipiv[2] + 7] = t_reg[1];

    t_reg[0] = b_reg[3], t_reg[1] = b_reg[11];
    b_reg[3] = b_reg[ipiv[3] - 1], b_reg[11] = b_reg[ipiv[3] + 7];
    b_reg[ipiv[3] - 1] = t_reg[0], b_reg[ipiv[3] + 7] = t_reg[1];

    t_reg[0] = b_reg[4], t_reg[1] = b_reg[12];
    b_reg[4] = b_reg[ipiv[4] - 1], b_reg[12] = b_reg[ipiv[4] + 7];
    b_reg[ipiv[4] - 1] = t_reg[0], b_reg[ipiv[4] + 7] = t_reg[1];
    // REGION - TRSM row 1 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    b_reg[9] = _mm256_sub_pd(b_reg[9], t_reg[1]);
    // REGION - TRSM row 2 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);
    b_reg[10] = _mm256_sub_pd(b_reg[10], t_reg[1]);
    // REGION - TRSM row 3 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    b_reg[3] = _mm256_sub_pd(b_reg[3], t_reg[0]);
    b_reg[11] = _mm256_sub_pd(b_reg[11], t_reg[1]);
    // REGION - TRSM row 4 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[8], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    b_reg[4] = _mm256_sub_pd(b_reg[4], t_reg[0]);
    b_reg[12] = _mm256_sub_pd(b_reg[12], t_reg[1]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 4 + (4 * (*lda))));
    b_reg[4] = _mm256_div_pd(b_reg[4], t_reg[2]);
    b_reg[12] = _mm256_div_pd(b_reg[12], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (4 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[4], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[12], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (3 * (*lda))));
    b_reg[3] = _mm256_sub_pd(b_reg[3], t_reg[0]);
    b_reg[3] = _mm256_div_pd(b_reg[3], t_reg[2]);
    b_reg[11] = _mm256_sub_pd(b_reg[11], t_reg[1]);
    b_reg[11] = _mm256_div_pd(b_reg[11], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (4 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[4], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[12], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (2 * (*lda))));
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);
    b_reg[2] = _mm256_div_pd(b_reg[2], t_reg[2]);
    b_reg[10] = _mm256_sub_pd(b_reg[10], t_reg[1]);
    b_reg[10] = _mm256_div_pd(b_reg[10], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (4 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[4], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[12], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (1 * (*lda))));
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    b_reg[1] = _mm256_div_pd(b_reg[1], t_reg[2]);
    b_reg[9] = _mm256_sub_pd(b_reg[9], t_reg[1]);
    b_reg[9] = _mm256_div_pd(b_reg[9], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (4 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[4], t_reg[2]);
    t_reg[1] = _mm256_mul_pd(b_reg[12], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (3 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[3], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[11], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[10], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[1] = _mm256_fmadd_pd(b_reg[9], t_reg[2], t_reg[1]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a));
    b_reg[0] = _mm256_sub_pd(b_reg[0], t_reg[0]);
    b_reg[0] = _mm256_div_pd(b_reg[0], t_reg[2]);
    b_reg[8] = _mm256_sub_pd(b_reg[8], t_reg[1]);
    b_reg[8] = _mm256_div_pd(b_reg[8], t_reg[2]);

    TRANSPOSE_8x8() for (i = 0; i < (*nrhs); ++i)
    {
        j = (int)(i / 4);
        j *= 4;
        _mm256_maskstore_pd((b + (i * (*ldb)) + 0), mask1, b_reg[i + 0 + j]);
        _mm256_maskstore_pd((b + (i * (*ldb)) + 4), mask2, b_reg[i + 4 + j]);
    }
}

static void n_4(integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *ipiv)
{
    int i;
    __m256d b_reg[4];
    __m256i mask1;
    __m256d t_reg[4];
    mask1 = mask_reg[4];
    for (i = 0; i < (*nrhs); ++i)
    {
        b_reg[i] = _mm256_maskload_pd((void const *)(b + (i * (*ldb))), mask1);
    }
    TRANSPOSE_4x4(0, 1, 2, 3)
    /* After tranpose, B matrix is stored in row major format in the b_reg registers*/
    /* so in order to swap row, we only need to swap registers                      */
    /* Swap Row [n] with Row [Ipiv[n]]                                              */
    t_reg[0] = b_reg[0];
    b_reg[0] = b_reg[ipiv[0] - 1];
    b_reg[ipiv[0] - 1] = t_reg[0];
    t_reg[1] = b_reg[1];
    b_reg[1] = b_reg[ipiv[1] - 1];
    b_reg[ipiv[1] - 1] = t_reg[1];
    t_reg[2] = b_reg[2];
    b_reg[2] = b_reg[ipiv[2] - 1];
    b_reg[ipiv[2] - 1] = t_reg[2];
    t_reg[3] = b_reg[3];
    b_reg[3] = b_reg[ipiv[3] - 1];
    b_reg[ipiv[3] - 1] = t_reg[3];
    // REGION - TRSM row 1 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    // REGION - TRSM row 2 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);
    // REGION - TRSM row 3 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 3 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    b_reg[3] = _mm256_sub_pd(b_reg[3], t_reg[0]);

    t_reg[3] = _mm256_broadcast_sd((double const *)(a + 3 + (3 * (*lda))));
    b_reg[3] = _mm256_div_pd(b_reg[3], t_reg[3]);

    t_reg[3] = _mm256_broadcast_sd((double const *)(a + 2 + (2 * (*lda))));
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (3 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[3], t_reg[2]);
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);
    b_reg[2] = _mm256_div_pd(b_reg[2], t_reg[3]);

    t_reg[3] = _mm256_broadcast_sd((double const *)(a + 1 + (1 * (*lda))));
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (3 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[3], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    b_reg[1] = _mm256_div_pd(b_reg[1], t_reg[3]);

    t_reg[3] = _mm256_broadcast_sd((double const *)(a + 0 + (0 * (*lda))));
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 0 + (3 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[3], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 0 + (2 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[2], t_reg[2], t_reg[0]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 0 + (1 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    b_reg[0] = _mm256_sub_pd(b_reg[0], t_reg[0]);
    b_reg[0] = _mm256_div_pd(b_reg[0], t_reg[3]);
    TRANSPOSE_4x4(0, 1, 2, 3) for (i = 0; i < (*nrhs); ++i)
    {
        _mm256_maskstore_pd((b + (i * (*ldb))), mask1, b_reg[i]);
    }
}

static void n_3(integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *ipiv)
{
    int i;
    __m256d b_reg[4];
    __m256i mask1;
    __m256d t_reg[4];
    mask1 = mask_reg[3];
    for (i = 0; i < (*nrhs); ++i)
    {
        b_reg[i] = _mm256_maskload_pd((void const *)(b + (i * (*ldb))), mask1);
    }
    TRANSPOSE_4x4(0, 1, 2, 3)
    /* After tranpose, B matrix is stored in row major format in the b_reg registers*/
    /* so in order to swap row, we only need to swap registers                      */
    /* Swap Row [n] with Row [Ipiv[n]]                                              */
    t_reg[0] = b_reg[0];
    b_reg[0] = b_reg[ipiv[0] - 1];
    b_reg[ipiv[0] - 1] = t_reg[0];
    t_reg[1] = b_reg[1];
    b_reg[1] = b_reg[ipiv[1] - 1];
    b_reg[ipiv[1] - 1] = t_reg[1];
    t_reg[2] = b_reg[2];
    b_reg[2] = b_reg[ipiv[2] - 1];
    b_reg[ipiv[2] - 1] = t_reg[2];
    // REGION - TRSM row 1 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    // REGION - TRSM row 2 computation
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2));
    t_reg[0] = _mm256_mul_pd(b_reg[0], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (*lda)));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    b_reg[2] = _mm256_sub_pd(b_reg[2], t_reg[0]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 2 + (2 * (*lda))));
    b_reg[2] = _mm256_div_pd(b_reg[2], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (2 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[2], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 1 + (1 * (*lda))));
    b_reg[1] = _mm256_sub_pd(b_reg[1], t_reg[0]);
    b_reg[1] = _mm256_div_pd(b_reg[1], t_reg[2]);

    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 0 + (2 * (*lda))));
    t_reg[0] = _mm256_mul_pd(b_reg[2], t_reg[2]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 0 + (1 * (*lda))));
    t_reg[0] = _mm256_fmadd_pd(b_reg[1], t_reg[2], t_reg[0]);
    t_reg[2] = _mm256_broadcast_sd((double const *)(a + 0 + (0 * (*lda))));
    b_reg[0] = _mm256_sub_pd(b_reg[0], t_reg[0]);
    b_reg[0] = _mm256_div_pd(b_reg[0], t_reg[2]);
    TRANSPOSE_4x4(0, 1, 2, 3) for (i = 0; i < (*nrhs); ++i)
    {
        _mm256_maskstore_pd((b + (i * (*ldb))), mask1, b_reg[i]);
    }
}

/* Small DGETRS path (NOTRANS) should only be used for size between 3 to 8 and NRHS <= N */
int fla_dgetrs_small_trsm_ll_avx2(char *trans, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info)
{
    // Array of function pointers
    void (*fp[6]) (integer *, integer *, doublereal *, integer *,doublereal *, integer *, integer *) = {n_3, n_4, n_5, n_6, n_7, n_8};
    fp[*n - 3](n, nrhs, a, lda, b, ldb, ipiv);
    return 0;
}
#endif
