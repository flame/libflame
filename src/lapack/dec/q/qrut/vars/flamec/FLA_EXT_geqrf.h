/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    May 09, 2021
*/

FLA_Error FLA_EXT_sgeqrf( integer m_A, integer n_A,
                          float* buff_A, integer cs_A,
                          float* buff_t,
                          float* buff_w,
                          integer* lwork,
                          integer* info );
FLA_Error FLA_EXT_dgeqrf( integer m_A, integer n_A,
                          double* buff_A, integer cs_A,
                          double* buff_t,
                          double* buff_w,
                          integer* lwork,
                          integer* info );
