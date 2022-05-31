/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

/* Update test api count */
#define test_api_count 11


/* Test API function declaration */
void fla_test_geevx(test_params_t *params);
void fla_test_gesdd(test_params_t *params);
void fla_test_geqrf(test_params_t *params);
void fla_test_gerqf(test_params_t *params);
void fla_test_gerq2(test_params_t *params);
void fla_test_potrf(test_params_t *params);
void fla_test_getrf(test_params_t *params);
void fla_test_getri(test_params_t *params);
void fla_test_orgqr(test_params_t *params);
void fla_test_potrs(test_params_t *params);
void fla_test_getrs(test_params_t* params);

/* Add test api function call entry below */
OPERATIONS API_test_functions[] =
{
    {"orgqr"               , fla_test_orgqr},
    {"potrs"               , fla_test_potrs},
    {"geevx"               , fla_test_geevx},
    {"gesdd"               , fla_test_gesdd},
    {"potrf"               , fla_test_potrf},
    {"geqrf"               , fla_test_geqrf},
    {"gerqf"               , fla_test_gerqf},
    {"gerq2"               , fla_test_gerq2},
    {"getrf"               , fla_test_getrf},
    {"getri"               , fla_test_getri},
    {"getrs"               , fla_test_getrs}
};
