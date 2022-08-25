/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

/* Test API function declaration */
void fla_test_steqr(test_params_t *params);
void fla_test_stevd(test_params_t *params);
void fla_test_geev(test_params_t *params);
void fla_test_geevx(test_params_t *params);
void fla_test_gesdd(test_params_t *params);
void fla_test_geqrf(test_params_t *params);
void fla_test_gerqf(test_params_t *params);
void fla_test_gerq2(test_params_t *params);
void fla_test_geqp3(test_params_t *params);
void fla_test_potrf(test_params_t *params);
void fla_test_getrf(test_params_t *params);
void fla_test_getri(test_params_t *params);
void fla_test_orgqr(test_params_t *params);
void fla_test_potrs(test_params_t *params);
void fla_test_getrs(test_params_t* params);
void fla_test_syevd(test_params_t *params);
void fla_test_gesvd(test_params_t *params);
void fla_test_ggevx(test_params_t *params);
void fla_test_gesv(test_params_t* params);
void fla_test_ggev(test_params_t *params);
void fla_test_stedc(test_params_t *params);

/* Add test api function call entry below */
OPERATIONS API_test_functions[] =
{
    {"orgqr"               , fla_test_orgqr},
    {"potrs"               , fla_test_potrs},
    {"geev"                , fla_test_geev},
    {"geevx"               , fla_test_geevx},
    {"gesdd"               , fla_test_gesdd},
    {"potrf"               , fla_test_potrf},
    {"geqrf"               , fla_test_geqrf},
    {"gerqf"               , fla_test_gerqf},
    {"gerq2"               , fla_test_gerq2},
    {"geqp3"               , fla_test_geqp3},
    {"getrf"               , fla_test_getrf},
    {"getri"               , fla_test_getri},
    {"getrs"               , fla_test_getrs},
    {"syevd"               , fla_test_syevd},
    {"gesvd"               , fla_test_gesvd},
    {"ggevx"               , fla_test_ggevx},
    {"gesv"                , fla_test_gesv},
    {"ggev"                , fla_test_ggev},
    {"steqr"               , fla_test_steqr},
    {"stevd"               , fla_test_stevd},
    {"stedc"               , fla_test_stedc}
};
