/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

/* Test API function declaration */
void fla_test_steqr(integer argc, char ** argv, test_params_t *params);
void fla_test_stevd(integer argc, char ** argv, test_params_t *params);
void fla_test_geev(integer argc, char ** argv, test_params_t *params);
void fla_test_geevx(integer argc, char ** argv, test_params_t *params);
void fla_test_gesdd(integer argc, char ** argv, test_params_t *params);
void fla_test_geqrf(integer argc, char ** argv, test_params_t *params);
void fla_test_gerqf(integer argc, char ** argv, test_params_t *params);
void fla_test_gerq2(integer argc, char ** argv, test_params_t *params);
void fla_test_geqp3(integer argc, char ** argv, test_params_t *params);
void fla_test_potrf(integer argc, char ** argv, test_params_t *params);
void fla_test_getrf(integer argc, char ** argv, test_params_t *params);
void fla_test_getri(integer argc, char ** argv, test_params_t *params);
void fla_test_orgqr(integer argc, char ** argv, test_params_t *params);
void fla_test_potrs(integer argc, char ** argv, test_params_t *params);
void fla_test_getrs(integer argc, char ** argv, test_params_t* params);
void fla_test_syevd(integer argc, char ** argv, test_params_t *params);
void fla_test_gesvd(integer argc, char ** argv, test_params_t *params);
void fla_test_ggevx(integer argc, char ** argv, test_params_t *params);
void fla_test_gesv(integer argc, char ** argv, test_params_t* params);
void fla_test_ggev(integer argc, char ** argv, test_params_t *params);
void fla_test_stedc(integer argc, char ** argv, test_params_t *params);

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
