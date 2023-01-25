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
void fla_test_syev(integer argc, char ** argv, test_params_t *params);

#define LIN_ID 0
#define EIG_ID 1
#define SVD_ID 2

/* Add test api function call entry below */
OPERATIONS API_test_functions[] =
{
    {LIN_ID,    "orgqr"               , fla_test_orgqr},
    {LIN_ID,    "potrs"               , fla_test_potrs},
    {EIG_ID,    "geev"                , fla_test_geev},
    {EIG_ID,    "geevx"               , fla_test_geevx},
    {SVD_ID,    "gesdd"               , fla_test_gesdd},
    {LIN_ID,    "potrf"               , fla_test_potrf},
    {LIN_ID,    "geqrf"               , fla_test_geqrf},
    {LIN_ID,    "gerqf"               , fla_test_gerqf},
    {LIN_ID,    "gerq2"               , fla_test_gerq2},
    {LIN_ID,    "geqp3"               , fla_test_geqp3},
    {LIN_ID,    "getrf"               , fla_test_getrf},
    {LIN_ID,    "getri"               , fla_test_getri},
    {LIN_ID,    "getrs"               , fla_test_getrs},
    {EIG_ID,    "syevd"               , fla_test_syevd},
    {SVD_ID,    "gesvd"               , fla_test_gesvd},
    {EIG_ID,    "ggevx"               , fla_test_ggevx},
    {LIN_ID,    "gesv"                , fla_test_gesv},
    {EIG_ID,    "ggev"                , fla_test_ggev},
    {EIG_ID,    "steqr"               , fla_test_steqr},
    {EIG_ID,    "stevd"               , fla_test_stevd},
    {EIG_ID,    "stedc"               , fla_test_stedc},
    {EIG_ID,    "syev"                , fla_test_syev}
};

/* Add test API's group entry below */
char *API_test_group[] =
{
    "LIN",
    "EIG",
    "SVD"
};
