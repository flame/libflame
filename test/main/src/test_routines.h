/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

// Update test api count
#define test_api_count 4


// Test API function declaration
void fla_test_gesvd(test_params_t params);
void fla_test_geqp3(test_params_t params);
void fla_test_gerqf(test_params_t params);
void fla_test_gerq2(test_params_t params);


// Add test api function call entry below
OPERATIONS API_test_functions[] =
{
	{"gesvd"               , fla_test_gesvd},
	{"geqp3"               , fla_test_geqp3},
	{"gerqf"               , fla_test_gerqf},
	{"gerq2"               , fla_test_gerq2}
};