/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_lartg_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, integer einfo, double* perf, double* t, double* residual);
void prepare_lartg_run(integer datatype, void *f, void *g, void *r, void *c, void *s, 
integer n_repeats, double* time_min_);
void invoke_lartg(integer datatype, void *f, void *g, void *c, void *s, void *r);

void fla_test_lartg(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Auxilary routines";
    char* front_str = "LARTG";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    integer i, num_types;
    integer datatype, n_repeats;
    double perf, time_min, residual;
    char stype, type_flag[4] = {0};
    char *endptr;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        num_types = params->aux_paramslist[0].num_data_types;
        n_repeats = params->aux_paramslist[0].num_repeats;

        if (n_repeats > 0)
        {
            /* Loop over the requested datatypes. */
            for ( i = 0; i < num_types; ++i )
            {
                datatype = params->datatype[i];
                stype    = params->datatype_char[i];

                /* Call the test code */
                fla_test_lartg_experiment(params, datatype,
                                          2, i_one,
                                          0,
                                          n_repeats, einfo,
                                          &perf, &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str,
                                     stype,
                                     RECT_INPUT,
                                     2, i_one,
                                     residual, params->aux_paramslist[0].aux_threshold,
                                     time_min, perf);
                tests_not_run = 0;
            }
        }
    }
    if (argc == 5)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[4]);
    }

    if (argc >= 4 && argc <= 5)
    {
        /* Test with parameters from commandline */
        /* Parse the arguments */
        num_types = strlen(argv[2]);
        
        n_repeats = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->aux_paramslist[0].aux_threshold = CLI_NORM_THRESH;

            for(i = 0; i < num_types; i++)
            {
                stype = argv[2][i];
                datatype = get_datatype(stype);

                /* Check for invalide dataype */
                if(datatype == INVALID_TYPE)
                {
                    invalid_dtype = 1;
                    continue;
                }

                /* Check for duplicate datatype presence */
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                /* Call the test code */
                fla_test_lartg_experiment(params, datatype,
                                          2, i_one,
                                          0,
                                          n_repeats, einfo,
                                          &perf, &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str,
                                      stype,
                                      RECT_INPUT,
                                      2, i_one,
                                      residual, params->aux_paramslist[0].aux_threshold,
                                      time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\n Illegal arguments for lartg \n");
        printf("./<EXE> lartg <precisions - sdcz> <repeats> [file] \n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
    if (g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
        g_ext_fptr = NULL;
    }

    return;
}

void fla_test_lartg_experiment(test_params_t *params,
    integer  datatype,
    integer  p_cur,
    integer  q_cur,
    integer pci,
    integer n_repeats,
    integer einfo,
    double* perf,
    double* t,
    double* residual)
{
    void  *s = NULL, *c = NULL;
    void *f = NULL, *g = NULL, *r = NULL;
    double time_min = 1e9;

    integer realtype;
    realtype = get_realtype(datatype);

    *residual = params->aux_paramslist[pci].aux_threshold;

    create_vector(realtype, &c, 1);
    create_vector(datatype, &s, 1);

    create_vector(datatype, &f, 1);
    create_vector(datatype, &g, 1);
    create_vector(datatype, &r, 1);
    
    if(g_ext_fptr != NULL)
    {
        init_vector_from_file(datatype, f, 1, 1, g_ext_fptr);
        init_vector_from_file(datatype, g, 1, 1, g_ext_fptr);
    }
    else
    {
        rand_vector(datatype, f, 1, 1);
        rand_vector(datatype, g, 1, 1);
    } 
    /* call to API */
    prepare_lartg_run(datatype, f, g, r, c, s, n_repeats, &time_min);

    /* execution time */
    *t = time_min;
    if(time_min == d_zero)
    {
        time_min = 1e-9;
        *t = time_min;
    }
    /* Compute the performance of the best experiment repeat */
    *perf = (double)(6.0) / time_min / FLOPS_PER_UNIT_PERF;
    
    /* output validation */
    validate_lartg(datatype, f, g, r, c, s, residual);

    /* Free up the buffers */
    free_vector(c);
    free_vector(s);
    free_vector(f);
    free_vector(g);
    free_vector(r);
}

void prepare_lartg_run(integer datatype,
    void *f,
    void *g,
    void *r,
    void *c,
    void *s,
    integer n_repeats,
    double* time_min_)
{
    integer i;
    double time_min = 1e9, exe_time;

    for (i = 0; i < n_repeats; ++i)
    {
        exe_time = fla_test_clock();

        /*  call  lartg API */
        invoke_lartg(datatype, f, g, c, s, r);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
}

void invoke_lartg(integer datatype, void *f, void *g, void *c, void *s, void *r)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_slartg(f, g, (float *) c, (float *)s, r);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dlartg(f, g, (double *) c, (double *)s, r);
            break;
        }
        case COMPLEX:
        {
           fla_lapack_clartg(f, g, (float *) c, ((scomplex *)s), r);
           break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zlartg(f, g, (double *) c, ((dcomplex *)s), r);
            break;
        }
    }
}