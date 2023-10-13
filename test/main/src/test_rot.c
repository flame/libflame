/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_rot_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer pci,
                                    integer n_repeats, integer einfo, double* perf, double* t, double* residual);
void prepare_rot_run(integer datatype, integer n_A, void* cx, integer incx, void* cy, integer incy, void *c, void *s, 
integer n_repeats, double* time_min_);
void invoke_rot(integer datatype, integer *n, void *cx, integer *incx, void *cy,integer *incy, void *c, void *s);
extern void invoke_lartg(integer datatype, void *f, void *g, void *c, void *s, void *r);

void fla_test_rot(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "Auxilary routines";
    char* front_str = "ROT";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, AUX, fla_test_rot_experiment);
        tests_not_run = 0;
    }
    if (argc == 8)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[7]);
    }

    if (argc >= 7 && argc <= 8)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        params->aux_paramslist[0].incx = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->aux_paramslist[0].incy = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_rot_experiment(params, datatype,
                                          N, N,
                                          0,
                                          n_repeats, einfo,
                                          &perf, &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str,
                                      stype,
                                      SQUARE_INPUT,
                                      N, N,
                                      residual, params->aux_paramslist[0].aux_threshold,
                                      time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for rot \n");
        printf("./<EXE> rot <precisions - sdcz>  <N> <incx>  <incy> <repeats> [file]\n");
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

void fla_test_rot_experiment(test_params_t *params,
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
    integer n, incx, incy;
    void *cx = NULL, *cy = NULL, *s = NULL, *c = NULL;;
    void *f = NULL, *g = NULL, *r = NULL;
    void *cx_test = NULL, *cy_test = NULL;
    double time_min = 1e9;

    integer realtype;
    realtype = get_realtype(datatype);

    *residual = params->aux_paramslist[pci].aux_threshold;
    incx = params->aux_paramslist[pci].incx;
    incy = params->aux_paramslist[pci].incy;
    /* Determine the dimensions*/    
    n = p_cur ;
    if (n <= 0)
    {
        *residual = DBL_MIN;
        return;
    }

    /* Create the vectors for the current operation*/
    create_vector(datatype, &cx, 1 + (n-1)*abs(incx));
    create_vector(datatype, &cy, 1 + (n-1)*abs(incy));

    create_vector(datatype, &cx_test, 1 + (n-1)*abs(incx));
    create_vector(datatype, &cy_test, 1 + (n-1)*abs(incy));

    create_vector(realtype, &c, 1);
    create_vector(datatype, &s, 1);

    create_vector(datatype, &f, 1);
    create_vector(datatype, &g, 1);
    create_vector(datatype, &r, 1);

    rand_vector(datatype, f, 1, 1);
    rand_vector(datatype, g, 1, 1);

    /*calling lartg api for getting c and s value for plane rotation*/
    invoke_lartg(datatype, f, g, c, s, r);

    if (g_ext_fptr != NULL)
    {
        /* Initialize input vectors with custom data */
        init_vector_from_file(datatype, cx, 1 + (n-1)*abs(incx), 1, g_ext_fptr);
        init_vector_from_file(datatype, cy, 1 + (n-1)*abs(incy), 1, g_ext_fptr);
    }
    else
    {
        /* Initialize input matrix with random numbers */
        rand_vector(datatype, cx, 1 + (n-1)*abs(incx), 1);
        rand_vector(datatype, cy, 1 + (n-1)*abs(incy), 1);
    }
    copy_vector(datatype, 1 + (n - 1)*abs(incx), cx, i_one, cx_test, i_one);
    copy_vector(datatype, 1 + (n - 1)*abs(incy), cy, i_one, cy_test, i_one);
    /* call to API */
    prepare_rot_run(datatype, n, cx, incx, cy, incy, c, s, n_repeats, &time_min);

    /* execution time */
    *t = time_min;
    if(time_min == d_zero)
    {
        time_min = 1e-9;
        *t = time_min;
    }
    /* Compute the performance of the best experiment repeat */
    /* 4*n */
    *perf = (double)(4.0 * n ) / time_min / FLOPS_PER_UNIT_PERF;
    /* output validation */
    validate_rot(datatype, n, cx, cx_test, incx, cy, cy_test, incy, c, s, residual);

    /* Free up the buffers */
    free_vector(cx);
    free_vector(cy);
    free_vector(cx_test);
    free_vector(cy_test);
    free_vector(c);
    free_vector(s);
    free_vector(f);
    free_vector(g);
    free_vector(r);
}


void prepare_rot_run(integer datatype,
    integer n_A,
    void* cx,
    integer incx,
    void* cy,
    integer incy,
    void *c,
    void *s,
    integer n_repeats,
    double* time_min_)
{
    integer i;
    void *cx_save = NULL, *cy_save = NULL;
    double time_min = 1e9, exe_time;

    create_vector(datatype, &cx_save, 1 + (n_A - 1)*abs(incx));
    create_vector(datatype, &cy_save, 1 + (n_A - 1)*abs(incy));

    for (i = 0; i < n_repeats; ++i)
    {
        /* Copy original input data */
        copy_vector(datatype, 1 + (n_A - 1)*abs(incx), cx, i_one, cx_save, i_one);
        copy_vector(datatype, 1 + (n_A - 1)*abs(incy), cy, i_one, cy_save, i_one);

        exe_time = fla_test_clock();

        /*  call  rot API */
        invoke_rot(datatype, &n_A, cx_save, &incx, cy_save, &incy, c, s);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    /*  Save the final result to A matrix*/
    copy_vector(datatype, 1 + (n_A - 1)*abs(incx), cx_save, i_one, cx, i_one);
    copy_vector(datatype, 1 + (n_A - 1)*abs(incy), cy_save, i_one, cy, i_one);

    free_vector(cx_save);
    free_vector(cy_save);
}


/*
 *  rot_API calls LAPACK interface
 *  */
void invoke_rot(integer datatype, integer *n, void *cx, integer *incx, void *cy,integer *incy, void *c, void *s)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_srot(n, cx, incx, cy, incy,((float *)c), ((float *)s));
            break;
        }

        case DOUBLE:
        {
            fla_lapack_drot(n, cx, incx, cy, incy,((double *)c), ((double *)s));
            break;
        }
        case COMPLEX:
        {
            fla_lapack_crot(n, cx, incx, cy, incy,((float *)c), ((scomplex *)s));
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zrot(n, cx, incx, cy, incy,((double *)c), ((dcomplex *)s));
            break;
        }
    }
}

