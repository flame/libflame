/*
    Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

extern double perf;
extern double time_min;

typedef enum
{
    FLA_LARF = 0,
    FLA_LARF1F = 1,
    FLA_LARF1L = 2
} fla_larf_variant_t;

/* Local prototypes */
void fla_test_larf_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo);
void prepare_larf_run(fla_larf_variant_t larf_variant, integer datatype, char side, integer m,
                      integer n, void *v, integer incv, void *tau, void *c__, integer ldc__,
                      void *c__out, integer ldc__out, void *work, integer interfacetype,
                      test_params_t *params);

void invoke_larf_common(fla_larf_variant_t larf_variant, integer datatype, char *side, integer *m,
                        integer *n, void *v, integer *incv, void *tau, void *c__, integer *ldc,
                        void *work);
void invoke_cpp_larf_common(fla_larf_variant_t larf_variant, integer datatype, char *side,
                            integer *m, integer *n, void *v, integer *incv, void *tau, void *c__,
                            integer *ldc, void *work);
void invoke_larfg(integer datatype, integer *n, void *alpha, void *x, integer *incx, void *tau);
void fla_test_larf_common(integer argc, char **argv, test_params_t *params, char *front_str);

void fla_test_larf(integer argc, char **argv, test_params_t *params)
{
    fla_test_larf_common(argc, argv, params, "LARF");
}

void fla_test_larf1f(integer argc, char **argv, test_params_t *params)
{
    fla_test_larf_common(argc, argv, params, "LARF1F");
}

void fla_test_larf1l(integer argc, char **argv, test_params_t *params)
{
    fla_test_larf_common(argc, argv, params, "LARF1L");
}

void fla_test_larf_common(integer argc, char **argv, test_params_t *params, char *front_str)
{
    char *op_str = "Auxilary routines";

    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, AUX, fla_test_larf_experiment);
        tests_not_run = 0;
    }
    if(argc == 10)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[9]);
    }

    if(argc >= 9 && argc <= 10)
    {
        /* Test with parameters from commandline */
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->aux_paramslist[0].side = argv[3][0];
        M = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->aux_paramslist[0].incv = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->aux_paramslist[0].ldc = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

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

                // Skip test for datatype and API combinations not supported.
                if(fla_skip_test(front_str, stype))
                {
                    printf("\nDatatype %c is not supported for %s, skipping the test.\n", stype,
                           front_str);
                    continue;
                }

                /* Check for duplicate datatype presence */
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                /* Call the test code */
                fla_test_larf_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for %s\n", front_str);
        printf("./<EXE> %s <precisions - sdcz>  <SIDE> <M> <N> <INCV> <LDC> <repeats>\n",
               front_str);
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
    if(g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
        g_ext_fptr = NULL;
    }

    return;
}

void fla_test_larf_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo)
{
    integer m, n;
    void *tau = NULL;
    integer v_length;
    void *work = NULL;
    void *v = NULL;
    void *v_tmp = NULL;
    void *c__ = NULL;
    void *c__out = NULL;
    double residual, err_thresh;
    void *filename = NULL;

    char side = params->aux_paramslist[pci].side;
    integer incv = params->aux_paramslist[pci].incv;
    integer ldc = params->aux_paramslist[pci].ldc;
    integer interfacetype = params->interfacetype;

    /* Choose the variant based on tst_api */
    fla_larf_variant_t larf_variant = FLA_LARF;

    if(same_string(tst_api, "LARF1F"))
    {
        larf_variant = FLA_LARF1F;
    }
    else if(same_string(tst_api, "LARF1L"))
    {
        larf_variant = FLA_LARF1L;
    }

    m = p_cur;
    n = q_cur;
    err_thresh = params->aux_paramslist[0].aux_threshold;

    if(m == 0 || n == 0)
        return;

    integer incv_abs = fla_i_abs(&incv);
    integer v_num_elements;
    integer work_num_elements;

    if(same_char(side, 'L'))
    {
        v_num_elements = m;
        work_num_elements = n;
    }
    else
    {
        v_num_elements = n;
        work_num_elements = m;
    }

    v_length = 1 + (v_num_elements - 1) * incv_abs;
    create_vector(datatype, &work, work_num_elements);
    create_vector(datatype, &tau, 1);
    create_vector(datatype, &v, v_length);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &c__, ldc);

    if(!FLA_BRT_VERIFICATION_RUN)
    {
        create_vector(datatype, &v_tmp, v_length);
        rand_vector(datatype, v_num_elements, v_tmp, incv_abs, d_zero, d_zero, 'R');

        /* To test correct handing of operations reduction due to zeroes in
         * v and C, we induce zeroes in v and C. A random length in 0 to 0.5 * length
         * is chosen and those many zeroes are induced in v and C.
         *
         * Do not set zeroes in benchmark mode or test-mode=perf as
         * it will not reflect the actual performance of the routine.
         */
        integer set_zero_flag = (params->benchmark_mode || params->test_mode == FLA_TEST_MODE_PERF
                                 || params->test_mode == FLA_TEST_MODE_RANDOM_PERF)
                                    ? 0
                                    : 1;
        integer num_zeroes_v = set_zero_flag ? (rand() % ((v_num_elements >> 1) + 1)) : 0;
        if(larf_variant != FLA_LARF1L)
        {
            /* put zeroes in the end of the vector */
            reset_vector(
                datatype,
                get_v_ptr(datatype, v_tmp, v_num_elements, v_num_elements - num_zeroes_v, incv_abs),
                num_zeroes_v, incv_abs);
        }
        else
        {
            /* put zeroes in the starting of the vector */
            reset_vector(datatype, v_tmp, num_zeroes_v, incv_abs);
        }

        /* Input generation (v_tmp and tau) for larf from larfg
        Increment of v_tmp for larfg must be positive. Hence calling larfg with incv_abs
        Increment of v for larf could be positive or negative. Hence copying
        from v_tmp using incv(which could be positive or negative)
        */
        void *alpha = larf_variant != FLA_LARF1L ? v_tmp
                                                 : get_v_ptr(datatype, v_tmp, v_num_elements,
                                                             v_num_elements - 1, incv_abs);
        void *x_for_larfg = larf_variant != FLA_LARF1L
                                ? get_v_ptr(datatype, v_tmp, v_num_elements, 1, incv_abs)
                                : v_tmp;
        invoke_larfg(datatype, &v_num_elements, alpha, x_for_larfg, &incv_abs, tau);
        if(larf_variant == FLA_LARF)
            assign_value(datatype, v_tmp, 1, 0);

        copy_vector(datatype, v_num_elements, v_tmp, incv_abs, v, incv);

        init_matrix(datatype, c__, m, n, ldc, g_ext_fptr, params->imatrix_char);

        integer num_zeroes_c = set_zero_flag ? (rand() % ((work_num_elements >> 1) + 1)) : 0;

        if(same_char(side, 'L'))
        {
            /* put zeroes in the last num_zeroes_c columns */
            reset_matrix(datatype, m, num_zeroes_c,
                         get_m_ptr(datatype, c__, 0, n - num_zeroes_c, ldc), ldc);
        }
        else
        {
            /* put zeroes in the last num_zeroes_c rows */
            reset_matrix(datatype, num_zeroes_c, n,
                         get_m_ptr(datatype, c__, m - num_zeroes_c, 0, ldc), ldc);
        }

        free_vector(v_tmp);
    }
    FLA_BRT_PROCESS_THREE_INPUT(datatype, m, n, c__, ldc, datatype, v_length, 1, v, v_length,
                                datatype, 1, 1, tau, 1, "cdddd", side, m, n, incv, ldc);

    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        scale_matrix_underflow_overflow_larf(datatype, m, n, c__, ldc, params->imatrix_char);
    }

    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &c__out, ldc);

    /* call to API */
    prepare_larf_run(larf_variant, datatype, side, m, n, v, incv, tau, c__, ldc, c__out, ldc, work,
                     interfacetype, params);

    /* execution time */
    if(time_min == d_zero)
    {
        time_min = 1e-9;
    }
    /* Performance Computation */
    perf = (double)(2.0 * m * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        perf *= 4.0;
    }

    /* Testing.
      LARF- Leave the vector as it is
      LARF1F- Set the first element of the vector to be 1
      LARF1L- Set the last element of the vector to be 1
    */

    if(larf_variant == FLA_LARF1F)
    {
        assign_value(datatype, get_v_ptr(datatype, v, v_num_elements, 0, incv), 1, 0);
    }
    else if(larf_variant == FLA_LARF1L)
    {
        assign_value(datatype, get_v_ptr(datatype, v, v_num_elements, v_num_elements - 1, incv), 1,
                     0);
    }

    /* Output Validation */
    IF_FLA_BRT_VALIDATION(
        m, n, store_outputs_base(filename, params, 1, 0, datatype, m, n, c__out, ldc),
        validate_larf(tst_api, datatype, side, m, n, v, incv, c__, ldc, c__out, ldc, tau,
                      err_thresh, params),
        check_reproducibility_base(filename, params, 1, 0, datatype, m, n, c__out, ldc))
    else if(!FLA_EXTREME_CASE_TEST)
    {
        validate_larf(tst_api, datatype, side, m, n, v, incv, c__, ldc, c__out, ldc, tau,
                      err_thresh, params);
    }
    else
    {
        if(!check_extreme_value(datatype, m, n, c__out, ldc, params->imatrix_char))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    }

    /* Free up the buffers */
    free_matrix(c__out);
free_buffers:
    FLA_FREE_FILENAME(filename)
    free_matrix(c__);
    free_vector(v);
    free_vector(work);
    free_vector(tau);
}

void prepare_larf_run(fla_larf_variant_t larf_variant, integer datatype, char side, integer m,
                      integer n, void *v, integer incv, void *tau, void *c__, integer ldc__,
                      void *c__out, integer ldc__out, void *work, integer interfacetype,
                      test_params_t *params)
{
    double exe_time;

    FLA_EXEC_LOOP_BEGIN
    {
        copy_matrix(datatype, "full", m, n, c__, ldc__, c__out, ldc__out);

#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call larf CPP API */
            invoke_cpp_larf_common(larf_variant, datatype, &side, &m, &n, v, &incv, tau, c__out,
                                   &ldc__out, work);
            exe_time = fla_test_clock() - exe_time;
        }
        else
#endif
        {
            exe_time = fla_test_clock();
            /* call larf API */
            invoke_larf_common(larf_variant, datatype, &side, &m, &n, v, &incv, tau, c__out,
                               &ldc__out, work);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_NO_INFO
    }
}

/* larf API call interface */

void invoke_larf(integer datatype, char *side, integer *m, integer *n, void *v, integer *incv,
                 void *tau, void *c__, integer *ldc__, void *work)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_slarf(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dlarf(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_clarf(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zlarf(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
    }
}

void invoke_larf1f(integer datatype, char *side, integer *m, integer *n, void *v, integer *incv,
                   void *tau, void *c__, integer *ldc__, void *work)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_slarf1f(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dlarf1f(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_clarf1f(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zlarf1f(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
    }
}

void invoke_larf1l(integer datatype, char *side, integer *m, integer *n, void *v, integer *incv,
                   void *tau, void *c__, integer *ldc__, void *work)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_slarf1l(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dlarf1l(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_clarf1l(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zlarf1l(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
    }
}

void invoke_larf_common(fla_larf_variant_t larf_variant, integer datatype, char *side, integer *m,
                        integer *n, void *v, integer *incv, void *tau, void *c__, integer *ldc__,
                        void *work)
{
    switch(larf_variant)
    {
        case FLA_LARF:
        {
            invoke_larf(datatype, side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case FLA_LARF1F:
        {
            invoke_larf1f(datatype, side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case FLA_LARF1L:
        {
            invoke_larf1l(datatype, side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
    }
}

void invoke_cpp_larf_common(fla_larf_variant_t larf_variant, integer datatype, char *side,
                            integer *m, integer *n, void *v, integer *incv, void *tau, void *c__,
                            integer *ldc__, void *work)
{
    switch(larf_variant)
    {
        case FLA_LARF:
        {
            invoke_cpp_larf(datatype, side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case FLA_LARF1F:
        {
            invoke_cpp_larf1f(datatype, side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case FLA_LARF1L:
        {
            invoke_cpp_larf1l(datatype, side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
    }
}
