/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "test_lapack2flame.h"

int verbose = 0;
double threshold = 0.0, dtime = 0.0, total_dtime = 0.0;

/*! File IO
 */
int read_to_endl(FILE* in, char* str)
{
    char c;
    if (in != NULL)
    {
        while ((c = (char)fgetc(in)) && c != EOF && c != '\n') // && c !=  COMMENT_CHAR)
        {
            *str = c;
            ++str;
        }
        *str = '\0';
        if (c == EOF) return 1;
    }
    return 0;
}
int move_to_endl(FILE* in)
{
    char c;
    if (in != NULL)
    {
        while ((c = (char)fgetc(in)) && c != EOF && c != '\n')
            ;
        if (c == EOF) return 1;
    }
    return 0;
}
int skip_comment(FILE* in)
{
    char c;
    if (in != NULL)
    {
        while((c = (char)fgetc(in)) && c != EOF && c ==  COMMENT_CHAR)
            move_to_endl(in);
        if (c == EOF) return 1;
        fseek(in , -1 , SEEK_CUR);
    }
    return 0;
}

/*! Environment setup
 */
int test_initialize(FLA_Datatype datatype)
{
    FLA_Init();
    switch (datatype)
    {
    case FLA_FLOAT:
    case FLA_COMPLEX:
        threshold = (FLT_EPSILON * EPSILON_RESCALE);
        break;
    case FLA_DOUBLE:
    case FLA_DOUBLE_COMPLEX:
        threshold = (DBL_EPSILON * EPSILON_RESCALE);
        break;
    default:
        CHKERR(-1);
    }

    fprintf(stdout, SINGLE_BAR);
    fprintf(stdout, ">>>> Test lapack2flame initialized\n");
    fprintf(stdout, "     Test threshold = %6.4e\n", threshold);
    fprintf(stdout, SINGLE_BAR);
    return 0;
}
int test_finalize()
{
    fprintf(stdout, SINGLE_BAR);
    fprintf(stdout, ">>>> Test lapack2flame finalized  \n");
    fprintf(stdout, "     Total time : %6f [sec]       \n", total_dtime);
    fprintf(stdout, SINGLE_BAR);
    FLA_Finalize();
    return 0;
}
int test_info( FILE* stream, const char* format, ...)
{
    //fprintf(stream, ">>   %s: ", testname);
    va_list arglist;
    va_start(arglist, format);
    vfprintf(stream, format, arglist);
    va_end(arglist);
    return 0;
}

int flush_test_parameter( param_t* param )
{
    int i;
    param->repeat = 0;
    param->datatype = 0;
    for ( i=0; i<MAX_NUM_PARAMS ; ++i)
    {
        param->trans[i] = 0;
        param->uplos[i] = 0;
        param->dims [i] = 0;
    }
    sprintf( param->testname, " " );
    return 0;
}

int flush_test_result( result_t* result )
{
    flush_test_parameter( &result->param );
    result->performance = 0.0;
    result->residual    = 0.0;
    return 0;
}

int get_io_filenames_from_list( FILE* fp, char* testkind, char* file_in, char* file_out )
{
    char
    str [MAX_STRING_LENGTH],
        kind[MAX_STRING_LENGTH],
        fin [MAX_STRING_LENGTH],
        fout [MAX_STRING_LENGTH];
    int ieof = skip_comment(fp);

    while (ieof == 0)
    {
        ieof = read_to_endl(fp, str) + skip_comment(fp);
        if ( strlen(str) )
        {
            sscanf (str,"%s %s %s", kind, fin, fout);
            if (!strcmp( kind, testkind ))
            {
                strcpy(file_in,  fin);
                strcpy(file_out, fout);
                return 0;
            }
        }
    }
    return 1;
}

int set_test_parameter( char* testname, FLA_Datatype datatype, dim_t repeat,
                        dim_t m, dim_t n, dim_t k,
                        FLA_Trans transa, FLA_Trans transb, FLA_Trans transc, FLA_Trans transd,
                        FLA_Uplo uploa, FLA_Uplo uplob, FLA_Uplo uploc, FLA_Uplo uplod,
                        param_t *param )
{
    param->datatype = datatype;
    param->repeat   = repeat;
    param->trans[0] = transa;
    param->trans[1] = transb;
    param->trans[2] = transc;
    param->trans[3] = transd;
    param->uplos[0] = uploa;
    param->uplos[1] = uplob;
    param->uplos[2] = uploc;
    param->uplos[3] = uplod;
    param->dims [0] = m;
    param->dims [1] = n;
    param->dims [2] = k;
    strcpy(param->testname, testname);

    return 0;
}


int run( int (*func)(FILE*,param_t,result_t*),
         FILE* stream, param_t param, result_t* result )
{
    flush_test_result(result);
    result->param = param;

    if (verbose)
    {
        fprintf( stream, SINGLE_BAR );
        fprintf( stream, ">>>> Testing ... %s, %s \n", datatype2str(param.datatype), param.testname );
        show_test_param( stream, param );
    }
    dtime = FLA_Clock();
    func( stream, param, result );
    dtime = FLA_Clock() - dtime;
    total_dtime += dtime;
    if (verbose)
    {
        fprintf( stream, "     Time : %6.4f [sec] \n", dtime);
        fprintf( stream, SINGLE_BAR );
    }
    if (result->residual > threshold)
    {
        show_test_param ( stream, result->param );
        show_test_result( stream, *result );
    }
    return 0;
}

int show_test_param( FILE* stream, param_t param )
{
    int i;
    fprintf( stream, SINGLE_BAR );
    fprintf( stream, ">>>> Parameters: %s, %s\n", datatype2str(param.datatype), param.testname );

    fprintf( stream, "TRANS    = ");
    for (i=0; i<MAX_NUM_PARAMS; ++i)
        fprintf( stream,"%5d", (int)param.trans[i]);
    fprintf( stream, "\n");

    fprintf( stream, "UPLOS    = ");
    for (i=0; i<MAX_NUM_PARAMS; ++i)
        fprintf( stream,"%5d", (int)param.uplos[i]);
    fprintf( stream, "\n");

    fprintf( stream, "DIMS     = ");
    for (i=0; i<MAX_NUM_PARAMS; ++i)
        fprintf( stream,"%5d", (int)param.dims[i]);
    fprintf( stream, "\n");

    fprintf( stream, SINGLE_BAR );
    return 0;
}

int show_test_result( FILE* stream, result_t result )
{
    fprintf( stream, SINGLE_BAR );
    fprintf( stream, ">>>> Test results: %s, %s\n", datatype2str(result.param.datatype), result.param.testname );
    fprintf( stream, "FLOPs    = %6.4e [GFLOPs]\n", result.performance );
    fprintf( stream, "Residual = %6.4e \n", result.residual );
    fprintf( stream, SINGLE_BAR );
    return 0;
}

int create_testitems( const char* name, testitem_t* item)
{
    FILE* fp;
    param_t* param;
    char str[MAX_STRING_LENGTH];
    int i, ieof = 0, k1, k2, k3, mm, nn, kk,
           m[MAX_NUM_TESTITEMS], n[MAX_NUM_TESTITEMS], k[MAX_NUM_TESTITEMS];
    fp = fopen(name, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error opening the file: %s\n", name);
        return -1;
    }
    else
    {
        ieof = skip_comment(fp);

        /* M */
        ieof = read_to_endl(fp, str) + skip_comment(fp);
        sscanf (str,"%d", &mm);
        for (i=0; i<mm; ++i)
            fscanf(fp,"%d", &m[i]);
        ieof = move_to_endl(fp) + skip_comment(fp);

        /* N */
        ieof = read_to_endl(fp, str) + skip_comment(fp);
        sscanf (str,"%d", &nn);
        for (i=0; i<nn; ++i)
            fscanf(fp,"%d", &n[i]);
        ieof = move_to_endl(fp) + skip_comment(fp);

        /* K */
        ieof = read_to_endl(fp, str) + skip_comment(fp);
        sscanf (str,"%d", &kk);
        for (i=0; i<kk; ++i)
            fscanf(fp,"%d", &k[i]);
        ieof = move_to_endl(fp) + skip_comment(fp);
    }
    fclose(fp);

    /* Generate parameters */
    item->size = mm*nn*kk;
    item->results = (result_t*) malloc(item->size*sizeof(result_t));

    i=0;
    for (k1=0; k1<kk; ++k1)
        for (k2=0; k2<nn; ++k2)
            for (k3=0; k3<mm; ++k3)
            {
                param = &item->results[i].param;
                flush_test_parameter( param );
                param->dims[0] = m[k3];
                param->dims[1] = n[k2];
                param->dims[2] = k[k1];
                ++i;
            }
    return 0;
}

int free_testitems( testitem_t* item)
{
    if (item->size)
        free(item->results);
    item->size = 0;
    return 0;
}

int write_test_result_to_file( testitem_t items, char* name )
{
    int i, size = items.size;
    result_t result;
    FILE* fp;

    fp = fopen(name, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Error opening the file: %s\n", name);
        return -1;
    }
    else
    {
        for (i=0; i<size; ++i)
        {
            result = items.results[i];
            fprintf(fp, "  %5d, %5d, %5d,    %10.6e\n",
                    (int)result.param.dims[0], (int)result.param.dims[1], (int)result.param.dims[2],
                    result.performance);
        }
    }
    fclose(fp);

    return 0;
}

