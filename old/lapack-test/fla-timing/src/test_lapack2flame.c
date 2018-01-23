/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "test_lapack2flame.h"
#define TESTLIST "test.list"
#define TESTUPLO FLA_LOWER_TRIANGULAR



#define RUNTEST(test_name, test_func,                                   \
                transa, transb, transc, transd,                         \
                uploa,  uplob,  uploc,  uplod)                          \
  rewind(fp);                                                           \
  while ( ! get_io_filenames_from_list(fp, #test_name, file_in, file_out) ) { \
    if (! create_testitems( file_in, &items )) {                        \
      for (i=0;i<items.size;++i) {                                      \
        param = items.results[i].param;                                 \
        ierr = set_test_parameter( #test_name, str2datatype(datatype), repeat, \
                                   param.dims[0], param.dims[1], param.dims[2], \
                                   transa, transb, transc, transd,      \
                                   uploa,  uplob,  uploc,  uplod,       \
                                   &param );CHKERR(ierr);               \
        ierr = run( test_func, stdout, param, &items.results[i] );CHKERR(ierr); \
        PROGRESS(param.testname, i, items.size);                        \
      }                                                                 \
      sprintf(file_name, "%s_%s", datatype, file_out);                  \
      ierr = write_test_result_to_file( items, file_name );CHKERR(ierr); \
      ierr = free_testitems( &items );CHKERR(ierr);                     \
    }                                                                   \
  }


int main( int argc, char** argv )
{
    testitem_t  items;
    param_t     param;
    char
    datatype [ MAX_NAME_LENGTH ],
             file_in  [ MAX_NAME_LENGTH ],
             file_out [ MAX_NAME_LENGTH ],
             file_name[ MAX_NAME_LENGTH ];
    int         i, ierr, repeat;
    FILE*       fp;

    // Display direction
    if (argc == 3)
    {
        strcpy(datatype, argv[1]);
        repeat = atoi(argv[2]);
    }
    else
    {
        fprintf(stderr, "       \n");
        fprintf(stderr, "Usage: %s datatype repeat\n", argv[0]);
        fprintf(stderr, "       datatype : single, double, scomplex, dcomplex\n");
        fprintf(stderr, "       repeat   : # of test repeat\n");
        fprintf(stderr, "       \n");
        return -1;
    }

    // Initialize libflame
    ierr = test_initialize(str2datatype(datatype));
    CHKERR(ierr);

    // Open test list
    fp = fopen(TESTLIST,"r");
    CHKERR(fp == NULL);

    // ** GEMM to get reference performance
    RUNTEST( GEMM, test_gemm,
             FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, NOT_DEFINED, NOT_DEFINED,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED );

    // ** CHOL to get reference performance
    RUNTEST( CHOL, test_chol,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED,
             TESTUPLO, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED );

    // ** TRINV to get reference performance
    RUNTEST( TRINV, test_trinv,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED,
             TESTUPLO, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED );

    // ** LUPIV to get reference performance
    RUNTEST( LUPIV, test_lu_piv,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED );

    // ** QR to get reference performance
    RUNTEST( QR, test_qr,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED );

    // ** LQ to get reference performance
    RUNTEST( LQ, test_lq,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED );

    // ** APPQ to get reference performance
    RUNTEST( APPQ, test_appq,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED );

    // ** SVD to get reference performance
    RUNTEST( SVD, test_svd,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED,
             NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED );

    /* // ** SYLV to get reference performance */
    /* RUNTEST( SYLV, test_sylv, */
    /*          FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, NOT_DEFINED, NOT_DEFINED, */
    /*          NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED ); */


    /* // ** HESS to get reference performance */
    /* RUNTEST( HESS, test_hess, */
    /*          NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, */
    /*          NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED ); */

    /* // ** TRIDIAG to get reference performance */
    /* RUNTEST( TRIDIAG, test_tridiag, */
    /*          NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, */
    /*          NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED ); */


    /* // ** HEVD to get reference performance */
    /* RUNTEST( HEVD, test_hevd, */
    /*          NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, */
    /*          NOT_DEFINED, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED ); */


    // Close the test list
    ierr = fclose(fp);
    CHKERR(ierr);

    // Finalize libflame
    ierr = test_finalize();
    CHKERR(ierr);

    return 0;
}



