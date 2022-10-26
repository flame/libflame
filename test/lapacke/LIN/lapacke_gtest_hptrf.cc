#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define hptrf_free() \
       free (a   ); \
       free (aref); \
       free (ipiv); \
       free (ipivref)


/* Begin hptrf_scomplex_parameters  class definition */
class hptrf_scomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      lapack_complex_float *a,*aref; //The array ab contains the matrix A
      /* Output parameters */
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      hptrf_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                 lapack_int n_i);
      ~hptrf_scomplex_parameters (); 
};  /* end of hptrf_scomplex_parameters  class definition */


/* Constructor hptrf_scomplex_parameters definition */
hptrf_scomplex_parameters:: hptrf_scomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i ){ 
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*(n+1)/2)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       hptrf_free();
       printf(" hptrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,(n*(n+1)/2));

   } /* end of Constructor  */

hptrf_scomplex_parameters:: ~hptrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hptrf_scomplex_parameters object: destructor invoked. \n");
#endif
   hptrf_free();
}

TEST(hptrf,chptrf1) {

    /* LAPACKE CHPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chptrf) ( int matrix_layout,char uplo,
                   lapack_int n,lapack_complex_float *a,lapack_int *ipiv);

    Fptr_NL_LAPACKE_chptrf CHPTRF;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    hptrf_scomplex_parameters   chptrf_obj(LAPACK_COL_MAJOR,'U',510);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CHPTRF = (Fptr_NL_LAPACKE_chptrf)dlsym(hModule,"LAPACKE_chptrf");
    if (NULL == CHPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chptrf_obj.inforef = CHPTRF( chptrf_obj.matrix_layout,chptrf_obj.uplo,
                                             chptrf_obj.n,chptrf_obj.aref,
                                                      chptrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    chptrf_obj.info = LAPACKE_chptrf(chptrf_obj.matrix_layout,chptrf_obj.uplo,
                                                    chptrf_obj.n,chptrf_obj.a,
                                                             chptrf_obj.ipiv );

    if( chptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chptrf is wrong\n",
                    chptrf_obj.info );
    }
    if( chptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chptrf is wrong\n",
        chptrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( chptrf_obj.n,chptrf_obj.ipiv,chptrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhptrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( chptrf_obj.n*(chptrf_obj.n + 1)/2,chptrf_obj.a,chptrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST(hptrf,chptrf2) {

    /* LAPACKE CHPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chptrf) ( int matrix_layout,char uplo,
                   lapack_int n,lapack_complex_float *a,lapack_int *ipiv);

    Fptr_NL_LAPACKE_chptrf CHPTRF;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    hptrf_scomplex_parameters   chptrf_obj(LAPACK_ROW_MAJOR,'L',210);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CHPTRF = (Fptr_NL_LAPACKE_chptrf)dlsym(hModule,"LAPACKE_chptrf");
    if (NULL == CHPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chptrf_obj.inforef = CHPTRF( chptrf_obj.matrix_layout,chptrf_obj.uplo,
                                             chptrf_obj.n,chptrf_obj.aref,
                                                      chptrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    chptrf_obj.info = LAPACKE_chptrf(chptrf_obj.matrix_layout,chptrf_obj.uplo,
                                                    chptrf_obj.n,chptrf_obj.a,
                                                             chptrf_obj.ipiv );

    if( chptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chptrf is wrong\n",
                    chptrf_obj.info );
    }
    if( chptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chptrf is wrong\n",
        chptrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( chptrf_obj.n,chptrf_obj.ipiv,chptrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhptrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( chptrf_obj.n*(chptrf_obj.n + 1)/2,chptrf_obj.a,chptrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}


/* Begin hptrf_dcomplex_parameters  class definition */
class hptrf_dcomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      lapack_complex_double *a,*aref; //The array ab contains the matrix A
      /* Output parameters */
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      hptrf_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i);
      ~hptrf_dcomplex_parameters (); 
};  /* end of hptrf_dcomplex_parameters  class definition */


/* Constructor hptrf_dcomplex_parameters definition */
hptrf_dcomplex_parameters:: hptrf_dcomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i ){
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*(n+1)/2)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       hptrf_free();
       printf(" hptrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,(n*(n+1)/2));

   } /* end of Constructor  */

hptrf_dcomplex_parameters:: ~hptrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hptrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   hptrf_free();
}

TEST(hptrf,zhptrf1) {

    /* LAPACKE ZHPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptrf) ( int matrix_layout,char uplo,
                                  lapack_int n,lapack_complex_double *a,
                                                       lapack_int *ipiv);

    Fptr_NL_LAPACKE_zhptrf ZHPTRF;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    hptrf_dcomplex_parameters   zhptrf_obj(LAPACK_ROW_MAJOR,'U',121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZHPTRF = (Fptr_NL_LAPACKE_zhptrf)dlsym(hModule,"LAPACKE_zhptrf");
    if (NULL == ZHPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhptrf_obj.inforef = ZHPTRF( zhptrf_obj.matrix_layout,zhptrf_obj.uplo,
                                             zhptrf_obj.n,zhptrf_obj.aref,
                                                     zhptrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zhptrf_obj.info = LAPACKE_zhptrf(zhptrf_obj.matrix_layout,zhptrf_obj.uplo,
                                                    zhptrf_obj.n,zhptrf_obj.a,
                                                             zhptrf_obj.ipiv);

    if( zhptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhptrf is wrong\n",
                    zhptrf_obj.info );
    }
    if( zhptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhptrf is wrong\n",
        zhptrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zhptrf_obj.n,zhptrf_obj.ipiv,zhptrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhptrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zhptrf_obj.n*(zhptrf_obj.n + 1)/2,zhptrf_obj.a,zhptrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}


TEST(hptrf,zhptrf2) {

    /* LAPACKE ZHPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptrf) ( int matrix_layout,char uplo,
                                  lapack_int n,lapack_complex_double *a,
                                                       lapack_int *ipiv);

    Fptr_NL_LAPACKE_zhptrf ZHPTRF;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    hptrf_dcomplex_parameters   zhptrf_obj(LAPACK_COL_MAJOR,'L',121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZHPTRF = (Fptr_NL_LAPACKE_zhptrf)dlsym(hModule,"LAPACKE_zhptrf");
    if (NULL == ZHPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhptrf_obj.inforef = ZHPTRF( zhptrf_obj.matrix_layout,zhptrf_obj.uplo,
                                             zhptrf_obj.n,zhptrf_obj.aref,
                                                     zhptrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zhptrf_obj.info = LAPACKE_zhptrf(zhptrf_obj.matrix_layout,zhptrf_obj.uplo,
                                                    zhptrf_obj.n,zhptrf_obj.a,
                                                             zhptrf_obj.ipiv);

    if( zhptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhptrf is wrong\n",
                    zhptrf_obj.info );
    }
    if( zhptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhptrf is wrong\n",
        zhptrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zhptrf_obj.n,zhptrf_obj.ipiv,zhptrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhptrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zhptrf_obj.n*(zhptrf_obj.n + 1)/2,zhptrf_obj.a,zhptrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
