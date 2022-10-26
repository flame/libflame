#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define hetrf_free() \
       free (a   ); \
       free (aref); \
       free (ipiv); \
       free (ipivref)


/* Begin hetrf_scomplex_parameters  class definition */
class hetrf_scomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_float *a,*aref; //The array ab contains the matrix A
      /* Output parameters */
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      hetrf_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                 lapack_int n_i,lapack_int lda_i);
      ~hetrf_scomplex_parameters (); 
};  /* end of hetrf_scomplex_parameters  class definition */


/* Constructor hetrf_scomplex_parameters definition */
hetrf_scomplex_parameters:: hetrf_scomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                            lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       hetrf_free();
       printf(" hetrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

hetrf_scomplex_parameters:: ~hetrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrf_scomplex_parameters object: destructor invoked. \n");
#endif
   hetrf_free();
}

TEST(hetrf,chetrf1) {

    /* LAPACKE CHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf) ( int matrix_layout ,char uplo ,
                   lapack_int n ,lapack_complex_float *a ,lapack_int lda,
                                                      lapack_int *ipiv);

    Fptr_NL_LAPACKE_chetrf CHETRF;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    hetrf_scomplex_parameters   chetrf_obj(LAPACK_COL_MAJOR,'U',510,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CHETRF = (Fptr_NL_LAPACKE_chetrf)dlsym(hModule,"LAPACKE_chetrf");
    if (NULL == CHETRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chetrf_obj.inforef = CHETRF( chetrf_obj.matrix_layout,chetrf_obj.uplo,
                            chetrf_obj.n,chetrf_obj.aref,chetrf_obj.lda,
                                                      chetrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    chetrf_obj.info = LAPACKE_chetrf(chetrf_obj.matrix_layout,chetrf_obj.uplo,
                                    chetrf_obj.n,chetrf_obj.a,chetrf_obj.lda,
                                                             chetrf_obj.ipiv );

    if( chetrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrf is wrong\n",
                    chetrf_obj.info );
    }
    if( chetrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrf is wrong\n",
        chetrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( chetrf_obj.n,chetrf_obj.ipiv,chetrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( chetrf_obj.n*chetrf_obj.lda,chetrf_obj.a,chetrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST(hetrf,chetrf2) {

    /* LAPACKE CHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf) ( int matrix_layout ,char uplo ,
                   lapack_int n ,lapack_complex_float *a ,lapack_int lda,
                                                      lapack_int *ipiv);

    Fptr_NL_LAPACKE_chetrf CHETRF;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    hetrf_scomplex_parameters   chetrf_obj(LAPACK_COL_MAJOR,'L',421,478);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CHETRF = (Fptr_NL_LAPACKE_chetrf)dlsym(hModule,"LAPACKE_chetrf");
    if (NULL == CHETRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chetrf_obj.inforef = CHETRF( chetrf_obj.matrix_layout,chetrf_obj.uplo,
                            chetrf_obj.n,chetrf_obj.aref,chetrf_obj.lda,
                                                      chetrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    chetrf_obj.info = LAPACKE_chetrf(chetrf_obj.matrix_layout,chetrf_obj.uplo,
                                    chetrf_obj.n,chetrf_obj.a,chetrf_obj.lda,
                                                             chetrf_obj.ipiv );

    if( chetrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrf is wrong\n",
                    chetrf_obj.info );
    }
    if( chetrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrf is wrong\n",
        chetrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( chetrf_obj.n,chetrf_obj.ipiv,chetrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( chetrf_obj.n*chetrf_obj.lda,chetrf_obj.a,chetrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin hetrf_dcomplex_parameters  class definition */
class hetrf_dcomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_double *a,*aref; //The array ab contains the matrix A
      /* Output parameters */
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      hetrf_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~hetrf_dcomplex_parameters (); 
};  /* end of hetrf_dcomplex_parameters  class definition */


/* Constructor hetrf_dcomplex_parameters definition */
hetrf_dcomplex_parameters:: hetrf_dcomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       hetrf_free();
       printf(" hetrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

hetrf_dcomplex_parameters:: ~hetrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetrf_free();
}

TEST(hetrf,zhetrf1) {

    /* LAPACKE ZHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf) ( int matrix_layout ,char uplo ,
                  lapack_int n ,lapack_complex_double *a ,lapack_int lda,
                        lapack_int *ipiv);

    Fptr_NL_LAPACKE_zhetrf ZHETRF;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    hetrf_dcomplex_parameters   zhetrf_obj(LAPACK_ROW_MAJOR,'U',100,121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZHETRF = (Fptr_NL_LAPACKE_zhetrf)dlsym(hModule,"LAPACKE_zhetrf");
    if (NULL == ZHETRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhetrf_obj.inforef = ZHETRF( zhetrf_obj.matrix_layout,zhetrf_obj.uplo,
                            zhetrf_obj.n,zhetrf_obj.aref,zhetrf_obj.lda,
                                                     zhetrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zhetrf_obj.info = LAPACKE_zhetrf(zhetrf_obj.matrix_layout,zhetrf_obj.uplo,
                                   zhetrf_obj.n,zhetrf_obj.a,zhetrf_obj.lda,
                                                             zhetrf_obj.ipiv);

    if( zhetrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrf is wrong\n",
                    zhetrf_obj.info );
    }
    if( zhetrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrf is wrong\n",
        zhetrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zhetrf_obj.n,zhetrf_obj.ipiv,zhetrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zhetrf_obj.n*zhetrf_obj.lda,zhetrf_obj.a,zhetrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}


TEST(hetrf,zhetrf2) {

    /* LAPACKE ZHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf) ( int matrix_layout ,char uplo ,
                  lapack_int n ,lapack_complex_double *a ,lapack_int lda,
                        lapack_int *ipiv);

    Fptr_NL_LAPACKE_zhetrf ZHETRF;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    hetrf_dcomplex_parameters   zhetrf_obj(LAPACK_ROW_MAJOR,'L',1000,1210);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZHETRF = (Fptr_NL_LAPACKE_zhetrf)dlsym(hModule,"LAPACKE_zhetrf");
    if (NULL == ZHETRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhetrf_obj.inforef = ZHETRF( zhetrf_obj.matrix_layout,zhetrf_obj.uplo,
                            zhetrf_obj.n,zhetrf_obj.aref,zhetrf_obj.lda,
                                                     zhetrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zhetrf_obj.info = LAPACKE_zhetrf(zhetrf_obj.matrix_layout,zhetrf_obj.uplo,
                                   zhetrf_obj.n,zhetrf_obj.a,zhetrf_obj.lda,
                                                             zhetrf_obj.ipiv);

    if( zhetrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrf is wrong\n",
                    zhetrf_obj.info );
    }
    if( zhetrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrf is wrong\n",
        zhetrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zhetrf_obj.n,zhetrf_obj.ipiv,zhetrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zhetrf_obj.n*zhetrf_obj.lda,zhetrf_obj.a,zhetrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
