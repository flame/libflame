#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define sytrf_free() \
       free (a   ); \
       free (aref); \
       free (ipiv); \
       free (ipivref)

/* Begin sytrf_double_parameters  class definition */
class sytrf_double_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      double *a,*aref; //The array ab contains the matrix A

      /* Output parameters */
      lapack_int *ipiv,*ipivref; // The ipivot indices
      /* Return Values */
      lapack_int info,inforef;

   public: 
      sytrf_double_parameters ( int matrix_layout_i,char uplo_i,
                    lapack_int n_i,lapack_int lda_i);
              
      ~sytrf_double_parameters (); 
};  /* end of sytrf_double_parameters  class definition */


/* Constructor sytrf_double_parameters definition */
sytrf_double_parameters:: sytrf_double_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (lda*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       sytrf_free();
       printf(" sytrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_double_parameters:: ~sytrf_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_double_parameters object: destructor invoked. \n");
#endif
   sytrf_free();
}

TEST(sytrf,dsytrf1) {

    /* LAPACKE DSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf) ( int matrix_layout,char uplo,
                                 lapack_int n,double *a,lapack_int lda,
                                                       lapack_int *ipiv);

    Fptr_NL_LAPACKE_dsytrf DSYTRF;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    sytrf_double_parameters   dsytrf_obj(LAPACK_ROW_MAJOR,'U',451,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DSYTRF = (Fptr_NL_LAPACKE_dsytrf)dlsym(hModule,"LAPACKE_dsytrf");
    if (NULL == DSYTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dsytrf_obj.inforef = DSYTRF( dsytrf_obj.matrix_layout,dsytrf_obj.uplo,
                             dsytrf_obj.n,dsytrf_obj.aref,dsytrf_obj.lda,
                                                       dsytrf_obj.ipivref);

    /* Compute libflame's Lapacke o/p  */
    dsytrf_obj.info = LAPACKE_dsytrf( dsytrf_obj.matrix_layout,dsytrf_obj.uplo,
                                     dsytrf_obj.n,dsytrf_obj.a,dsytrf_obj.lda,
                                                               dsytrf_obj.ipiv);

    if( dsytrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrf is wrong\n",
                    dsytrf_obj.info );
    }
    if( dsytrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrf is wrong\n",
        dsytrf_obj.inforef );
    }
    ipiv_diff = computeDiff_i( dsytrf_obj.n,dsytrf_obj.ipiv,dsytrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dsytrf_obj.n*dsytrf_obj.lda,dsytrf_obj.a,dsytrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_float_parameters  class definition */
class sytrf_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      float *a,*aref; //The array ab contains the matrix A

      /* Output parameters */
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      sytrf_float_parameters ( int matrix_layout_i,char uplo_i,
                  lapack_int n_i,lapack_int lda_i);
      ~sytrf_float_parameters (); 
};  /* end of sytrf_float_parameters  class definition */


/* Constructor sytrf_float_parameters definition */
sytrf_float_parameters:: sytrf_float_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       sytrf_free();
       printf(" sytrf_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_float_parameters:: ~sytrf_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_float_parameters object: destructor invoked. \n");
#endif
   sytrf_free();
}

TEST(sytrf,ssytrf1) {

    /* LAPACKE SSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf) ( int matrix_layout ,char uplo,
                                 lapack_int n ,float *a ,lapack_int lda,
                                                        lapack_int *ipiv);

    Fptr_NL_LAPACKE_ssytrf SSYTRF;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;

    sytrf_float_parameters   ssytrf_obj(LAPACK_COL_MAJOR,'U',1020,1041);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SSYTRF = (Fptr_NL_LAPACKE_ssytrf)dlsym(hModule,"LAPACKE_ssytrf");
    if (NULL == SSYTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    ssytrf_obj.inforef = SSYTRF( ssytrf_obj.matrix_layout,ssytrf_obj.uplo,
                            ssytrf_obj.n,ssytrf_obj.aref,ssytrf_obj.lda,
                                                     ssytrf_obj.ipivref );

        /* Compute libflame's Lapacke o/p  */
    ssytrf_obj.info = LAPACKE_ssytrf(ssytrf_obj.matrix_layout,ssytrf_obj.uplo,
                                    ssytrf_obj.n,ssytrf_obj.a,ssytrf_obj.lda,
                                                              ssytrf_obj.ipiv);

    if( ssytrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrf is wrong\n",
                    ssytrf_obj.info );
    }
    if( ssytrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrf is wrong\n",
        ssytrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( ssytrf_obj.n,ssytrf_obj.ipiv,ssytrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( ssytrf_obj.n*ssytrf_obj.lda,ssytrf_obj.a,
                                                   ssytrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_scomplex_parameters  class definition */
class sytrf_scomplex_parameters{

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
      sytrf_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                 lapack_int n_i,lapack_int lda_i);
      ~sytrf_scomplex_parameters (); 
};  /* end of sytrf_scomplex_parameters  class definition */


/* Constructor sytrf_scomplex_parameters definition */
sytrf_scomplex_parameters:: sytrf_scomplex_parameters ( int matrix_layout_i,
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
       sytrf_free();
       printf(" sytrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_scomplex_parameters:: ~sytrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrf_free();
}

TEST(sytrf,csytrf1) {

    /* LAPACKE CSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf) ( int matrix_layout ,char uplo ,
                   lapack_int n ,lapack_complex_float *a ,lapack_int lda,
                                                      lapack_int *ipiv);

    Fptr_NL_LAPACKE_csytrf CSYTRF;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    sytrf_scomplex_parameters   csytrf_obj(LAPACK_COL_MAJOR,'U',510,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CSYTRF = (Fptr_NL_LAPACKE_csytrf)dlsym(hModule,"LAPACKE_csytrf");
    if (NULL == CSYTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    csytrf_obj.inforef = CSYTRF( csytrf_obj.matrix_layout,csytrf_obj.uplo,
                            csytrf_obj.n,csytrf_obj.aref,csytrf_obj.lda,
                                                      csytrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    csytrf_obj.info = LAPACKE_csytrf(csytrf_obj.matrix_layout,csytrf_obj.uplo,
                                    csytrf_obj.n,csytrf_obj.a,csytrf_obj.lda,
                                                             csytrf_obj.ipiv );

    if( csytrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrf is wrong\n",
                    csytrf_obj.info );
    }
    if( csytrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrf is wrong\n",
        csytrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( csytrf_obj.n,csytrf_obj.ipiv,csytrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( csytrf_obj.n*csytrf_obj.lda,csytrf_obj.a,csytrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_dcomplex_parameters  class definition */
class sytrf_dcomplex_parameters{

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
      sytrf_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~sytrf_dcomplex_parameters (); 
};  /* end of sytrf_dcomplex_parameters  class definition */


/* Constructor sytrf_dcomplex_parameters definition */
sytrf_dcomplex_parameters:: sytrf_dcomplex_parameters ( int matrix_layout_i,
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
       sytrf_free();
       printf(" sytrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_dcomplex_parameters:: ~sytrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrf_free();
}

TEST(sytrf,zsytrf1) {

    /* LAPACKE ZSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf) ( int matrix_layout ,char uplo ,
                  lapack_int n ,lapack_complex_double *a ,lapack_int lda,
                        lapack_int *ipiv);

    Fptr_NL_LAPACKE_zsytrf ZSYTRF;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    sytrf_dcomplex_parameters   zsytrf_obj(LAPACK_ROW_MAJOR,'U',100,121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZSYTRF = (Fptr_NL_LAPACKE_zsytrf)dlsym(hModule,"LAPACKE_zsytrf");
    if (NULL == ZSYTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zsytrf_obj.inforef = ZSYTRF( zsytrf_obj.matrix_layout,zsytrf_obj.uplo,
                            zsytrf_obj.n,zsytrf_obj.aref,zsytrf_obj.lda,
                                                     zsytrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zsytrf_obj.info = LAPACKE_zsytrf(zsytrf_obj.matrix_layout,zsytrf_obj.uplo,
                                   zsytrf_obj.n,zsytrf_obj.a,zsytrf_obj.lda,
                                                             zsytrf_obj.ipiv);

    if( zsytrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrf is wrong\n",
                    zsytrf_obj.info );
    }
    if( zsytrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrf is wrong\n",
        zsytrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zsytrf_obj.n,zsytrf_obj.ipiv,zsytrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zsytrf_obj.n*zsytrf_obj.lda,zsytrf_obj.a,zsytrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
