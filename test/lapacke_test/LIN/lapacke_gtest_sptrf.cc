#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define sptrf_free() \
       free (a   ); \
       free (aref); \
       free (ipiv); \
       free (ipivref)

/* Begin sptrf_double_parameters  class definition */
class sptrf_double_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      double *a,*aref; //The array ab contains the matrix A

      /* Output parameters */
      lapack_int *ipiv,*ipivref; // The ipivot indices
      /* Return Values */
      lapack_int info,inforef;

   public: 
      sptrf_double_parameters ( int matrix_layout_i,char uplo_i,
                    lapack_int n_i);
              
      ~sptrf_double_parameters (); 
};  /* end of sptrf_double_parameters  class definition */


/* Constructor sptrf_double_parameters definition */
sptrf_double_parameters:: sptrf_double_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*(n+1)/2));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       sptrf_free();
       printf(" sptrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,(n*(n+1)/2));

   } /* end of Constructor  */

sptrf_double_parameters:: ~sptrf_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptrf_double_parameters object: destructor invoked. \n");
#endif
   sptrf_free();
}

TEST(sptrf,dsptrf1) {

    /* LAPACKE DSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dsptrf) ( int matrix_layout,char uplo,
                                 lapack_int n,double *ap,lapack_int *ipiv);

    Fptr_NL_LAPACKE_dsptrf DSPTRF;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    sptrf_double_parameters   dsptrf_obj(LAPACK_ROW_MAJOR,'U',521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DSPTRF = (Fptr_NL_LAPACKE_dsptrf)dlsym(hModule,"LAPACKE_dsptrf");
    if (NULL == DSPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dsptrf_obj.inforef = DSPTRF( dsptrf_obj.matrix_layout,dsptrf_obj.uplo,
                                             dsptrf_obj.n,dsptrf_obj.aref,
                                                       dsptrf_obj.ipivref);

    /* Compute libflame's Lapacke o/p  */
    dsptrf_obj.info = LAPACKE_dsptrf( dsptrf_obj.matrix_layout,dsptrf_obj.uplo,
                                                     dsptrf_obj.n,dsptrf_obj.a,
                                                               dsptrf_obj.ipiv);

    if( dsptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsptrf is wrong\n",
                    dsptrf_obj.info );
    }
    if( dsptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsptrf is wrong\n",
        dsptrf_obj.inforef );
    }
    ipiv_diff = computeDiff_i( dsptrf_obj.n,dsptrf_obj.ipiv,dsptrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsptrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dsptrf_obj.n*(dsptrf_obj.n + 1)/2,dsptrf_obj.a,dsptrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sptrf_float_parameters  class definition */
class sptrf_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      float *a,*aref; //The array ab contains the matrix A

      /* Output parameters */
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      sptrf_float_parameters ( int matrix_layout_i,char uplo_i,
                  lapack_int n_i);
      ~sptrf_float_parameters (); 
};  /* end of sptrf_float_parameters  class definition */


/* Constructor sptrf_float_parameters definition */
sptrf_float_parameters:: sptrf_float_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i ){ 
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*(n+1)/2)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       sptrf_free();
       printf(" sptrf_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,(n*(n+1)/2));

   } /* end of Constructor  */

sptrf_float_parameters:: ~sptrf_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptrf_float_parameters object: destructor invoked. \n");
#endif
   sptrf_free();
}

TEST(sptrf,ssptrf1) {

    /* LAPACKE SSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_ssptrf) ( int matrix_layout,char uplo,
                                 lapack_int n,float *a,lapack_int *ipiv);

    Fptr_NL_LAPACKE_ssptrf SSPTRF;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;

    sptrf_float_parameters   ssptrf_obj(LAPACK_COL_MAJOR,'U',1020);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SSPTRF = (Fptr_NL_LAPACKE_ssptrf)dlsym(hModule,"LAPACKE_ssptrf");
    if (NULL == SSPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    ssptrf_obj.inforef = SSPTRF( ssptrf_obj.matrix_layout,ssptrf_obj.uplo,
                            ssptrf_obj.n,ssptrf_obj.aref,ssptrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    ssptrf_obj.info = LAPACKE_ssptrf(ssptrf_obj.matrix_layout,ssptrf_obj.uplo,
                                  ssptrf_obj.n,ssptrf_obj.a,ssptrf_obj.ipiv);

    if( ssptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssptrf is wrong\n",
                    ssptrf_obj.info );
    }
    if( ssptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssptrf is wrong\n",
        ssptrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( ssptrf_obj.n,ssptrf_obj.ipiv,ssptrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsptrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( ssptrf_obj.n*(ssptrf_obj.n + 1)/2,ssptrf_obj.a,
                                                   ssptrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sptrf_scomplex_parameters  class definition */
class sptrf_scomplex_parameters{

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
      sptrf_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                 lapack_int n_i);
      ~sptrf_scomplex_parameters (); 
};  /* end of sptrf_scomplex_parameters  class definition */


/* Constructor sptrf_scomplex_parameters definition */
sptrf_scomplex_parameters:: sptrf_scomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i ){ 
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*(n+1)/2)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       sptrf_free();
       printf(" sptrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,(n*(n+1)/2));

   } /* end of Constructor  */

sptrf_scomplex_parameters:: ~sptrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptrf_scomplex_parameters object: destructor invoked. \n");
#endif
   sptrf_free();
}

TEST(sptrf,csptrf1) {

    /* LAPACKE CSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_csptrf) ( int matrix_layout,char uplo,
                   lapack_int n,lapack_complex_float *a,lapack_int *ipiv);

    Fptr_NL_LAPACKE_csptrf CSPTRF;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    sptrf_scomplex_parameters   csptrf_obj(LAPACK_COL_MAJOR,'U',510);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CSPTRF = (Fptr_NL_LAPACKE_csptrf)dlsym(hModule,"LAPACKE_csptrf");
    if (NULL == CSPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    csptrf_obj.inforef = CSPTRF( csptrf_obj.matrix_layout,csptrf_obj.uplo,
                                             csptrf_obj.n,csptrf_obj.aref,
                                                      csptrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    csptrf_obj.info = LAPACKE_csptrf(csptrf_obj.matrix_layout,csptrf_obj.uplo,
                                                    csptrf_obj.n,csptrf_obj.a,
                                                             csptrf_obj.ipiv );

    if( csptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csptrf is wrong\n",
                    csptrf_obj.info );
    }
    if( csptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csptrf is wrong\n",
        csptrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( csptrf_obj.n,csptrf_obj.ipiv,csptrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsptrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( csptrf_obj.n*(csptrf_obj.n + 1)/2,csptrf_obj.a,csptrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sptrf_dcomplex_parameters  class definition */
class sptrf_dcomplex_parameters{

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
      sptrf_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i);
      ~sptrf_dcomplex_parameters (); 
};  /* end of sptrf_dcomplex_parameters  class definition */


/* Constructor sptrf_dcomplex_parameters definition */
sptrf_dcomplex_parameters:: sptrf_dcomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i ){
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*(n+1)/2)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       sptrf_free();
       printf(" sptrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,(n*(n+1)/2));

   } /* end of Constructor  */

sptrf_dcomplex_parameters:: ~sptrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   sptrf_free();
}

TEST(sptrf,zsptrf1) {

    /* LAPACKE ZSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zsptrf) ( int matrix_layout,char uplo,
                                  lapack_int n,lapack_complex_double *a,
                                                       lapack_int *ipiv);

    Fptr_NL_LAPACKE_zsptrf ZSPTRF;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    sptrf_dcomplex_parameters   zsptrf_obj(LAPACK_ROW_MAJOR,'U',121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZSPTRF = (Fptr_NL_LAPACKE_zsptrf)dlsym(hModule,"LAPACKE_zsptrf");
    if (NULL == ZSPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zsptrf_obj.inforef = ZSPTRF( zsptrf_obj.matrix_layout,zsptrf_obj.uplo,
                                             zsptrf_obj.n,zsptrf_obj.aref,
                                                     zsptrf_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zsptrf_obj.info = LAPACKE_zsptrf(zsptrf_obj.matrix_layout,zsptrf_obj.uplo,
                                                    zsptrf_obj.n,zsptrf_obj.a,
                                                             zsptrf_obj.ipiv);

    if( zsptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsptrf is wrong\n",
                    zsptrf_obj.info );
    }
    if( zsptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsptrf is wrong\n",
        zsptrf_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zsptrf_obj.n,zsptrf_obj.ipiv,zsptrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsptrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zsptrf_obj.n*(zsptrf_obj.n + 1)/2,zsptrf_obj.a,zsptrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
