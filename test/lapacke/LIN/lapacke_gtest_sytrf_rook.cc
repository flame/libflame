#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define sytrf_rook_free() \
       free (a   ); \
       free (aref); \
       free (ipiv); \
       free (ipivref)

/* Begin sytrf_rook_double_parameters  class definition */
class sytrf_rook_double_parameters{

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
      sytrf_rook_double_parameters ( int matrix_layout_i,char uplo_i,
                    lapack_int n_i,lapack_int lda_i);
              
      ~sytrf_rook_double_parameters (); 
};  /* end of sytrf_rook_double_parameters  class definition */


/* Constructor sytrf_rook_double_parameters definition */
sytrf_rook_double_parameters:: sytrf_rook_double_parameters ( int matrix_layout_i,
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
       sytrf_rook_free();
       printf(" sytrf_rook_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_rook_double_parameters:: ~sytrf_rook_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_rook_double_parameters object: destructor invoked. \n");
#endif
   sytrf_rook_free();
}

TEST(sytrf_rook,dsytrf_rook1) {

    /* LAPACKE DSYTRF_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf_rook) ( int matrix_layout,char uplo,
                                 lapack_int n,double *a,lapack_int lda,
                                                       lapack_int *ipiv);

    Fptr_NL_LAPACKE_dsytrf_rook DSYTRF_ROOK;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    sytrf_rook_double_parameters   dsytrf_rook_obj(LAPACK_ROW_MAJOR,'U',451,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DSYTRF_ROOK = (Fptr_NL_LAPACKE_dsytrf_rook)dlsym(hModule,"LAPACKE_dsytrf_rook");
    if (NULL == DSYTRF_ROOK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dsytrf_rook_obj.inforef = DSYTRF_ROOK( dsytrf_rook_obj.matrix_layout,dsytrf_rook_obj.uplo,
                             dsytrf_rook_obj.n,dsytrf_rook_obj.aref,dsytrf_rook_obj.lda,
                                                       dsytrf_rook_obj.ipivref);

    /* Compute libflame's Lapacke o/p  */
    dsytrf_rook_obj.info = LAPACKE_dsytrf_rook( dsytrf_rook_obj.matrix_layout,dsytrf_rook_obj.uplo,
                                     dsytrf_rook_obj.n,dsytrf_rook_obj.a,dsytrf_rook_obj.lda,
                                                               dsytrf_rook_obj.ipiv);

    if( dsytrf_rook_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrf_rook is wrong\n",
                    dsytrf_rook_obj.info );
    }
    if( dsytrf_rook_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrf_rook is wrong\n",
        dsytrf_rook_obj.inforef );
    }
    ipiv_diff = computeDiff_i( dsytrf_rook_obj.n,dsytrf_rook_obj.ipiv,dsytrf_rook_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_rook1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dsytrf_rook_obj.n*dsytrf_rook_obj.lda,dsytrf_rook_obj.a,dsytrf_rook_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_rook_float_parameters  class definition */
class sytrf_rook_float_parameters{

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
      sytrf_rook_float_parameters ( int matrix_layout_i,char uplo_i,
                  lapack_int n_i,lapack_int lda_i);
      ~sytrf_rook_float_parameters (); 
};  /* end of sytrf_rook_float_parameters  class definition */


/* Constructor sytrf_rook_float_parameters definition */
sytrf_rook_float_parameters:: sytrf_rook_float_parameters ( int matrix_layout_i,
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
       sytrf_rook_free();
       printf(" sytrf_rook_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_rook_float_parameters:: ~sytrf_rook_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_rook_float_parameters object: destructor invoked. \n");
#endif
   sytrf_rook_free();
}

TEST(sytrf_rook,ssytrf_rook1) {

    /* LAPACKE SSYTRF_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf_rook) ( int matrix_layout ,char uplo,
                                 lapack_int n ,float *a ,lapack_int lda,
                                                        lapack_int *ipiv);

    Fptr_NL_LAPACKE_ssytrf_rook SSYTRF_ROOK;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;

    sytrf_rook_float_parameters   ssytrf_rook_obj(LAPACK_COL_MAJOR,'U',1020,1041);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SSYTRF_ROOK = (Fptr_NL_LAPACKE_ssytrf_rook)dlsym(hModule,"LAPACKE_ssytrf_rook");
    if (NULL == SSYTRF_ROOK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    ssytrf_rook_obj.inforef = SSYTRF_ROOK( ssytrf_rook_obj.matrix_layout,ssytrf_rook_obj.uplo,
                            ssytrf_rook_obj.n,ssytrf_rook_obj.aref,ssytrf_rook_obj.lda,
                                                     ssytrf_rook_obj.ipivref );

        /* Compute libflame's Lapacke o/p  */
    ssytrf_rook_obj.info = LAPACKE_ssytrf_rook(ssytrf_rook_obj.matrix_layout,ssytrf_rook_obj.uplo,
                                    ssytrf_rook_obj.n,ssytrf_rook_obj.a,ssytrf_rook_obj.lda,
                                                              ssytrf_rook_obj.ipiv);

    if( ssytrf_rook_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrf_rook is wrong\n",
                    ssytrf_rook_obj.info );
    }
    if( ssytrf_rook_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrf_rook is wrong\n",
        ssytrf_rook_obj.inforef );
    }

    ipiv_diff = computeDiff_i( ssytrf_rook_obj.n,ssytrf_rook_obj.ipiv,ssytrf_rook_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_rook1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( ssytrf_rook_obj.n*ssytrf_rook_obj.lda,ssytrf_rook_obj.a,
                                                   ssytrf_rook_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_rook_scomplex_parameters  class definition */
class sytrf_rook_scomplex_parameters{

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
      sytrf_rook_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                 lapack_int n_i,lapack_int lda_i);
      ~sytrf_rook_scomplex_parameters (); 
};  /* end of sytrf_rook_scomplex_parameters  class definition */


/* Constructor sytrf_rook_scomplex_parameters definition */
sytrf_rook_scomplex_parameters:: sytrf_rook_scomplex_parameters ( int matrix_layout_i,
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
       sytrf_rook_free();
       printf(" sytrf_rook_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_rook_scomplex_parameters:: ~sytrf_rook_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_rook_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrf_rook_free();
}

TEST(sytrf_rook,csytrf_rook1) {

    /* LAPACKE CSYTRF_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf_rook) ( int matrix_layout ,char uplo ,
                   lapack_int n ,lapack_complex_float *a ,lapack_int lda,
                                                      lapack_int *ipiv);

    Fptr_NL_LAPACKE_csytrf_rook CSYTRF_ROOK;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    sytrf_rook_scomplex_parameters   csytrf_rook_obj(LAPACK_COL_MAJOR,'U',510,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CSYTRF_ROOK = (Fptr_NL_LAPACKE_csytrf_rook)dlsym(hModule,"LAPACKE_csytrf_rook");
    if (NULL == CSYTRF_ROOK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    csytrf_rook_obj.inforef = CSYTRF_ROOK( csytrf_rook_obj.matrix_layout,csytrf_rook_obj.uplo,
                            csytrf_rook_obj.n,csytrf_rook_obj.aref,csytrf_rook_obj.lda,
                                                      csytrf_rook_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    csytrf_rook_obj.info = LAPACKE_csytrf_rook(csytrf_rook_obj.matrix_layout,csytrf_rook_obj.uplo,
                                    csytrf_rook_obj.n,csytrf_rook_obj.a,csytrf_rook_obj.lda,
                                                             csytrf_rook_obj.ipiv );

    if( csytrf_rook_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrf_rook is wrong\n",
                    csytrf_rook_obj.info );
    }
    if( csytrf_rook_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrf_rook is wrong\n",
        csytrf_rook_obj.inforef );
    }

    ipiv_diff = computeDiff_i( csytrf_rook_obj.n,csytrf_rook_obj.ipiv,csytrf_rook_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_rook1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( csytrf_rook_obj.n*csytrf_rook_obj.lda,csytrf_rook_obj.a,csytrf_rook_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_rook_dcomplex_parameters  class definition */
class sytrf_rook_dcomplex_parameters{

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
      sytrf_rook_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~sytrf_rook_dcomplex_parameters (); 
};  /* end of sytrf_rook_dcomplex_parameters  class definition */


/* Constructor sytrf_rook_dcomplex_parameters definition */
sytrf_rook_dcomplex_parameters:: sytrf_rook_dcomplex_parameters ( int matrix_layout_i,
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
       sytrf_rook_free();
       printf(" sytrf_rook_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_rook_dcomplex_parameters:: ~sytrf_rook_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_rook_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrf_rook_free();
}

TEST(sytrf_rook,zsytrf_rook1) {

    /* LAPACKE ZSYTRF_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf_rook) ( int matrix_layout ,char uplo ,
                  lapack_int n ,lapack_complex_double *a ,lapack_int lda,
                        lapack_int *ipiv);

    Fptr_NL_LAPACKE_zsytrf_rook ZSYTRF_ROOK;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    sytrf_rook_dcomplex_parameters   zsytrf_rook_obj(LAPACK_ROW_MAJOR,'U',100,121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZSYTRF_ROOK = (Fptr_NL_LAPACKE_zsytrf_rook)dlsym(hModule,"LAPACKE_zsytrf_rook");
    if (NULL == ZSYTRF_ROOK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zsytrf_rook_obj.inforef = ZSYTRF_ROOK( zsytrf_rook_obj.matrix_layout,zsytrf_rook_obj.uplo,
                            zsytrf_rook_obj.n,zsytrf_rook_obj.aref,zsytrf_rook_obj.lda,
                                                     zsytrf_rook_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zsytrf_rook_obj.info = LAPACKE_zsytrf_rook(zsytrf_rook_obj.matrix_layout,zsytrf_rook_obj.uplo,
                                   zsytrf_rook_obj.n,zsytrf_rook_obj.a,zsytrf_rook_obj.lda,
                                                             zsytrf_rook_obj.ipiv);

    if( zsytrf_rook_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrf_rook is wrong\n",
                    zsytrf_rook_obj.info );
    }
    if( zsytrf_rook_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrf_rook is wrong\n",
        zsytrf_rook_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zsytrf_rook_obj.n,zsytrf_rook_obj.ipiv,zsytrf_rook_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_rook1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zsytrf_rook_obj.n*zsytrf_rook_obj.lda,zsytrf_rook_obj.a,zsytrf_rook_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
