#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define sytrf_aa_free() \
       free (a   ); \
       free (aref); \
       free (ipiv); \
       free (ipivref)

/* Begin sytrf_aa_double_parameters  class definition */
class sytrf_aa_double_parameters{

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
      sytrf_aa_double_parameters ( int matrix_layout_i,char uplo_i,
                    lapack_int n_i,lapack_int lda_i);
              
      ~sytrf_aa_double_parameters (); 
};  /* end of sytrf_aa_double_parameters  class definition */


/* Constructor sytrf_aa_double_parameters definition */
sytrf_aa_double_parameters:: sytrf_aa_double_parameters ( int matrix_layout_i,
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
       sytrf_aa_free();
       printf(" sytrf_aa_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_aa_double_parameters:: ~sytrf_aa_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_aa_double_parameters object: destructor invoked. \n");
#endif
   sytrf_aa_free();
}

TEST(sytrf_aa,dsytrf_aa1) {

    /* LAPACKE DSYTRF_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf_aa) ( int matrix_layout,char uplo,
                                 lapack_int n,double *a,lapack_int lda,
                                                       lapack_int *ipiv);

    Fptr_NL_LAPACKE_dsytrf_aa DSYTRF_AA;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    sytrf_aa_double_parameters   dsytrf_aa_obj(LAPACK_ROW_MAJOR,'U',451,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DSYTRF_AA = (Fptr_NL_LAPACKE_dsytrf_aa)dlsym(hModule,"LAPACKE_dsytrf_aa");
    if (NULL == DSYTRF_AA)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dsytrf_aa_obj.inforef = DSYTRF_AA( dsytrf_aa_obj.matrix_layout,dsytrf_aa_obj.uplo,
                             dsytrf_aa_obj.n,dsytrf_aa_obj.aref,dsytrf_aa_obj.lda,
                                                       dsytrf_aa_obj.ipivref);

    /* Compute libflame's Lapacke o/p  */
    dsytrf_aa_obj.info = LAPACKE_dsytrf_aa( dsytrf_aa_obj.matrix_layout,dsytrf_aa_obj.uplo,
                                     dsytrf_aa_obj.n,dsytrf_aa_obj.a,dsytrf_aa_obj.lda,
                                                               dsytrf_aa_obj.ipiv);

    if( dsytrf_aa_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrf_aa is wrong\n",
                    dsytrf_aa_obj.info );
    }
    if( dsytrf_aa_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrf_aa is wrong\n",
        dsytrf_aa_obj.inforef );
    }
    ipiv_diff = computeDiff_i( dsytrf_aa_obj.n,dsytrf_aa_obj.ipiv,dsytrf_aa_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_aa1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dsytrf_aa_obj.n*dsytrf_aa_obj.lda,dsytrf_aa_obj.a,dsytrf_aa_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_aa_float_parameters  class definition */
class sytrf_aa_float_parameters{

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
      sytrf_aa_float_parameters ( int matrix_layout_i,char uplo_i,
                  lapack_int n_i,lapack_int lda_i);
      ~sytrf_aa_float_parameters (); 
};  /* end of sytrf_aa_float_parameters  class definition */


/* Constructor sytrf_aa_float_parameters definition */
sytrf_aa_float_parameters:: sytrf_aa_float_parameters ( int matrix_layout_i,
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
       sytrf_aa_free();
       printf(" sytrf_aa_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_aa_float_parameters:: ~sytrf_aa_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_aa_float_parameters object: destructor invoked. \n");
#endif
   sytrf_aa_free();
}

TEST(sytrf_aa,ssytrf_aa1) {

    /* LAPACKE SSYTRF_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf_aa) ( int matrix_layout ,char uplo,
                                 lapack_int n ,float *a ,lapack_int lda,
                                                        lapack_int *ipiv);

    Fptr_NL_LAPACKE_ssytrf_aa SSYTRF_AA;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;

    sytrf_aa_float_parameters   ssytrf_aa_obj(LAPACK_COL_MAJOR,'U',1020,1041);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SSYTRF_AA = (Fptr_NL_LAPACKE_ssytrf_aa)dlsym(hModule,"LAPACKE_ssytrf_aa");
    if (NULL == SSYTRF_AA)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    ssytrf_aa_obj.inforef = SSYTRF_AA( ssytrf_aa_obj.matrix_layout,ssytrf_aa_obj.uplo,
                            ssytrf_aa_obj.n,ssytrf_aa_obj.aref,ssytrf_aa_obj.lda,
                                                     ssytrf_aa_obj.ipivref );

        /* Compute libflame's Lapacke o/p  */
    ssytrf_aa_obj.info = LAPACKE_ssytrf_aa(ssytrf_aa_obj.matrix_layout,ssytrf_aa_obj.uplo,
                                    ssytrf_aa_obj.n,ssytrf_aa_obj.a,ssytrf_aa_obj.lda,
                                                              ssytrf_aa_obj.ipiv);

    if( ssytrf_aa_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrf_aa is wrong\n",
                    ssytrf_aa_obj.info );
    }
    if( ssytrf_aa_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrf_aa is wrong\n",
        ssytrf_aa_obj.inforef );
    }

    ipiv_diff = computeDiff_i( ssytrf_aa_obj.n,ssytrf_aa_obj.ipiv,ssytrf_aa_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_aa1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( ssytrf_aa_obj.n*ssytrf_aa_obj.lda,ssytrf_aa_obj.a,
                                                   ssytrf_aa_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_aa_scomplex_parameters  class definition */
class sytrf_aa_scomplex_parameters{

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
      sytrf_aa_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                 lapack_int n_i,lapack_int lda_i);
      ~sytrf_aa_scomplex_parameters (); 
};  /* end of sytrf_aa_scomplex_parameters  class definition */


/* Constructor sytrf_aa_scomplex_parameters definition */
sytrf_aa_scomplex_parameters:: sytrf_aa_scomplex_parameters ( int matrix_layout_i,
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
       sytrf_aa_free();
       printf(" sytrf_aa_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_aa_scomplex_parameters:: ~sytrf_aa_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_aa_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrf_aa_free();
}

TEST(sytrf_aa,csytrf_aa1) {

    /* LAPACKE CSYTRF_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf_aa) ( int matrix_layout ,char uplo ,
                   lapack_int n ,lapack_complex_float *a ,lapack_int lda,
                                                      lapack_int *ipiv);

    Fptr_NL_LAPACKE_csytrf_aa CSYTRF_AA;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    sytrf_aa_scomplex_parameters   csytrf_aa_obj(LAPACK_COL_MAJOR,'U',510,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CSYTRF_AA = (Fptr_NL_LAPACKE_csytrf_aa)dlsym(hModule,"LAPACKE_csytrf_aa");
    if (NULL == CSYTRF_AA)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    csytrf_aa_obj.inforef = CSYTRF_AA( csytrf_aa_obj.matrix_layout,csytrf_aa_obj.uplo,
                            csytrf_aa_obj.n,csytrf_aa_obj.aref,csytrf_aa_obj.lda,
                                                      csytrf_aa_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    csytrf_aa_obj.info = LAPACKE_csytrf_aa(csytrf_aa_obj.matrix_layout,csytrf_aa_obj.uplo,
                                    csytrf_aa_obj.n,csytrf_aa_obj.a,csytrf_aa_obj.lda,
                                                             csytrf_aa_obj.ipiv );

    if( csytrf_aa_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrf_aa is wrong\n",
                    csytrf_aa_obj.info );
    }
    if( csytrf_aa_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrf_aa is wrong\n",
        csytrf_aa_obj.inforef );
    }

    ipiv_diff = computeDiff_i( csytrf_aa_obj.n,csytrf_aa_obj.ipiv,csytrf_aa_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_aa1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( csytrf_aa_obj.n*csytrf_aa_obj.lda,csytrf_aa_obj.a,csytrf_aa_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_aa_dcomplex_parameters  class definition */
class sytrf_aa_dcomplex_parameters{

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
      sytrf_aa_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~sytrf_aa_dcomplex_parameters (); 
};  /* end of sytrf_aa_dcomplex_parameters  class definition */


/* Constructor sytrf_aa_dcomplex_parameters definition */
sytrf_aa_dcomplex_parameters:: sytrf_aa_dcomplex_parameters ( int matrix_layout_i,
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
       sytrf_aa_free();
       printf(" sytrf_aa_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_aa_dcomplex_parameters:: ~sytrf_aa_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_aa_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrf_aa_free();
}

TEST(sytrf_aa,zsytrf_aa1) {

    /* LAPACKE ZSYTRF_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf_aa) ( int matrix_layout ,char uplo ,
                  lapack_int n ,lapack_complex_double *a ,lapack_int lda,
                        lapack_int *ipiv);

    Fptr_NL_LAPACKE_zsytrf_aa ZSYTRF_AA;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    sytrf_aa_dcomplex_parameters   zsytrf_aa_obj(LAPACK_ROW_MAJOR,'U',100,121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZSYTRF_AA = (Fptr_NL_LAPACKE_zsytrf_aa)dlsym(hModule,"LAPACKE_zsytrf_aa");
    if (NULL == ZSYTRF_AA)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zsytrf_aa_obj.inforef = ZSYTRF_AA( zsytrf_aa_obj.matrix_layout,zsytrf_aa_obj.uplo,
                            zsytrf_aa_obj.n,zsytrf_aa_obj.aref,zsytrf_aa_obj.lda,
                                                     zsytrf_aa_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zsytrf_aa_obj.info = LAPACKE_zsytrf_aa(zsytrf_aa_obj.matrix_layout,zsytrf_aa_obj.uplo,
                                   zsytrf_aa_obj.n,zsytrf_aa_obj.a,zsytrf_aa_obj.lda,
                                                             zsytrf_aa_obj.ipiv);

    if( zsytrf_aa_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrf_aa is wrong\n",
                    zsytrf_aa_obj.info );
    }
    if( zsytrf_aa_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrf_aa is wrong\n",
        zsytrf_aa_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zsytrf_aa_obj.n,zsytrf_aa_obj.ipiv,zsytrf_aa_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_aa1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zsytrf_aa_obj.n*zsytrf_aa_obj.lda,zsytrf_aa_obj.a,zsytrf_aa_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
