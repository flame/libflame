#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define sytrf_rk_free() \
       free (a   ); \
       free (aref); \
       free (e   ); \
       free (eref); \
       free (ipiv); \
       free (ipivref)

/* Begin sytrf_rk_double_parameters  class definition */
class sytrf_rk_double_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      double *a,*aref; //The array ab contains the matrix A
 
      /* Output parameters */
      /*  superdiagonal (or subdiagonal) elements of the symmetric block diagonal 
           matrix D with 1-by-1 or 2-by-2 diagonal blocks.  */
      double *e,*eref;
      lapack_int *ipiv,*ipivref; // The ipivot indices
      /* Return Values */
      lapack_int info,inforef;

   public: 
      sytrf_rk_double_parameters ( int matrix_layout_i,char uplo_i,
                    lapack_int n_i,lapack_int lda_i);
              
      ~sytrf_rk_double_parameters (); 
};  /* end of sytrf_rk_double_parameters  class definition */


/* Constructor sytrf_rk_double_parameters definition */
sytrf_rk_double_parameters:: sytrf_rk_double_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (lda*n));
    lapacke_gtest_alloc_double_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       sytrf_rk_free();
       printf(" sytrf_rk_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_rk_double_parameters:: ~sytrf_rk_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_rk_double_parameters object: destructor invoked. \n");
#endif
   sytrf_rk_free();
}

TEST(sytrf_rk,dsytrf_rk1) {

    /* LAPACKE DSYTRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf_rk) ( int matrix_layout,char uplo,
                                 lapack_int n,double *a,lapack_int lda,
                                              double *e,lapack_int *ipiv);

    Fptr_NL_LAPACKE_dsytrf_rk DSYTRF_RK;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    sytrf_rk_double_parameters   dsytrf_rk_obj(LAPACK_ROW_MAJOR,'U',451,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DSYTRF_RK = (Fptr_NL_LAPACKE_dsytrf_rk)dlsym(hModule,"LAPACKE_dsytrf_rk");
    if (NULL == DSYTRF_RK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dsytrf_rk_obj.inforef = DSYTRF_RK( dsytrf_rk_obj.matrix_layout,dsytrf_rk_obj.uplo,
                             dsytrf_rk_obj.n,dsytrf_rk_obj.aref,dsytrf_rk_obj.lda,
                                           dsytrf_rk_obj.eref,dsytrf_rk_obj.ipivref);

    /* Compute libflame's Lapacke o/p  */
    dsytrf_rk_obj.info = LAPACKE_dsytrf_rk( dsytrf_rk_obj.matrix_layout,dsytrf_rk_obj.uplo,
                                     dsytrf_rk_obj.n,dsytrf_rk_obj.a,dsytrf_rk_obj.lda,
                                                    dsytrf_rk_obj.e,dsytrf_rk_obj.ipiv);

    if( dsytrf_rk_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrf_rk is wrong\n",
                    dsytrf_rk_obj.info );
    }
    if( dsytrf_rk_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrf_rk is wrong\n",
        dsytrf_rk_obj.inforef );
    }
    ipiv_diff = computeDiff_i( dsytrf_rk_obj.n,dsytrf_rk_obj.ipiv,dsytrf_rk_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_rk1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dsytrf_rk_obj.n*dsytrf_rk_obj.lda,dsytrf_rk_obj.a,
                                                         dsytrf_rk_obj.aref );
    diff +=  computeDiff_d( dsytrf_rk_obj.n,dsytrf_rk_obj.e,dsytrf_rk_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_rk_float_parameters  class definition */
class sytrf_rk_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      float *a,*aref; //The array ab contains the matrix A

      /* Output parameters */
      /*  superdiagonal (or subdiagonal) elements of the symmetric block diagonal 
           matrix D with 1-by-1 or 2-by-2 diagonal blocks.  */
      float *e,*eref;
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      sytrf_rk_float_parameters ( int matrix_layout_i,char uplo_i,
                  lapack_int n_i,lapack_int lda_i);
      ~sytrf_rk_float_parameters (); 
};  /* end of sytrf_rk_float_parameters  class definition */


/* Constructor sytrf_rk_float_parameters definition */
sytrf_rk_float_parameters:: sytrf_rk_float_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);
    lapacke_gtest_alloc_float_buffer_pair( &e, &eref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       sytrf_rk_free();
       printf(" sytrf_rk_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_rk_float_parameters:: ~sytrf_rk_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_rk_float_parameters object: destructor invoked. \n");
#endif
   sytrf_rk_free();
}

TEST(sytrf_rk,ssytrf_rk1) {

    /* LAPACKE SSYTRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf_rk) ( int matrix_layout,char uplo,
                                 lapack_int n,float *a,lapack_int lda,
                                             float *e,lapack_int *ipiv);

    Fptr_NL_LAPACKE_ssytrf_rk SSYTRF_RK;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;

    sytrf_rk_float_parameters   ssytrf_rk_obj(LAPACK_COL_MAJOR,'U',1020,1041);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SSYTRF_RK = (Fptr_NL_LAPACKE_ssytrf_rk)dlsym(hModule,"LAPACKE_ssytrf_rk");
    if (NULL == SSYTRF_RK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    ssytrf_rk_obj.inforef = SSYTRF_RK( ssytrf_rk_obj.matrix_layout,ssytrf_rk_obj.uplo,
                            ssytrf_rk_obj.n,ssytrf_rk_obj.aref,ssytrf_rk_obj.lda,
                                       ssytrf_rk_obj.eref,ssytrf_rk_obj.ipivref );

        /* Compute libflame's Lapacke o/p  */
    ssytrf_rk_obj.info = LAPACKE_ssytrf_rk(ssytrf_rk_obj.matrix_layout,ssytrf_rk_obj.uplo,
                                    ssytrf_rk_obj.n,ssytrf_rk_obj.a,ssytrf_rk_obj.lda,
                                                   ssytrf_rk_obj.e,ssytrf_rk_obj.ipiv);

    if( ssytrf_rk_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrf_rk is wrong\n",
                    ssytrf_rk_obj.info );
    }
    if( ssytrf_rk_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrf_rk is wrong\n",
        ssytrf_rk_obj.inforef );
    }

    ipiv_diff = computeDiff_i( ssytrf_rk_obj.n,ssytrf_rk_obj.ipiv,ssytrf_rk_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_rk1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( ssytrf_rk_obj.n*ssytrf_rk_obj.lda,ssytrf_rk_obj.a,
                                                   ssytrf_rk_obj.aref );
    diff +=  computeDiff_s( ssytrf_rk_obj.n,ssytrf_rk_obj.e,ssytrf_rk_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_rk_scomplex_parameters  class definition */
class sytrf_rk_scomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_float *a,*aref; //The array ab contains the matrix A
      /* Output parameters */
      /*  superdiagonal (or subdiagonal) elements of the symmetric block diagonal 
           matrix D with 1-by-1 or 2-by-2 diagonal blocks.  */
      lapack_complex_float *e,*eref;
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      sytrf_rk_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                 lapack_int n_i,lapack_int lda_i);
      ~sytrf_rk_scomplex_parameters (); 
};  /* end of sytrf_rk_scomplex_parameters  class definition */


/* Constructor sytrf_rk_scomplex_parameters definition */
sytrf_rk_scomplex_parameters:: sytrf_rk_scomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                            lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e, &eref, n); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       sytrf_rk_free();
       printf(" sytrf_rk_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_rk_scomplex_parameters:: ~sytrf_rk_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_rk_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrf_rk_free();
}

TEST(sytrf_rk,csytrf_rk1) {

    /* LAPACKE CSYTRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf_rk) ( int matrix_layout,char uplo,
                   lapack_int n,lapack_complex_float *a,lapack_int lda,
                              lapack_complex_float *e,lapack_int *ipiv);

    Fptr_NL_LAPACKE_csytrf_rk CSYTRF_RK;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    sytrf_rk_scomplex_parameters   csytrf_rk_obj(LAPACK_COL_MAJOR,'U',510,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CSYTRF_RK = (Fptr_NL_LAPACKE_csytrf_rk)dlsym(hModule,"LAPACKE_csytrf_rk");
    if (NULL == CSYTRF_RK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    csytrf_rk_obj.inforef = CSYTRF_RK( csytrf_rk_obj.matrix_layout,csytrf_rk_obj.uplo,
                            csytrf_rk_obj.n,csytrf_rk_obj.aref,csytrf_rk_obj.lda,
                                        csytrf_rk_obj.eref,csytrf_rk_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    csytrf_rk_obj.info = LAPACKE_csytrf_rk(csytrf_rk_obj.matrix_layout,csytrf_rk_obj.uplo,
                                    csytrf_rk_obj.n,csytrf_rk_obj.a,csytrf_rk_obj.lda,
                                                  csytrf_rk_obj.e,csytrf_rk_obj.ipiv );

    if( csytrf_rk_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrf_rk is wrong\n",
                    csytrf_rk_obj.info );
    }
    if( csytrf_rk_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrf_rk is wrong\n",
        csytrf_rk_obj.inforef );
    }

    ipiv_diff = computeDiff_i( csytrf_rk_obj.n,csytrf_rk_obj.ipiv,csytrf_rk_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_rk1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( csytrf_rk_obj.n*csytrf_rk_obj.lda,csytrf_rk_obj.a,
                                                         csytrf_rk_obj.aref );
    diff +=  computeDiff_c( csytrf_rk_obj.n,csytrf_rk_obj.e,csytrf_rk_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrf_rk_dcomplex_parameters  class definition */
class sytrf_rk_dcomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_double *a,*aref; //The array ab contains the matrix A
      /* Output parameters */
      /*  superdiagonal (or subdiagonal) elements of the symmetric block diagonal 
           matrix D with 1-by-1 or 2-by-2 diagonal blocks.  */
      lapack_complex_double *e,*eref;
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      sytrf_rk_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~sytrf_rk_dcomplex_parameters (); 
};  /* end of sytrf_rk_dcomplex_parameters  class definition */


/* Constructor sytrf_rk_dcomplex_parameters definition */
sytrf_rk_dcomplex_parameters:: sytrf_rk_dcomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e, &eref, n); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       sytrf_rk_free();
       printf(" sytrf_rk_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

sytrf_rk_dcomplex_parameters:: ~sytrf_rk_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrf_rk_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrf_rk_free();
}

TEST(sytrf_rk,zsytrf_rk1) {

    /* LAPACKE ZSYTRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf_rk) ( int matrix_layout,char uplo,
                  lapack_int n,lapack_complex_double *a,lapack_int lda,
                             lapack_complex_double *e,lapack_int *ipiv);

    Fptr_NL_LAPACKE_zsytrf_rk ZSYTRF_RK;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    sytrf_rk_dcomplex_parameters   zsytrf_rk_obj(LAPACK_ROW_MAJOR,'U',100,121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZSYTRF_RK = (Fptr_NL_LAPACKE_zsytrf_rk)dlsym(hModule,"LAPACKE_zsytrf_rk");
    if (NULL == ZSYTRF_RK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zsytrf_rk_obj.inforef = ZSYTRF_RK( zsytrf_rk_obj.matrix_layout,zsytrf_rk_obj.uplo,
                            zsytrf_rk_obj.n,zsytrf_rk_obj.aref,zsytrf_rk_obj.lda,
                                        zsytrf_rk_obj.eref,zsytrf_rk_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zsytrf_rk_obj.info = LAPACKE_zsytrf_rk(zsytrf_rk_obj.matrix_layout,zsytrf_rk_obj.uplo,
                                   zsytrf_rk_obj.n,zsytrf_rk_obj.a,zsytrf_rk_obj.lda,
                                                  zsytrf_rk_obj.e,zsytrf_rk_obj.ipiv);

    if( zsytrf_rk_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrf_rk is wrong\n",
                    zsytrf_rk_obj.info );
    }
    if( zsytrf_rk_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrf_rk is wrong\n",
        zsytrf_rk_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zsytrf_rk_obj.n,zsytrf_rk_obj.ipiv,zsytrf_rk_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dsytrf_rk1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff  =  computeDiff_z( zsytrf_rk_obj.n*zsytrf_rk_obj.lda,zsytrf_rk_obj.a,
                                                          zsytrf_rk_obj.aref );
    diff +=  computeDiff_z( zsytrf_rk_obj.n,zsytrf_rk_obj.e,zsytrf_rk_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
