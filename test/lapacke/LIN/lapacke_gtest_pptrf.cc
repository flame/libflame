#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define pptrf_free() \
       free (ap   ); \
       free (apref)

/* Begin pptrf_double_parameters  class definition */
class pptrf_double_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      double *ap,*apref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pptrf_double_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i);
      ~pptrf_double_parameters (); 
};  /* end of pptrf_double_parameters  class definition */


/* Constructor pptrf_double_parameters definition */
pptrf_double_parameters:: pptrf_double_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &ap, &apref, (n*(n+1)/2));

    if( (ap==NULL) || (apref==NULL) ){
       pptrf_free();
       printf(" pptrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ap,apref,(n*(n+1)/2));

   } /* end of Constructor  */

pptrf_double_parameters:: ~pptrf_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptrf_double_parameters object: destructor invoked. \n");
#endif
   pptrf_free();
}

TEST(pptrf,dpptrf1) {

    /* LAPACKE DPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpptrf) ( int matrix_layout ,char uplo ,
                             lapack_int n ,double *ap );
                    
    Fptr_NL_LAPACKE_dpptrf DPPTRF;
    void *hModule,*dModule;
    double diff;
    pptrf_double_parameters   dpptrf_obj(LAPACK_ROW_MAJOR,'U',451);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DPPTRF = (Fptr_NL_LAPACKE_dpptrf)dlsym(hModule,"LAPACKE_dpptrf");
    if (NULL == DPPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }
    
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    dpptrf_obj.inforef = DPPTRF( dpptrf_obj.matrix_layout,dpptrf_obj.uplo,
                                           dpptrf_obj.n,dpptrf_obj.apref);

        /* Compute libflame's Lapacke o/p  */
    dpptrf_obj.info     = LAPACKE_dpptrf( dpptrf_obj.matrix_layout,
                                     dpptrf_obj.uplo,dpptrf_obj.n,
                                                    dpptrf_obj.ap);
    if( dpptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpptrf is wrong\n",
                    dpptrf_obj.info );
    }
    if( dpptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpptrf is wrong\n",
        dpptrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dpptrf_obj.n*(dpptrf_obj.n +1)/2,dpptrf_obj.ap,dpptrf_obj.apref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pptrf_float_parameters  class definition */
class pptrf_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      float *ap,*apref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pptrf_float_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i);
      ~pptrf_float_parameters (); 
};  /* end of pptrf_float_parameters  class definition */


/* Constructor pptrf_float_parameters definition */
pptrf_float_parameters:: pptrf_float_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &ap, &apref, (n*(n+1)/2));

    if( (ap==NULL) || (apref==NULL) ){
       pptrf_free();
       printf(" pptrf_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ap,apref,(n*(n+1)/2));

   } /* end of Constructor  */

pptrf_float_parameters:: ~pptrf_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptrf_float_parameters object: destructor invoked. \n");
#endif
   pptrf_free();
}

TEST(pptrf,spptrf1) {

    /* LAPACKE SPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spptrf) ( int matrix_layout ,char uplo ,
                             lapack_int n ,float *ap );
                    
    Fptr_NL_LAPACKE_spptrf SPPTRF;
    void *hModule,*dModule;
    float diff;
    pptrf_float_parameters   spptrf_obj(LAPACK_COL_MAJOR,'U',1020);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SPPTRF = (Fptr_NL_LAPACKE_spptrf)dlsym(hModule,"LAPACKE_spptrf");
    if (NULL == SPPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }
    
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    spptrf_obj.inforef = SPPTRF( spptrf_obj.matrix_layout,spptrf_obj.uplo,
                            spptrf_obj.n,spptrf_obj.apref);

        /* Compute libflame's Lapacke o/p  */
    spptrf_obj.info     = LAPACKE_spptrf( spptrf_obj.matrix_layout,
                         spptrf_obj.uplo,spptrf_obj.n,spptrf_obj.ap);

    if( spptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spptrf is wrong\n",
                    spptrf_obj.info );
    }
    if( spptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spptrf is wrong\n",
        spptrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( spptrf_obj.n*(spptrf_obj.n +1)/2,spptrf_obj.ap,spptrf_obj.apref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pptrf_scomplex_parameters  class definition */
class pptrf_scomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      lapack_complex_float *ap,*apref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pptrf_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i);
      ~pptrf_scomplex_parameters (); 
};  /* end of pptrf_scomplex_parameters  class definition */


/* Constructor pptrf_scomplex_parameters definition */
pptrf_scomplex_parameters:: pptrf_scomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ap, &apref, (n*(n+1)/2)); 

    if( (ap==NULL) || (apref==NULL) ){
       pptrf_free();
       printf(" pptrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ap,apref,(n*(n+1)/2));

   } /* end of Constructor  */

pptrf_scomplex_parameters:: ~pptrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptrf_scomplex_parameters object: destructor invoked. \n");
#endif
   pptrf_free();
}

TEST(pptrf,cpptrf1) {

    /* LAPACKE CPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpptrf) ( int matrix_layout ,char uplo ,
                             lapack_int n ,lapack_complex_float *ap );

    Fptr_NL_LAPACKE_cpptrf CPPTRF;
    void *hModule,*dModule;
    float diff;
    pptrf_scomplex_parameters   cpptrf_obj(LAPACK_COL_MAJOR,'U',510);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CPPTRF = (Fptr_NL_LAPACKE_cpptrf)dlsym(hModule,"LAPACKE_cpptrf");
    if (NULL == CPPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cpptrf_obj.inforef = CPPTRF( cpptrf_obj.matrix_layout,cpptrf_obj.uplo,
                            cpptrf_obj.n,cpptrf_obj.apref);

        /* Compute libflame's Lapacke o/p  */
    cpptrf_obj.info     = LAPACKE_cpptrf( cpptrf_obj.matrix_layout,
                                     cpptrf_obj.uplo,cpptrf_obj.n,
                                     cpptrf_obj.ap);
    if( cpptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpptrf is wrong\n",
                    cpptrf_obj.info );
    }
    if( cpptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpptrf is wrong\n",
        cpptrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( cpptrf_obj.n*(cpptrf_obj.n + 1)/2,cpptrf_obj.ap,
                                                          cpptrf_obj.apref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pptrf_dcomplex_parameters  class definition */
class pptrf_dcomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      lapack_complex_double *ap,*apref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pptrf_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i);
      ~pptrf_dcomplex_parameters (); 
};  /* end of pptrf_dcomplex_parameters  class definition */


/* Constructor pptrf_dcomplex_parameters definition */
pptrf_dcomplex_parameters:: pptrf_dcomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ap, &apref, (n*(n+1)/2)); 

    if( (ap==NULL) || (apref==NULL) ){
       pptrf_free();
       printf(" pptrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ap,apref,(n*(n+1)/2));

   } /* end of Constructor  */

pptrf_dcomplex_parameters:: ~pptrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   pptrf_free();
}

TEST(pptrf,zpptrf1) {

    /* LAPACKE ZPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpptrf) ( int matrix_layout ,char uplo ,
                             lapack_int n ,lapack_complex_double *ap );

    Fptr_NL_LAPACKE_zpptrf ZPPTRF;
    void *hModule,*dModule;
    double diff;
    pptrf_dcomplex_parameters   zpptrf_obj(LAPACK_ROW_MAJOR,'U',100);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZPPTRF = (Fptr_NL_LAPACKE_zpptrf)dlsym(hModule,"LAPACKE_zpptrf");
    if (NULL == ZPPTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zpptrf_obj.inforef = ZPPTRF( zpptrf_obj.matrix_layout,zpptrf_obj.uplo,
                            zpptrf_obj.n,zpptrf_obj.apref);

        /* Compute libflame's Lapacke o/p  */
    zpptrf_obj.info     = LAPACKE_zpptrf( zpptrf_obj.matrix_layout,
                                     zpptrf_obj.uplo,zpptrf_obj.n,
                                     zpptrf_obj.ap);
    if( zpptrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpptrf is wrong\n",
                    zpptrf_obj.info );
    }
    if( zpptrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpptrf is wrong\n",
        zpptrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zpptrf_obj.n*(zpptrf_obj.n +1)/2,zpptrf_obj.ap,zpptrf_obj.apref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
