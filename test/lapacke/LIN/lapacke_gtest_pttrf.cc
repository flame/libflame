#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define pttrf_free() \
       free (d); \
       free (dref); \
       free (e   ); \
       free (eref)


/* Begin pttrf_double_parameters  class definition */
class pttrf_double_parameters{

   public:
      /* Input parameters */
      lapack_int n; // No of rows,Columns
    
      /* Input/ Output parameters */
      double *d,*dref; // diagonal elements of A
      double *e,*eref; // subdiagonal elements of A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pttrf_double_parameters ( lapack_int n_i);
      ~pttrf_double_parameters (); 
};  /* end of pttrf_double_parameters  class definition */

/* Constructor pttrf_double_parameters definition */
pttrf_double_parameters:: pttrf_double_parameters (lapack_int n_i) {

    n = n_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_double_buffer_pair( &e, &eref, n-1);

    if((d==NULL) || (dref==NULL) ||  \
       (e==NULL) || (eref==NULL)){
       pttrf_free();
       printf(" pttrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( d,dref,n);
    lapacke_gtest_init_double_buffer_pair_rand( e,eref,n-1);

   } /* end of Constructor  */

pttrf_double_parameters:: ~pttrf_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pttrf_double_parameters object: destructor invoked. \n");
#endif
   pttrf_free();
}

TEST(pttrf,dpttrf1) {

    /* LAPACKE dpttrf prototype */
    typedef int (*Fptr_NL_LAPACKE_dpttrf) (lapack_int n,double *d,double *e);

    Fptr_NL_LAPACKE_dpttrf DPTTRF;
    double diff;
    void *hModule,*dModule;

    pttrf_double_parameters dpttrf_obj(140 );

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    DPTTRF = (Fptr_NL_LAPACKE_dpttrf)dlsym(hModule,"LAPACKE_dpttrf");
    if (NULL == DPTTRF)
    {
        printf("Could not get the symbol LAPACKE_dpttrf. Exiting...\n");
        dlclose(hModule);
        dlclose(dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    dpttrf_obj.info     = LAPACKE_dpttrf( dpttrf_obj.n,dpttrf_obj.d,
                                                       dpttrf_obj.e);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dpttrf_obj.inforef = DPTTRF( dpttrf_obj.n,dpttrf_obj.dref,
                                              dpttrf_obj.eref);

    if( dpttrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpttrf \
                  is wrong\n",dpttrf_obj.info );
    }
    if( dpttrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpttrf is wrong\n",
        dpttrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dpttrf_obj.n,dpttrf_obj.d,dpttrf_obj.dref );
    diff +=  computeDiff_d( dpttrf_obj.n-1,dpttrf_obj.e,dpttrf_obj.e );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pttrf_float_parameters  class definition */
class pttrf_float_parameters{

   public:
      /* Input parameters */
      lapack_int n; // No of rows,Columns
    
      /* Input/ Output parameters */
      float *e,*eref;// subdiagonal elements of A
      float *d,*dref; // diagonal elements of A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pttrf_float_parameters ( lapack_int n_i);
      ~pttrf_float_parameters (); 
};  /* end of pttrf_float_parameters  class definition */

/* Constructor pttrf_float_parameters definition */
pttrf_float_parameters:: pttrf_float_parameters (lapack_int n_i) {

    n = n_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_float_buffer_pair( &e, &eref, n-1);

    if((d==NULL) || (dref==NULL) ||  \
       (e==NULL) || (e==NULL)){
       pttrf_free();
       printf(" pttrf_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( d,dref,n);
    lapacke_gtest_init_float_buffer_pair_rand( e,eref,n-1);

   } /* end of Constructor  */

pttrf_float_parameters:: ~pttrf_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pttrf_float_parameters object: destructor invoked. \n");
#endif
   pttrf_free();
}

TEST(pttrf,spttrf1) {

    /* LAPACKE spttrf prototype */
    typedef int (*Fptr_NL_LAPACKE_spttrf) ( lapack_int n,float *d,
                                                          float *e );

    Fptr_NL_LAPACKE_spttrf SPTTRF;
    float diff;
    void *hModule,*dModule;

    pttrf_float_parameters spttrf_obj(453 );

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    SPTTRF = (Fptr_NL_LAPACKE_spttrf)dlsym(hModule,"LAPACKE_spttrf");
    if (NULL == SPTTRF)
    {
        printf("Could not get the symbol LAPACKE_spttrf. Exiting...\n");
        dlclose(hModule);
        dlclose(dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    spttrf_obj.info     = LAPACKE_spttrf( spttrf_obj.n,spttrf_obj.d,
                                                       spttrf_obj.e);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    spttrf_obj.inforef = SPTTRF( spttrf_obj.n,spttrf_obj.dref,
                                              spttrf_obj.eref);

    if( spttrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spttrf \
                  is wrong\n",spttrf_obj.info );
    }
    if( spttrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spttrf is wrong\n",
        spttrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff  =  computeDiff_s( spttrf_obj.n,spttrf_obj.d,spttrf_obj.dref );
    diff +=  computeDiff_s( spttrf_obj.n-1,spttrf_obj.e,spttrf_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pttrf_scomplex_parameters  class definition */
class pttrf_scomplex_parameters{

   public:
      /* Input parameters */
      lapack_int n; // No of rows,Columns
    
      /* Input/ Output parameters */
      lapack_complex_float *e,*eref;// subdiagonal elements of A
      float *d,*dref; // diagonal elements of A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pttrf_scomplex_parameters ( lapack_int n_i);
      ~pttrf_scomplex_parameters (); 
};  /* end of pttrf_scomplex_parameters  class definition */

/* Constructor pttrf_scomplex_parameters definition */
pttrf_scomplex_parameters:: pttrf_scomplex_parameters (lapack_int n_i) {

    n = n_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e, &eref, n-1);
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n);

    if((d==NULL) || (dref==NULL) ||  \
       (e==NULL) || (eref==NULL)){
       pttrf_free();
       printf(" pttrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( d,dref,n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( e,eref,n-1);

   } /* end of Constructor  */

pttrf_scomplex_parameters:: ~pttrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pttrf_scomplex_parameters object: destructor invoked. \n");
#endif
   pttrf_free();
}

TEST(pttrf,cpttrf1) {

    /* LAPACKE cpttrf prototype */
    typedef int (*Fptr_NL_LAPACKE_cpttrf) ( lapack_int n,float *d,
                                                   lapack_complex_float *e);

    Fptr_NL_LAPACKE_cpttrf CPTTRF;
    float diff;
    void *hModule,*dModule;

    pttrf_scomplex_parameters cpttrf_obj(350 );

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    CPTTRF = (Fptr_NL_LAPACKE_cpttrf)dlsym(hModule,"LAPACKE_cpttrf");
    if (NULL == CPTTRF)
    {
        printf("Could not get the symbol LAPACKE_cpttrf. Exiting...\n");
        dlclose(hModule);
        dlclose(dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    cpttrf_obj.info     = LAPACKE_cpttrf( cpttrf_obj.n,cpttrf_obj.d,
                                                      cpttrf_obj.e);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cpttrf_obj.inforef = CPTTRF( cpttrf_obj.n,cpttrf_obj.dref,
                                                  cpttrf_obj.eref);

    if( cpttrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpttrf \
                  is wrong\n",cpttrf_obj.info );
    }
    if( cpttrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpttrf is wrong\n",
        cpttrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff  =  computeDiff_s( cpttrf_obj.n,cpttrf_obj.d,cpttrf_obj.dref );
    diff +=  computeDiff_c( cpttrf_obj.n-1,cpttrf_obj.e,cpttrf_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pttrf_dcomplex_parameters  class definition */
class pttrf_dcomplex_parameters{

   public:
      /* Input parameters */
      lapack_int n; // No of rows,Columns
    
      /* Input/ Output parameters */
      lapack_complex_double *e,*eref;// subdiagonal elements of A
      double *d,*dref; // diagonal elements of A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pttrf_dcomplex_parameters ( lapack_int n_i);
      ~pttrf_dcomplex_parameters (); 
};  /* end of pttrf_dcomplex_parameters  class definition */

/* Constructor pttrf_dcomplex_parameters definition */
pttrf_dcomplex_parameters:: pttrf_dcomplex_parameters (lapack_int n_i) {

    n = n_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e, &eref, n-1);

    if((d==NULL) || (dref==NULL) ||  \
       (e==NULL) || (eref==NULL)){
       pttrf_free();
       printf(" pttrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( d,dref,n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( e,eref,n-1);

   } /* end of Constructor  */

pttrf_dcomplex_parameters:: ~pttrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pttrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   pttrf_free();
}

TEST(pttrf,zpttrf1) {

    /* LAPACKE zpttrf prototype */
    typedef int (*Fptr_NL_LAPACKE_zpttrf) ( lapack_int n,
                 double *d,lapack_complex_double *e);

    Fptr_NL_LAPACKE_zpttrf ZPTTRF;
    double diff;
    void *hModule,*dModule;

    pttrf_dcomplex_parameters zpttrf_obj(415 );

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    ZPTTRF = (Fptr_NL_LAPACKE_zpttrf)dlsym(hModule,"LAPACKE_zpttrf");
    if (NULL == ZPTTRF)
    {
        printf("Could not get the symbol LAPACKE_zpttrf. Exiting...\n");
        dlclose(hModule);
        dlclose(dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    zpttrf_obj.info     = LAPACKE_zpttrf( zpttrf_obj.n,zpttrf_obj.d,
                           zpttrf_obj.e);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zpttrf_obj.inforef = ZPTTRF( zpttrf_obj.n,zpttrf_obj.dref,
                           zpttrf_obj.eref);

    if( zpttrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpttrf \
                  is wrong\n",zpttrf_obj.info );
    }
    if( zpttrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpttrf is wrong\n",
        zpttrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff  =  computeDiff_d( zpttrf_obj.n,zpttrf_obj.d,zpttrf_obj.dref );
    diff +=  computeDiff_z( zpttrf_obj.n-1,zpttrf_obj.e,zpttrf_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}