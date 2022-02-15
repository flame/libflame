#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define pftrf_free() \
       free (a   ); \
       free (aref)

/* Begin pftrf_double_parameters  class definition */
class pftrf_double_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char transr; // Must be 'N','T' (for real data) or 'C' (for complex data)
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      double *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pftrf_double_parameters ( int matrix_layout_i,char uplo_i,
                                    char transr_i,lapack_int n_i);
      ~pftrf_double_parameters (); 
};  /* end of pftrf_double_parameters  class definition */


/* Constructor pftrf_double_parameters definition */
pftrf_double_parameters:: pftrf_double_parameters ( int matrix_layout_i,
                              char uplo_i,char transr_i,lapack_int n_i){
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    transr = transr_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*(n+1)/2));

    if( (a==NULL) || (aref==NULL) ){
       pftrf_free();
       printf(" pftrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,(n*(n+1)/2));

   } /* end of Constructor  */

pftrf_double_parameters:: ~pftrf_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftrf_double_parameters object: destructor invoked. \n");
#endif
   pftrf_free();
}

TEST(pftrf,dpftrf1) {

    /* LAPACKE DPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpftrf) ( int matrix_layout,
             char transr,char uplo,lapack_int n,double *a  );
                    
    Fptr_NL_LAPACKE_dpftrf DPFTRF;
    void *hModule,*dModule;
    double diff;
    pftrf_double_parameters   dpftrf_obj(LAPACK_ROW_MAJOR,'U','N',521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DPFTRF = (Fptr_NL_LAPACKE_dpftrf)dlsym(hModule,"LAPACKE_dpftrf");
    if (NULL == DPFTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }
    
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    dpftrf_obj.inforef = DPFTRF( dpftrf_obj.matrix_layout,dpftrf_obj.transr,
                            dpftrf_obj.uplo,dpftrf_obj.n,dpftrf_obj.aref);

        /* Compute libflame's Lapacke o/p  */
    dpftrf_obj.info     = LAPACKE_dpftrf( dpftrf_obj.matrix_layout,
                                 dpftrf_obj.transr,dpftrf_obj.uplo,
                         dpftrf_obj.n,dpftrf_obj.a);
    if( dpftrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpftrf is wrong\n",
                    dpftrf_obj.info );
    }
    if( dpftrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpftrf is wrong\n",
        dpftrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dpftrf_obj.n*(dpftrf_obj.n+1)/2,dpftrf_obj.a,dpftrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pftrf_float_parameters  class definition */
class pftrf_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char transr; // Must be 'N','T' (for real data) or 'C' (for complex data)
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      float *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pftrf_float_parameters ( int matrix_layout_i,char uplo_i,
                                    char transr_i,lapack_int n_i);
      ~pftrf_float_parameters (); 
};  /* end of pftrf_float_parameters  class definition */


/* Constructor pftrf_float_parameters definition */
pftrf_float_parameters:: pftrf_float_parameters ( int matrix_layout_i,
                            char uplo_i,char transr_i,lapack_int n_i){
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    transr = transr_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*(n+1)/2));

    if( (a==NULL) || (aref==NULL) ){
       pftrf_free();
       printf(" pftrf_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,(n*(n+1)/2));

   } /* end of Constructor  */

pftrf_float_parameters:: ~pftrf_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftrf_float_parameters object: destructor invoked. \n");
#endif
   pftrf_free();
}

TEST(pftrf,spftrf1) {

    /* LAPACKE SPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spftrf) ( int matrix_layout,
              char transr,char uplo,lapack_int n,float *a  );
                    
    Fptr_NL_LAPACKE_spftrf SPFTRF;
    void *hModule,*dModule;
    float diff;
    pftrf_float_parameters   spftrf_obj(LAPACK_COL_MAJOR,'U','T',1041);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SPFTRF = (Fptr_NL_LAPACKE_spftrf)dlsym(hModule,"LAPACKE_spftrf");
    if (NULL == SPFTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }
    
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    spftrf_obj.inforef = SPFTRF( spftrf_obj.matrix_layout,spftrf_obj.transr,
                            spftrf_obj.uplo,spftrf_obj.n,spftrf_obj.aref);

        /* Compute libflame's Lapacke o/p  */
    spftrf_obj.info     = LAPACKE_spftrf( spftrf_obj.matrix_layout,
                                spftrf_obj.transr,spftrf_obj.uplo,
                                       spftrf_obj.n,spftrf_obj.a);
    if( spftrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spftrf is wrong\n",
                    spftrf_obj.info );
    }
    if( spftrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spftrf is wrong\n",
        spftrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( spftrf_obj.n*(spftrf_obj.n+1)/2,spftrf_obj.a,spftrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pftrf_scomplex_parameters  class definition */
class pftrf_scomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      char transr; // Must be 'N','T' (for real data) or 'C' (for complex data)
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_float *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pftrf_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                    char transr_i,lapack_int n_i);
      ~pftrf_scomplex_parameters (); 
};  /* end of pftrf_scomplex_parameters  class definition */


/* Constructor pftrf_scomplex_parameters definition */
pftrf_scomplex_parameters:: pftrf_scomplex_parameters ( int matrix_layout_i,
                                  char uplo_i,char transr_i,lapack_int n_i){
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    transr = transr_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*(n+1)/2)); 

    if( (a==NULL) || (aref==NULL) ){
       pftrf_free();
       printf(" pftrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,(n*(n+1)/2));

   } /* end of Constructor  */

pftrf_scomplex_parameters:: ~pftrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftrf_scomplex_parameters object: destructor invoked. \n");
#endif
   pftrf_free();
}

TEST(pftrf,cpftrf1) {

    /* LAPACKE CPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpftrf) ( int matrix_layout,
             char transr,char uplo,lapack_int n,lapack_complex_float *a  );

    Fptr_NL_LAPACKE_cpftrf CPFTRF;
    void *hModule,*dModule;
    float diff;
    pftrf_scomplex_parameters   cpftrf_obj(LAPACK_COL_MAJOR,'U','C',521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CPFTRF = (Fptr_NL_LAPACKE_cpftrf)dlsym(hModule,"LAPACKE_cpftrf");
    if (NULL == CPFTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cpftrf_obj.inforef = CPFTRF( cpftrf_obj.matrix_layout,cpftrf_obj.transr,
                            cpftrf_obj.uplo,cpftrf_obj.n,cpftrf_obj.aref);

        /* Compute libflame's Lapacke o/p  */
    cpftrf_obj.info     = LAPACKE_cpftrf( cpftrf_obj.matrix_layout,
                                     cpftrf_obj.transr,cpftrf_obj.uplo,
                                     cpftrf_obj.n,cpftrf_obj.a);
    if( cpftrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpftrf is wrong\n",
                    cpftrf_obj.info );
    }
    if( cpftrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpftrf is wrong\n",
        cpftrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( cpftrf_obj.n*(cpftrf_obj.n+1)/2,cpftrf_obj.a,cpftrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pftrf_dcomplex_parameters  class definition */
class pftrf_dcomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      char transr; // Must be 'N','T' (for real data) or 'C' (for complex data)
      lapack_int n; // No of rows,Columns

      /* Input/ Output parameters */
      lapack_complex_double *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pftrf_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                    char transr_i,lapack_int n_i);
      ~pftrf_dcomplex_parameters (); 
};  /* end of pftrf_dcomplex_parameters  class definition */


/* Constructor pftrf_dcomplex_parameters definition */
pftrf_dcomplex_parameters:: pftrf_dcomplex_parameters ( int matrix_layout_i,
                                  char uplo_i,char transr_i,lapack_int n_i){
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    transr = transr_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*(n+1)/2)); 

    if( (a==NULL) || (aref==NULL) ){
       pftrf_free();
       printf(" pftrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,(n*(n+1)/2));

   } /* end of Constructor  */

pftrf_dcomplex_parameters:: ~pftrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   pftrf_free();
}

TEST(pftrf,zpftrf1) {

    /* LAPACKE ZPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpftrf) ( int matrix_layout,
             char transr,char uplo,lapack_int n,lapack_complex_double *a  );

    Fptr_NL_LAPACKE_zpftrf ZPFTRF;
    void *hModule,*dModule;
    double diff;
    pftrf_dcomplex_parameters   zpftrf_obj(LAPACK_ROW_MAJOR,'U','C',121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZPFTRF = (Fptr_NL_LAPACKE_zpftrf)dlsym(hModule,"LAPACKE_zpftrf");
    if (NULL == ZPFTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zpftrf_obj.inforef = ZPFTRF( zpftrf_obj.matrix_layout,zpftrf_obj.transr,
                            zpftrf_obj.uplo,zpftrf_obj.n,zpftrf_obj.aref);

        /* Compute libflame's Lapacke o/p  */
    zpftrf_obj.info     = LAPACKE_zpftrf( zpftrf_obj.matrix_layout,
                                     zpftrf_obj.transr,zpftrf_obj.uplo,
                                     zpftrf_obj.n,zpftrf_obj.a);
    if( zpftrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpftrf is wrong\n",
                    zpftrf_obj.info );
    }
    if( zpftrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpftrf is wrong\n",
        zpftrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zpftrf_obj.n*(zpftrf_obj.n+1)/2,zpftrf_obj.a,zpftrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
