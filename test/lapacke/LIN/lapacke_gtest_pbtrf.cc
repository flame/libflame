#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define pbtrf_free() \
       free (ab   ); \
       free (abref)

/* Begin pbtrf_double_parameters  class definition */
class pbtrf_double_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int kd; // No. of superdiagonals or subdiagonals in the matrix A
      /* Input/ Output parameters */
      double *ab,*abref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pbtrf_double_parameters(int matrix_layout_i,char uplo_i,lapack_int n_i,
                                             lapack_int kd_i,lapack_int ldab_i);
      ~pbtrf_double_parameters (); 
};  /* end of pbtrf_double_parameters  class definition */


/* Constructor pbtrf_double_parameters definition */
pbtrf_double_parameters:: pbtrf_double_parameters ( int matrix_layout_i,
                           char uplo_i,lapack_int n_i,lapack_int kd_i,
                                                     lapack_int ldab_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldab = ldab_i;
    kd = kd_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &ab, &abref, (ldab*n));

    if( (ab==NULL) || (abref==NULL) ){
       pbtrf_free();
       printf(" pbtrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ab, abref, (ldab*n));

   } /* end of Constructor  */

pbtrf_double_parameters:: ~pbtrf_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbtrf_double_parameters object: destructor invoked. \n");
#endif
   pbtrf_free();
}

TEST(pbtrf,dpbtrf1) {

    /* LAPACKE DPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpbtrf) ( int matrix_layout,char uplo,
              lapack_int n,lapack_int kd,double *ab,lapack_int ldab  );
                    
    Fptr_NL_LAPACKE_dpbtrf DPBTRF;
    void *hModule,*dModule;
    double diff;
    //pbtrf_double_parameters   dpbtrf_obj(LAPACK_ROW_MAJOR,'U',451,21,30);
    pbtrf_double_parameters   dpbtrf_obj(LAPACK_ROW_MAJOR,'U',1020,104,105  );

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DPBTRF = (Fptr_NL_LAPACKE_dpbtrf)dlsym(hModule,"LAPACKE_dpbtrf");
    if (NULL == DPBTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }
    
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    dpbtrf_obj.inforef = DPBTRF( dpbtrf_obj.matrix_layout,dpbtrf_obj.uplo,
                              dpbtrf_obj.n,dpbtrf_obj.kd,dpbtrf_obj.abref,
                                                           dpbtrf_obj.ldab);

        /* Compute libflame's Lapacke o/p  */
    dpbtrf_obj.info     = LAPACKE_dpbtrf( dpbtrf_obj.matrix_layout,
                      dpbtrf_obj.uplo,dpbtrf_obj.n,dpbtrf_obj.kd,
                                   dpbtrf_obj.ab,dpbtrf_obj.ldab);

    if( dpbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpbtrf is wrong\n",
                    dpbtrf_obj.info );
    }
    if( dpbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpbtrf is wrong\n",
        dpbtrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dpbtrf_obj.n*dpbtrf_obj.ldab,dpbtrf_obj.ab,dpbtrf_obj.abref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pbtrf_float_parameters  class definition */
class pbtrf_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int kd; // No. of superdiagonals or subdiagonals in the matrix A

      /* Input/ Output parameters */
      float *ab,*abref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pbtrf_float_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int kd_i, lapack_int ldab_i);
      ~pbtrf_float_parameters (); 
};  /* end of pbtrf_float_parameters  class definition */


/* Constructor pbtrf_float_parameters definition */
pbtrf_float_parameters:: pbtrf_float_parameters ( int matrix_layout_i,
                         char uplo_i,lapack_int n_i,lapack_int kd_i,
                                       lapack_int ldab_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldab = ldab_i;
    kd = kd_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &ab, &abref, (ldab*n));

    if( (ab==NULL) || (abref==NULL) ){
       pbtrf_free();
       printf(" pbtrf_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ab, abref, (ldab*n));

   } /* end of Constructor  */

pbtrf_float_parameters:: ~pbtrf_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbtrf_float_parameters object: destructor invoked. \n");
#endif
   pbtrf_free();
}

TEST(pbtrf,spbtrf1) {

    /* LAPACKE SPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spbtrf) ( int matrix_layout,char uplo,
                             lapack_int n,lapack_int kd,float *ab,
                             lapack_int ldab  );
                    
    Fptr_NL_LAPACKE_spbtrf SPBTRF;
    void *hModule,*dModule;
    float diff;
    pbtrf_float_parameters   spbtrf_obj(LAPACK_COL_MAJOR,'L',1020,104,134);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SPBTRF = (Fptr_NL_LAPACKE_spbtrf)dlsym(hModule,"LAPACKE_spbtrf");
    if (NULL == SPBTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }
    
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    spbtrf_obj.inforef = SPBTRF( spbtrf_obj.matrix_layout,spbtrf_obj.uplo,
                             spbtrf_obj.n,spbtrf_obj.kd,spbtrf_obj.abref,
                                                           spbtrf_obj.ldab);
        /* Compute libflame's Lapacke o/p  */
    spbtrf_obj.info     = LAPACKE_spbtrf( spbtrf_obj.matrix_layout,
                             spbtrf_obj.uplo,spbtrf_obj.n,spbtrf_obj.kd,
                                          spbtrf_obj.ab,spbtrf_obj.ldab);
    if( spbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spbtrf is wrong\n",
                    spbtrf_obj.info );
    }
    if( spbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spbtrf is wrong\n",
        spbtrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( spbtrf_obj.n*spbtrf_obj.ldab,spbtrf_obj.ab,spbtrf_obj.abref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pbtrf_scomplex_parameters  class definition */
class pbtrf_scomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int kd; // No. of superdiagonals or subdiagonals in the matrix A

      /* Input/ Output parameters */
      lapack_complex_float *ab,*abref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pbtrf_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int kd_i, lapack_int ldab_i);
      ~pbtrf_scomplex_parameters (); 
};  /* end of pbtrf_scomplex_parameters  class definition */


/* Constructor pbtrf_scomplex_parameters definition */
pbtrf_scomplex_parameters:: pbtrf_scomplex_parameters ( int matrix_layout_i,
                                 char uplo_i,lapack_int n_i,lapack_int kd_i,
                                                         lapack_int ldab_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldab = ldab_i;
    kd = kd_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ab, &abref, (ldab*n)); 

    if( (ab==NULL) || (abref==NULL) ){
       pbtrf_free();
       printf(" pbtrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ab, abref, (ldab*n));

   } /* end of Constructor  */

pbtrf_scomplex_parameters:: ~pbtrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbtrf_scomplex_parameters object: destructor invoked. \n");
#endif
   pbtrf_free();
}

TEST(pbtrf,cpbtrf1) {

    /* LAPACKE CPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbtrf) ( int matrix_layout,char uplo,
                                             lapack_int n,lapack_int kd,
                              lapack_complex_float *ab,lapack_int ldab );

    Fptr_NL_LAPACKE_cpbtrf CPBTRF;
    void *hModule,*dModule;
    float diff;
    pbtrf_scomplex_parameters   cpbtrf_obj(LAPACK_COL_MAJOR,'U',510,221,300);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CPBTRF = (Fptr_NL_LAPACKE_cpbtrf)dlsym(hModule,"LAPACKE_cpbtrf");
    if (NULL == CPBTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cpbtrf_obj.inforef = CPBTRF( cpbtrf_obj.matrix_layout,cpbtrf_obj.uplo,
                             cpbtrf_obj.n,cpbtrf_obj.kd,cpbtrf_obj.abref,
                                                          cpbtrf_obj.ldab);

        /* Compute libflame's Lapacke o/p  */
    cpbtrf_obj.info = LAPACKE_cpbtrf(cpbtrf_obj.matrix_layout,cpbtrf_obj.uplo,
                                                   cpbtrf_obj.n,cpbtrf_obj.kd,
                                               cpbtrf_obj.ab,cpbtrf_obj.ldab);
    if( cpbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpbtrf is wrong\n",
                    cpbtrf_obj.info );
    }
    if( cpbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpbtrf is wrong\n",
        cpbtrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( cpbtrf_obj.n*cpbtrf_obj.ldab,cpbtrf_obj.ab,cpbtrf_obj.abref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pbtrf_dcomplex_parameters  class definition */
class pbtrf_dcomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int kd; // No. of superdiagonals or subdiagonals in the matrix A

      /* Input/ Output parameters */
      lapack_complex_double *ab,*abref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pbtrf_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                           lapack_int n_i,lapack_int kd_i,lapack_int ldab_i);
      ~pbtrf_dcomplex_parameters (); 
};  /* end of pbtrf_dcomplex_parameters  class definition */


/* Constructor pbtrf_dcomplex_parameters definition */
pbtrf_dcomplex_parameters:: pbtrf_dcomplex_parameters ( int matrix_layout_i,
                               char uplo_i,lapack_int n_i,lapack_int kd_i,
                                                        lapack_int ldab_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldab = ldab_i;
    kd = kd_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ab, &abref, (ldab*n)); 

    if( (ab==NULL) || (abref==NULL) ){
       pbtrf_free();
       printf(" pbtrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ab, abref, (ldab*n));

   } /* end of Constructor  */

pbtrf_dcomplex_parameters:: ~pbtrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbtrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   pbtrf_free();
}

TEST(pbtrf,zpbtrf1) {

    /* LAPACKE ZPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbtrf) ( int matrix_layout,char uplo,
                                             lapack_int n,lapack_int kd,
                             lapack_complex_double *ab,lapack_int ldab );

    Fptr_NL_LAPACKE_zpbtrf ZPBTRF;
    void *hModule,*dModule;
    double diff;
    pbtrf_dcomplex_parameters   zpbtrf_obj(LAPACK_ROW_MAJOR,'L',510,221,300);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZPBTRF = (Fptr_NL_LAPACKE_zpbtrf)dlsym(hModule,"LAPACKE_zpbtrf");
    if (NULL == ZPBTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zpbtrf_obj.inforef = ZPBTRF( zpbtrf_obj.matrix_layout,zpbtrf_obj.uplo,
                             zpbtrf_obj.n,zpbtrf_obj.kd,zpbtrf_obj.abref,
                                                          zpbtrf_obj.ldab);

        /* Compute libflame's Lapacke o/p  */
    zpbtrf_obj.info     = LAPACKE_zpbtrf( zpbtrf_obj.matrix_layout,
                      zpbtrf_obj.uplo,zpbtrf_obj.n,zpbtrf_obj.kd,
                                     zpbtrf_obj.ab,zpbtrf_obj.ldab);
    if( zpbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpbtrf is wrong\n",
                    zpbtrf_obj.info );
    }
    if( zpbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpbtrf is wrong\n",
        zpbtrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zpbtrf_obj.n*zpbtrf_obj.ldab,zpbtrf_obj.ab,zpbtrf_obj.abref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
