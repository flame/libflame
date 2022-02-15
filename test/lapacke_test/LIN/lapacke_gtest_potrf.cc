#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define potrf_free() \
       free (a   ); \
       free (aref)

/* Begin potrf_double_parameters  class definition */
class potrf_double_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      double *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      potrf_double_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~potrf_double_parameters (); 
};  /* end of potrf_double_parameters  class definition */


/* Constructor potrf_double_parameters definition */
potrf_double_parameters:: potrf_double_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (lda*n));

    if( (a==NULL) || (aref==NULL) ){
       potrf_free();
       printf(" potrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

potrf_double_parameters:: ~potrf_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrf_double_parameters object: destructor invoked. \n");
#endif
   potrf_free();
}

TEST(potrf,dpotrf1) {

    /* LAPACKE DPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpotrf) ( int matrix_layout ,char uplo ,
                             lapack_int n ,double *a ,lapack_int lda  );
                    
    Fptr_NL_LAPACKE_dpotrf DPOTRF;
    void *hModule,*dModule;
    double diff;
    potrf_double_parameters   dpotrf_obj(LAPACK_ROW_MAJOR,'U',451,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DPOTRF = (Fptr_NL_LAPACKE_dpotrf)dlsym(hModule,"LAPACKE_dpotrf");
    if (NULL == DPOTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }
    
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    dpotrf_obj.inforef = DPOTRF( dpotrf_obj.matrix_layout,dpotrf_obj.uplo,
                            dpotrf_obj.n,dpotrf_obj.aref,dpotrf_obj.lda);

        /* Compute libflame's Lapacke o/p  */
    dpotrf_obj.info     = LAPACKE_dpotrf( dpotrf_obj.matrix_layout,
                                 dpotrf_obj.uplo,dpotrf_obj.n,
                         dpotrf_obj.a,dpotrf_obj.lda);
    if( dpotrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpotrf is wrong\n",
                    dpotrf_obj.info );
    }
    if( dpotrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpotrf is wrong\n",
        dpotrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dpotrf_obj.n*dpotrf_obj.lda,dpotrf_obj.a,dpotrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin potrf_float_parameters  class definition */
class potrf_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      float *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      potrf_float_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~potrf_float_parameters (); 
};  /* end of potrf_float_parameters  class definition */


/* Constructor potrf_float_parameters definition */
potrf_float_parameters:: potrf_float_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (lda*n));

    if( (a==NULL) || (aref==NULL) ){
       potrf_free();
       printf(" potrf_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

potrf_float_parameters:: ~potrf_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrf_float_parameters object: destructor invoked. \n");
#endif
   potrf_free();
}

TEST(potrf,spotrf1) {

    /* LAPACKE SPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spotrf) ( int matrix_layout ,char uplo ,
                             lapack_int n ,float *a ,lapack_int lda  );
                    
    Fptr_NL_LAPACKE_spotrf SPOTRF;
    void *hModule,*dModule;
    float diff;
    potrf_float_parameters   spotrf_obj(LAPACK_COL_MAJOR,'U',1020,1041);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SPOTRF = (Fptr_NL_LAPACKE_spotrf)dlsym(hModule,"LAPACKE_spotrf");
    if (NULL == SPOTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }
    
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    spotrf_obj.inforef = SPOTRF( spotrf_obj.matrix_layout,spotrf_obj.uplo,
                            spotrf_obj.n,spotrf_obj.aref,spotrf_obj.lda);

        /* Compute libflame's Lapacke o/p  */
    spotrf_obj.info     = LAPACKE_spotrf( spotrf_obj.matrix_layout,
                                 spotrf_obj.uplo,spotrf_obj.n,
                         spotrf_obj.a,spotrf_obj.lda);
    if( spotrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spotrf is wrong\n",
                    spotrf_obj.info );
    }
    if( spotrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spotrf is wrong\n",
        spotrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( spotrf_obj.n*spotrf_obj.lda,spotrf_obj.a,spotrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin potrf_scomplex_parameters  class definition */
class potrf_scomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_float *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      potrf_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~potrf_scomplex_parameters (); 
};  /* end of potrf_scomplex_parameters  class definition */


/* Constructor potrf_scomplex_parameters definition */
potrf_scomplex_parameters:: potrf_scomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (lda*n)); 

    if( (a==NULL) || (aref==NULL) ){
       potrf_free();
       printf(" potrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

potrf_scomplex_parameters:: ~potrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrf_scomplex_parameters object: destructor invoked. \n");
#endif
   potrf_free();
}

TEST(potrf,cpotrf1) {

    /* LAPACKE CPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpotrf) ( int matrix_layout ,char uplo ,
                             lapack_int n ,lapack_complex_float *a ,lapack_int lda  );

    Fptr_NL_LAPACKE_cpotrf CPOTRF;
    void *hModule,*dModule;
    float diff;
    potrf_scomplex_parameters   cpotrf_obj(LAPACK_COL_MAJOR,'U',510,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CPOTRF = (Fptr_NL_LAPACKE_cpotrf)dlsym(hModule,"LAPACKE_cpotrf");
    if (NULL == CPOTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cpotrf_obj.inforef = CPOTRF( cpotrf_obj.matrix_layout,cpotrf_obj.uplo,
                            cpotrf_obj.n,cpotrf_obj.aref,cpotrf_obj.lda);

        /* Compute libflame's Lapacke o/p  */
    cpotrf_obj.info     = LAPACKE_cpotrf( cpotrf_obj.matrix_layout,
                                     cpotrf_obj.uplo,cpotrf_obj.n,
                                     cpotrf_obj.a,cpotrf_obj.lda);
    if( cpotrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpotrf is wrong\n",
                    cpotrf_obj.info );
    }
    if( cpotrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpotrf is wrong\n",
        cpotrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( cpotrf_obj.n*cpotrf_obj.lda,cpotrf_obj.a,cpotrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin potrf_dcomplex_parameters  class definition */
class potrf_dcomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_double *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      potrf_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~potrf_dcomplex_parameters (); 
};  /* end of potrf_dcomplex_parameters  class definition */


/* Constructor potrf_dcomplex_parameters definition */
potrf_dcomplex_parameters:: potrf_dcomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (lda*n)); 

    if( (a==NULL) || (aref==NULL) ){
       potrf_free();
       printf(" potrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

potrf_dcomplex_parameters:: ~potrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   potrf_free();
}

TEST(potrf,zpotrf1) {

    /* LAPACKE ZPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpotrf) ( int matrix_layout ,char uplo ,
                             lapack_int n ,lapack_complex_double *a ,lapack_int lda  );

    Fptr_NL_LAPACKE_zpotrf ZPOTRF;
    void *hModule,*dModule;
    double diff;
    potrf_dcomplex_parameters   zpotrf_obj(LAPACK_ROW_MAJOR,'U',100,121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZPOTRF = (Fptr_NL_LAPACKE_zpotrf)dlsym(hModule,"LAPACKE_zpotrf");
    if (NULL == ZPOTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zpotrf_obj.inforef = ZPOTRF( zpotrf_obj.matrix_layout,zpotrf_obj.uplo,
                            zpotrf_obj.n,zpotrf_obj.aref,zpotrf_obj.lda);

        /* Compute libflame's Lapacke o/p  */
    zpotrf_obj.info     = LAPACKE_zpotrf( zpotrf_obj.matrix_layout,
                                     zpotrf_obj.uplo,zpotrf_obj.n,
                                     zpotrf_obj.a,zpotrf_obj.lda);
    if( zpotrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpotrf is wrong\n",
                    zpotrf_obj.info );
    }
    if( zpotrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpotrf is wrong\n",
        zpotrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zpotrf_obj.n*zpotrf_obj.lda,zpotrf_obj.a,zpotrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
