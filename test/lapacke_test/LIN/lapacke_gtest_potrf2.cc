#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define potrf2_free() \
       free (a   ); \
       free (aref)

/* Begin potrf2_double_parameters  class definition */
class potrf2_double_parameters{

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
      potrf2_double_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~potrf2_double_parameters (); 
};  /* end of potrf2_double_parameters  class definition */


/* Constructor potrf2_double_parameters definition */
potrf2_double_parameters:: potrf2_double_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (lda*n));

    if( (a==NULL) || (aref==NULL) ){
       potrf2_free();
       printf(" potrf2_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

potrf2_double_parameters:: ~potrf2_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrf2_double_parameters object: destructor invoked. \n");
#endif
   potrf2_free();
}

TEST(potrf2,dpotrf21) {

    /* LAPACKE DPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpotrf2) ( int matrix_layout ,char uplo ,
                             lapack_int n ,double *a ,lapack_int lda  );

    Fptr_NL_LAPACKE_dpotrf2 DPOTRF;
    void *hModule,*dModule;
    double diff;
    potrf2_double_parameters   dpotrf2_obj(LAPACK_ROW_MAJOR,'U',451,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DPOTRF = (Fptr_NL_LAPACKE_dpotrf2)dlsym(hModule,"LAPACKE_dpotrf2");
    if (NULL == DPOTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dpotrf2_obj.inforef = DPOTRF( dpotrf2_obj.matrix_layout,dpotrf2_obj.uplo,
                            dpotrf2_obj.n,dpotrf2_obj.aref,dpotrf2_obj.lda);

        /* Compute libflame's Lapacke o/p  */
    dpotrf2_obj.info     = LAPACKE_dpotrf2( dpotrf2_obj.matrix_layout,
                                 dpotrf2_obj.uplo,dpotrf2_obj.n,
                         dpotrf2_obj.a,dpotrf2_obj.lda);
    if( dpotrf2_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpotrf2 is wrong\n",
                    dpotrf2_obj.info );
    }
    if( dpotrf2_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpotrf2 is wrong\n",
        dpotrf2_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dpotrf2_obj.n*dpotrf2_obj.lda,dpotrf2_obj.a,dpotrf2_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin potrf2_float_parameters  class definition */
class potrf2_float_parameters{

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
      potrf2_float_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~potrf2_float_parameters (); 
};  /* end of potrf2_float_parameters  class definition */


/* Constructor potrf2_float_parameters definition */
potrf2_float_parameters:: potrf2_float_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (lda*n)); 

    if( (a==NULL) || (aref==NULL) ){
       potrf2_free();
       printf(" potrf2_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

potrf2_float_parameters:: ~potrf2_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrf2_float_parameters object: destructor invoked. \n");
#endif
   potrf2_free();
}

TEST(potrf2,spotrf21) {

    /* LAPACKE SPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spotrf2) ( int matrix_layout ,char uplo ,
                             lapack_int n ,float *a ,lapack_int lda  );

    Fptr_NL_LAPACKE_spotrf2 SPOTRF;
    void *hModule,*dModule;
    float diff;
    potrf2_float_parameters   spotrf2_obj(LAPACK_COL_MAJOR,'U',1020,1041);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SPOTRF = (Fptr_NL_LAPACKE_spotrf2)dlsym(hModule,"LAPACKE_spotrf2");
    if (NULL == SPOTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    spotrf2_obj.inforef = SPOTRF( spotrf2_obj.matrix_layout,spotrf2_obj.uplo,
                            spotrf2_obj.n,spotrf2_obj.aref,spotrf2_obj.lda);

        /* Compute libflame's Lapacke o/p  */
    spotrf2_obj.info     = LAPACKE_spotrf2( spotrf2_obj.matrix_layout,
                                 spotrf2_obj.uplo,spotrf2_obj.n,
                         spotrf2_obj.a,spotrf2_obj.lda);
    if( spotrf2_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spotrf2 is wrong\n",
                    spotrf2_obj.info );
    }
    if( spotrf2_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spotrf2 is wrong\n",
        spotrf2_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( spotrf2_obj.n*spotrf2_obj.lda,spotrf2_obj.a,spotrf2_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin potrf2_scomplex_parameters  class definition */
class potrf2_scomplex_parameters{

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
      potrf2_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~potrf2_scomplex_parameters (); 
};  /* end of potrf2_scomplex_parameters  class definition */


/* Constructor potrf2_scomplex_parameters definition */
potrf2_scomplex_parameters:: potrf2_scomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (lda*n)); 

    if( (a==NULL) || (aref==NULL) ){
       potrf2_free();
       printf(" potrf2_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

potrf2_scomplex_parameters:: ~potrf2_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrf2_scomplex_parameters object: destructor invoked. \n");
#endif
   potrf2_free();
}

TEST(potrf2,cpotrf21) {

    /* LAPACKE CPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpotrf2) ( int matrix_layout ,char uplo ,
                             lapack_int n ,lapack_complex_float *a ,lapack_int lda  );

    Fptr_NL_LAPACKE_cpotrf2 CPOTRF;
    void *hModule,*dModule;
    float diff;
    potrf2_scomplex_parameters   cpotrf2_obj(LAPACK_COL_MAJOR,'U',510,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CPOTRF = (Fptr_NL_LAPACKE_cpotrf2)dlsym(hModule,"LAPACKE_cpotrf2");
    if (NULL == CPOTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cpotrf2_obj.inforef = CPOTRF( cpotrf2_obj.matrix_layout,cpotrf2_obj.uplo,
                            cpotrf2_obj.n,cpotrf2_obj.aref,cpotrf2_obj.lda);

        /* Compute libflame's Lapacke o/p  */
    cpotrf2_obj.info     = LAPACKE_cpotrf2( cpotrf2_obj.matrix_layout,
                                     cpotrf2_obj.uplo,cpotrf2_obj.n,
                                     cpotrf2_obj.a,cpotrf2_obj.lda);
    if( cpotrf2_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpotrf2 is wrong\n",
                    cpotrf2_obj.info );
    }
    if( cpotrf2_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpotrf2 is wrong\n",
        cpotrf2_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( cpotrf2_obj.n*cpotrf2_obj.lda,cpotrf2_obj.a,cpotrf2_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin potrf2_dcomplex_parameters  class definition */
class potrf2_dcomplex_parameters{

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
      potrf2_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~potrf2_dcomplex_parameters (); 
};  /* end of potrf2_dcomplex_parameters  class definition */


/* Constructor potrf2_dcomplex_parameters definition */
potrf2_dcomplex_parameters:: potrf2_dcomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (lda*n)); 

    if( (a==NULL) || (aref==NULL) ){
       potrf2_free();
       printf(" potrf2_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

potrf2_dcomplex_parameters:: ~potrf2_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrf2_dcomplex_parameters object: destructor invoked. \n");
#endif
   potrf2_free();
}

TEST(potrf2,zpotrf21) {

    /* LAPACKE ZPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpotrf2) ( int matrix_layout ,char uplo ,
                             lapack_int n ,lapack_complex_double *a ,lapack_int lda  );

    Fptr_NL_LAPACKE_zpotrf2 ZPOTRF;
    void *hModule,*dModule;
    double diff;
    potrf2_dcomplex_parameters   zpotrf2_obj(LAPACK_ROW_MAJOR,'U',100,121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZPOTRF = (Fptr_NL_LAPACKE_zpotrf2)dlsym(hModule,"LAPACKE_zpotrf2");
    if (NULL == ZPOTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zpotrf2_obj.inforef = ZPOTRF( zpotrf2_obj.matrix_layout,zpotrf2_obj.uplo,
                            zpotrf2_obj.n,zpotrf2_obj.aref,zpotrf2_obj.lda);

        /* Compute libflame's Lapacke o/p  */
    zpotrf2_obj.info     = LAPACKE_zpotrf2( zpotrf2_obj.matrix_layout,
                                     zpotrf2_obj.uplo,zpotrf2_obj.n,
                                     zpotrf2_obj.a,zpotrf2_obj.lda);
    if( zpotrf2_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpotrf2 is wrong\n",
                    zpotrf2_obj.info );
    }
    if( zpotrf2_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpotrf2 is wrong\n",
        zpotrf2_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zpotrf2_obj.n*zpotrf2_obj.lda,zpotrf2_obj.a,zpotrf2_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
