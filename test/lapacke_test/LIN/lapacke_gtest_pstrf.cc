#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define pstrf_free() \
       free (a   ); \
       free (aref); \
       free (piv); \
       free (pivref)

/* Begin pstrf_double_parameters  class definition */
class pstrf_double_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'
      double tol; // User defined tolerance

      /* Input/ Output parameters */
      double *a,*aref; //The array ab contains the matrix A

      /* Output parameters */
      lapack_int *piv,*pivref; // The pivot indices
      lapack_int rank; // Rank of a given by the number of steps
      /* Return Values */
      lapack_int info,inforef;

   public: 
      pstrf_double_parameters ( int matrix_layout_i,char uplo_i,
                    lapack_int n_i,lapack_int lda_i,double tol_i);
              
      ~pstrf_double_parameters (); 
};  /* end of pstrf_double_parameters  class definition */


/* Constructor pstrf_double_parameters definition */
pstrf_double_parameters:: pstrf_double_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i,double tol_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    tol = tol_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (lda*n));
    lapacke_gtest_alloc_int_buffer_pair ( &piv,&pivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (piv==NULL) || (pivref==NULL) ){
       pstrf_free();
       printf(" pstrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

pstrf_double_parameters:: ~pstrf_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pstrf_double_parameters object: destructor invoked. \n");
#endif
   pstrf_free();
}

TEST(pstrf,dpstrf1) {

    /* LAPACKE DPSTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpstrf) ( int matrix_layout,char uplo,
                                 lapack_int n,double *a,lapack_int lda,
                        lapack_int *piv,lapack_int* rank,double tol  );

    Fptr_NL_LAPACKE_dpstrf DPSTRF;
    void *hModule,*dModule;
    double diff;
    int piv_diff;
    pstrf_double_parameters   dpstrf_obj(LAPACK_ROW_MAJOR,'U',451,521,2.0);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DPSTRF = (Fptr_NL_LAPACKE_dpstrf)dlsym(hModule,"LAPACKE_dpstrf");
    if (NULL == DPSTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dpstrf_obj.inforef = DPSTRF( dpstrf_obj.matrix_layout,dpstrf_obj.uplo,
                             dpstrf_obj.n,dpstrf_obj.aref,dpstrf_obj.lda,
                      dpstrf_obj.pivref,&dpstrf_obj.rank,dpstrf_obj.tol);

    /* Compute libflame's Lapacke o/p  */
    dpstrf_obj.info = LAPACKE_dpstrf( dpstrf_obj.matrix_layout,dpstrf_obj.uplo,
                                     dpstrf_obj.n,dpstrf_obj.a,dpstrf_obj.lda,
                           dpstrf_obj.piv,&dpstrf_obj.rank,dpstrf_obj.tol);

    if( dpstrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpstrf is wrong\n",
                    dpstrf_obj.info );
    }
    if( dpstrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpstrf is wrong\n",
        dpstrf_obj.inforef );
    }
    piv_diff = computeDiff_i( dpstrf_obj.n,dpstrf_obj.piv,dpstrf_obj.pivref);
    if( piv_diff >0){
        printf("\n warning: pivot computation in dpstrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dpstrf_obj.n*dpstrf_obj.lda,dpstrf_obj.a,dpstrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pstrf_float_parameters  class definition */
class pstrf_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'
      float tol; // User defined tolerance

      /* Input/ Output parameters */
      float *a,*aref; //The array ab contains the matrix A

      /* Output parameters */
      lapack_int *piv,*pivref; // The pivot indices
      lapack_int rank; // Rank of a given by the number of steps

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pstrf_float_parameters ( int matrix_layout_i,char uplo_i,
                  lapack_int n_i,lapack_int lda_i,float tol_i);
      ~pstrf_float_parameters (); 
};  /* end of pstrf_float_parameters  class definition */


/* Constructor pstrf_float_parameters definition */
pstrf_float_parameters:: pstrf_float_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i,float tol_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    tol = tol_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_int_buffer_pair ( &piv,&pivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (piv==NULL) || (pivref==NULL) ){
       pstrf_free();
       printf(" pstrf_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

pstrf_float_parameters:: ~pstrf_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pstrf_float_parameters object: destructor invoked. \n");
#endif
   pstrf_free();
}

TEST(pstrf,spstrf1) {

    /* LAPACKE SPSTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spstrf) ( int matrix_layout ,char uplo,
                                 lapack_int n ,float *a ,lapack_int lda,
                           lapack_int *piv,lapack_int* rank,float tol);

    Fptr_NL_LAPACKE_spstrf SPSTRF;
    void *hModule,*dModule;
    float diff;
    int piv_diff;

    pstrf_float_parameters   spstrf_obj(LAPACK_COL_MAJOR,'U',1020,1041,1.0);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    SPSTRF = (Fptr_NL_LAPACKE_spstrf)dlsym(hModule,"LAPACKE_spstrf");
    if (NULL == SPSTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    spstrf_obj.inforef = SPSTRF( spstrf_obj.matrix_layout,spstrf_obj.uplo,
                            spstrf_obj.n,spstrf_obj.aref,spstrf_obj.lda,
                      spstrf_obj.pivref,&spstrf_obj.rank,spstrf_obj.tol);

        /* Compute libflame's Lapacke o/p  */
    spstrf_obj.info = LAPACKE_spstrf(spstrf_obj.matrix_layout,spstrf_obj.uplo,
                                    spstrf_obj.n,spstrf_obj.a,spstrf_obj.lda,
                             spstrf_obj.piv,&spstrf_obj.rank,spstrf_obj.tol);

    if( spstrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spstrf is wrong\n",
                    spstrf_obj.info );
    }
    if( spstrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spstrf is wrong\n",
        spstrf_obj.inforef );
    }

    piv_diff = computeDiff_i( spstrf_obj.n,spstrf_obj.piv,spstrf_obj.pivref);
    if( piv_diff >0){
        printf("\n warning: pivot computation in dpstrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( spstrf_obj.n*spstrf_obj.lda,spstrf_obj.a,
                                                   spstrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pstrf_scomplex_parameters  class definition */
class pstrf_scomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'
      float tol; // User defined tolerance

      /* Input/ Output parameters */
      lapack_complex_float *a,*aref; //The array ab contains the matrix A
      /* Output parameters */
      lapack_int *piv,*pivref; // The pivot indices
      lapack_int rank; // Rank of a given by the number of steps

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pstrf_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i,float tol_i);
      ~pstrf_scomplex_parameters (); 
};  /* end of pstrf_scomplex_parameters  class definition */


/* Constructor pstrf_scomplex_parameters definition */
pstrf_scomplex_parameters:: pstrf_scomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i,float tol_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    tol = tol_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_int_buffer_pair ( &piv,&pivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (piv==NULL) || (pivref==NULL) ){
       pstrf_free();
       printf(" pstrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

pstrf_scomplex_parameters:: ~pstrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pstrf_scomplex_parameters object: destructor invoked. \n");
#endif
   pstrf_free();
}

TEST(pstrf,cpstrf1) {

    /* LAPACKE CPSTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpstrf) ( int matrix_layout ,char uplo ,
                   lapack_int n ,lapack_complex_float *a ,lapack_int lda,
                        lapack_int *piv,lapack_int* rank,float tol  );

    Fptr_NL_LAPACKE_cpstrf CPSTRF;
    void *hModule,*dModule;
    float diff;
    int piv_diff;
    pstrf_scomplex_parameters   cpstrf_obj(LAPACK_COL_MAJOR,'U',510,521,1.0);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CPSTRF = (Fptr_NL_LAPACKE_cpstrf)dlsym(hModule,"LAPACKE_cpstrf");
    if (NULL == CPSTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cpstrf_obj.inforef = CPSTRF( cpstrf_obj.matrix_layout,cpstrf_obj.uplo,
                            cpstrf_obj.n,cpstrf_obj.aref,cpstrf_obj.lda,
                      cpstrf_obj.pivref,&cpstrf_obj.rank,cpstrf_obj.tol);

        /* Compute libflame's Lapacke o/p  */
    cpstrf_obj.info = LAPACKE_cpstrf(cpstrf_obj.matrix_layout,cpstrf_obj.uplo,
                                    cpstrf_obj.n,cpstrf_obj.a,cpstrf_obj.lda,
                             cpstrf_obj.piv,&cpstrf_obj.rank,cpstrf_obj.tol);

    if( cpstrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpstrf is wrong\n",
                    cpstrf_obj.info );
    }
    if( cpstrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpstrf is wrong\n",
        cpstrf_obj.inforef );
    }

    piv_diff = computeDiff_i( cpstrf_obj.n,cpstrf_obj.piv,cpstrf_obj.pivref);
    if( piv_diff >0){
        printf("\n warning: pivot computation in dpstrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( cpstrf_obj.n*cpstrf_obj.lda,cpstrf_obj.a,cpstrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin pstrf_dcomplex_parameters  class definition */
class pstrf_dcomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'
      double tol; // User defined tolerance

      /* Input/ Output parameters */
      lapack_complex_double *a,*aref; //The array ab contains the matrix A
      /* Output parameters */
      lapack_int *piv,*pivref; // The pivot indices
      lapack_int rank; // Rank of a given by the number of steps

      /* Return Values */
      lapack_int info,inforef;

   public: 
      pstrf_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i,double tol_i);
      ~pstrf_dcomplex_parameters (); 
};  /* end of pstrf_dcomplex_parameters  class definition */


/* Constructor pstrf_dcomplex_parameters definition */
pstrf_dcomplex_parameters:: pstrf_dcomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i,double tol_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    tol = tol_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_int_buffer_pair ( &piv,&pivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (piv==NULL) || (pivref==NULL) ){
       pstrf_free();
       printf(" pstrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

pstrf_dcomplex_parameters:: ~pstrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pstrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   pstrf_free();
}

TEST(pstrf,zpstrf1) {

    /* LAPACKE ZPSTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpstrf) ( int matrix_layout ,char uplo ,
                  lapack_int n ,lapack_complex_double *a ,lapack_int lda,
                        lapack_int *piv,lapack_int* rank,double tol  );

    Fptr_NL_LAPACKE_zpstrf ZPSTRF;
    void *hModule,*dModule;
    double diff;
    int piv_diff;
    pstrf_dcomplex_parameters   zpstrf_obj(LAPACK_ROW_MAJOR,'U',100,121,0.5);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZPSTRF = (Fptr_NL_LAPACKE_zpstrf)dlsym(hModule,"LAPACKE_zpstrf");
    if (NULL == ZPSTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zpstrf_obj.inforef = ZPSTRF( zpstrf_obj.matrix_layout,zpstrf_obj.uplo,
                            zpstrf_obj.n,zpstrf_obj.aref,zpstrf_obj.lda,
                      zpstrf_obj.pivref,&zpstrf_obj.rank,zpstrf_obj.tol);

        /* Compute libflame's Lapacke o/p  */
    zpstrf_obj.info  = LAPACKE_zpstrf( zpstrf_obj.matrix_layout,zpstrf_obj.uplo,
                                      zpstrf_obj.n,zpstrf_obj.a,zpstrf_obj.lda,
                               zpstrf_obj.piv,&zpstrf_obj.rank,zpstrf_obj.tol);

    if( zpstrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpstrf is wrong\n",
                    zpstrf_obj.info );
    }
    if( zpstrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpstrf is wrong\n",
        zpstrf_obj.inforef );
    }

    piv_diff = computeDiff_i( zpstrf_obj.n,zpstrf_obj.piv,zpstrf_obj.pivref);
    if( piv_diff >0){
        printf("\n warning: pivot computation in dpstrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zpstrf_obj.n*zpstrf_obj.lda,zpstrf_obj.a,zpstrf_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
