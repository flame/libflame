#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define hetrf_rk_free() \
       free (a   ); \
       free (aref); \
       free (e   ); \
       free (eref); \
       free (ipiv); \
       free (ipivref)


/* Begin hetrf_rk_scomplex_parameters  class definition */
class hetrf_rk_scomplex_parameters{

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
      hetrf_rk_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                 lapack_int n_i,lapack_int lda_i);
      ~hetrf_rk_scomplex_parameters (); 
};  /* end of hetrf_rk_scomplex_parameters  class definition */


/* Constructor hetrf_rk_scomplex_parameters definition */
hetrf_rk_scomplex_parameters:: hetrf_rk_scomplex_parameters ( int matrix_layout_i,
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
       hetrf_rk_free();
       printf(" hetrf_rk_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

hetrf_rk_scomplex_parameters:: ~hetrf_rk_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrf_rk_scomplex_parameters object: destructor invoked. \n");
#endif
   hetrf_rk_free();
}

TEST(hetrf_rk,chetrf_rk1) {

    /* LAPACKE CHETRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf_rk) ( int matrix_layout,char uplo,
                   lapack_int n,lapack_complex_float *a,lapack_int lda,
                              lapack_complex_float *e,lapack_int *ipiv);

    Fptr_NL_LAPACKE_chetrf_rk CHETRF_RK;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    hetrf_rk_scomplex_parameters   chetrf_rk_obj(LAPACK_COL_MAJOR,'U',510,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CHETRF_RK = (Fptr_NL_LAPACKE_chetrf_rk)dlsym(hModule,"LAPACKE_chetrf_rk");
    if (NULL == CHETRF_RK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chetrf_rk_obj.inforef = CHETRF_RK( chetrf_rk_obj.matrix_layout,chetrf_rk_obj.uplo,
                            chetrf_rk_obj.n,chetrf_rk_obj.aref,chetrf_rk_obj.lda,
                                        chetrf_rk_obj.eref,chetrf_rk_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    chetrf_rk_obj.info = LAPACKE_chetrf_rk(chetrf_rk_obj.matrix_layout,chetrf_rk_obj.uplo,
                                    chetrf_rk_obj.n,chetrf_rk_obj.a,chetrf_rk_obj.lda,
                                                  chetrf_rk_obj.e,chetrf_rk_obj.ipiv );

    if( chetrf_rk_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrf_rk is wrong\n",
                    chetrf_rk_obj.info );
    }
    if( chetrf_rk_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrf_rk is wrong\n",
        chetrf_rk_obj.inforef );
    }

    ipiv_diff = computeDiff_i( chetrf_rk_obj.n,chetrf_rk_obj.ipiv,chetrf_rk_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf_rk1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( chetrf_rk_obj.n*chetrf_rk_obj.lda,chetrf_rk_obj.a,
                                                         chetrf_rk_obj.aref );
    diff +=  computeDiff_c( chetrf_rk_obj.n,chetrf_rk_obj.e,chetrf_rk_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST(hetrf_rk,chetrf_rk2) {

    /* LAPACKE CHETRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf_rk) ( int matrix_layout,char uplo,
                   lapack_int n,lapack_complex_float *a,lapack_int lda,
                              lapack_complex_float *e,lapack_int *ipiv);

    Fptr_NL_LAPACKE_chetrf_rk CHETRF_RK;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    hetrf_rk_scomplex_parameters   chetrf_rk_obj(LAPACK_COL_MAJOR,'L',310,421);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CHETRF_RK = (Fptr_NL_LAPACKE_chetrf_rk)dlsym(hModule,"LAPACKE_chetrf_rk");
    if (NULL == CHETRF_RK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chetrf_rk_obj.inforef = CHETRF_RK( chetrf_rk_obj.matrix_layout,chetrf_rk_obj.uplo,
                            chetrf_rk_obj.n,chetrf_rk_obj.aref,chetrf_rk_obj.lda,
                                        chetrf_rk_obj.eref,chetrf_rk_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    chetrf_rk_obj.info = LAPACKE_chetrf_rk(chetrf_rk_obj.matrix_layout,chetrf_rk_obj.uplo,
                                    chetrf_rk_obj.n,chetrf_rk_obj.a,chetrf_rk_obj.lda,
                                                  chetrf_rk_obj.e,chetrf_rk_obj.ipiv );

    if( chetrf_rk_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrf_rk is wrong\n",
                    chetrf_rk_obj.info );
    }
    if( chetrf_rk_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrf_rk is wrong\n",
        chetrf_rk_obj.inforef );
    }

    ipiv_diff = computeDiff_i( chetrf_rk_obj.n,chetrf_rk_obj.ipiv,chetrf_rk_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf_rk1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( chetrf_rk_obj.n*chetrf_rk_obj.lda,chetrf_rk_obj.a,
                                                         chetrf_rk_obj.aref );
    diff +=  computeDiff_c( chetrf_rk_obj.n,chetrf_rk_obj.e,chetrf_rk_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin hetrf_rk_dcomplex_parameters  class definition */
class hetrf_rk_dcomplex_parameters{

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
      hetrf_rk_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~hetrf_rk_dcomplex_parameters (); 
};  /* end of hetrf_rk_dcomplex_parameters  class definition */


/* Constructor hetrf_rk_dcomplex_parameters definition */
hetrf_rk_dcomplex_parameters:: hetrf_rk_dcomplex_parameters ( int matrix_layout_i,
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
       hetrf_rk_free();
       printf(" hetrf_rk_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

hetrf_rk_dcomplex_parameters:: ~hetrf_rk_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrf_rk_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetrf_rk_free();
}

TEST(hetrf_rk,zhetrf_rk1) {

    /* LAPACKE ZHETRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf_rk) ( int matrix_layout,char uplo,
                  lapack_int n,lapack_complex_double *a,lapack_int lda,
                             lapack_complex_double *e,lapack_int *ipiv);

    Fptr_NL_LAPACKE_zhetrf_rk ZHETRF_RK;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    hetrf_rk_dcomplex_parameters   zhetrf_rk_obj(LAPACK_ROW_MAJOR,'U',100,121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZHETRF_RK = (Fptr_NL_LAPACKE_zhetrf_rk)dlsym(hModule,"LAPACKE_zhetrf_rk");
    if (NULL == ZHETRF_RK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhetrf_rk_obj.inforef = ZHETRF_RK( zhetrf_rk_obj.matrix_layout,zhetrf_rk_obj.uplo,
                            zhetrf_rk_obj.n,zhetrf_rk_obj.aref,zhetrf_rk_obj.lda,
                                        zhetrf_rk_obj.eref,zhetrf_rk_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zhetrf_rk_obj.info = LAPACKE_zhetrf_rk(zhetrf_rk_obj.matrix_layout,zhetrf_rk_obj.uplo,
                                   zhetrf_rk_obj.n,zhetrf_rk_obj.a,zhetrf_rk_obj.lda,
                                                  zhetrf_rk_obj.e,zhetrf_rk_obj.ipiv);

    if( zhetrf_rk_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrf_rk is wrong\n",
                    zhetrf_rk_obj.info );
    }
    if( zhetrf_rk_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrf_rk is wrong\n",
        zhetrf_rk_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zhetrf_rk_obj.n,zhetrf_rk_obj.ipiv,zhetrf_rk_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf_rk1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff  =  computeDiff_z( zhetrf_rk_obj.n*zhetrf_rk_obj.lda,zhetrf_rk_obj.a,
                                                          zhetrf_rk_obj.aref );
    diff +=  computeDiff_z( zhetrf_rk_obj.n,zhetrf_rk_obj.e,zhetrf_rk_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST(hetrf_rk,zhetrf_rk2) {

    /* LAPACKE ZHETRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf_rk) ( int matrix_layout,char uplo,
                  lapack_int n,lapack_complex_double *a,lapack_int lda,
                             lapack_complex_double *e,lapack_int *ipiv);

    Fptr_NL_LAPACKE_zhetrf_rk ZHETRF_RK;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    hetrf_rk_dcomplex_parameters   zhetrf_rk_obj(LAPACK_ROW_MAJOR,'L',1000,1021);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZHETRF_RK = (Fptr_NL_LAPACKE_zhetrf_rk)dlsym(hModule,"LAPACKE_zhetrf_rk");
    if (NULL == ZHETRF_RK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhetrf_rk_obj.inforef = ZHETRF_RK( zhetrf_rk_obj.matrix_layout,zhetrf_rk_obj.uplo,
                            zhetrf_rk_obj.n,zhetrf_rk_obj.aref,zhetrf_rk_obj.lda,
                                        zhetrf_rk_obj.eref,zhetrf_rk_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zhetrf_rk_obj.info = LAPACKE_zhetrf_rk(zhetrf_rk_obj.matrix_layout,zhetrf_rk_obj.uplo,
                                   zhetrf_rk_obj.n,zhetrf_rk_obj.a,zhetrf_rk_obj.lda,
                                                  zhetrf_rk_obj.e,zhetrf_rk_obj.ipiv);

    if( zhetrf_rk_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrf_rk is wrong\n",
                    zhetrf_rk_obj.info );
    }
    if( zhetrf_rk_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrf_rk is wrong\n",
        zhetrf_rk_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zhetrf_rk_obj.n,zhetrf_rk_obj.ipiv,zhetrf_rk_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf_rk1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff  =  computeDiff_z( zhetrf_rk_obj.n*zhetrf_rk_obj.lda,zhetrf_rk_obj.a,
                                                          zhetrf_rk_obj.aref );
    diff +=  computeDiff_z( zhetrf_rk_obj.n,zhetrf_rk_obj.e,zhetrf_rk_obj.eref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
