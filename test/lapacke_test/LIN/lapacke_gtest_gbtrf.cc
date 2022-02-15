#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define gbtrf_free() \
       free (ab   ); \
       free (abref); \
       free (ipiv); \
       free (ipivref)

/* Begin gbtrf_double_parameters  class definition */
class gbtrf_double_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int m; // No of rows 
      lapack_int n; // No of Columns
      lapack_int kl; // No of subdiagonals within the band of A
      lapack_int ku; // No of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of the array ab

      /* Input/ Output parameters */
      double * ab,*abref; //The array ab contains the matrix A
      lapack_int *ipiv,*ipivref; // The pivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      gbtrf_double_parameters ( int matrix_layout_i,lapack_int m_i,
                                lapack_int n_i,lapack_int kl_i,
                                lapack_int ku_i,lapack_int ldab_i);
      ~gbtrf_double_parameters (); 
};  /* end of gbtrf_double_parameters  class definition */

/* Constructor gbtrf_double_parameters definition */
gbtrf_double_parameters:: gbtrf_double_parameters ( int matrix_layout_i,
                                       lapack_int m_i,lapack_int n_i,
                                       lapack_int kl_i,lapack_int ku_i,
                                       lapack_int ldab_i) {
    int arr_size;
    m = m_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    ldab = ldab_i;
    matrix_layout = matrix_layout_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &ab, &abref, (ldab*max(m,n)));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,max(1,min(m,n)));

    if( (ab==NULL) || (abref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL)){
       gbtrf_free();
       printf(" gbtrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ab,abref,(ldab*max(m,n)));

   } /* end of Constructor  */

gbtrf_double_parameters:: ~gbtrf_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbtrf_double_parameters object: destructor invoked. \n");
#endif
   gbtrf_free();
}

TEST(gbtrf,dgbtrf1) {

    /* LAPACKE dgbtrf prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbtrf) ( int matrix_layout ,lapack_int m ,
                                lapack_int n ,lapack_int kl ,lapack_int ku ,
                                double * ab ,lapack_int ldab ,
                                lapack_int * ipiv );

    Fptr_NL_LAPACKE_dgbtrf DGBTRF;
    double diff;
    int ipiv_diff;
    void *hModule,*dModule;

    gbtrf_double_parameters dgbtrf_obj(LAPACK_ROW_MAJOR,456,300,120,130,420);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    DGBTRF = (Fptr_NL_LAPACKE_dgbtrf)dlsym(hModule,"LAPACKE_dgbtrf");
    if (NULL == DGBTRF)
    {
        printf("Could not get the symbol LAPACKE_dgbtrf. Exiting...\n");
        dlclose(hModule);
        dlclose(dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    dgbtrf_obj.info     = LAPACKE_dgbtrf( LAPACK_COL_MAJOR,dgbtrf_obj.m,
                           dgbtrf_obj.n,dgbtrf_obj.kl,dgbtrf_obj.ku,
                           dgbtrf_obj.ab,dgbtrf_obj.ldab,dgbtrf_obj.ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgbtrf_obj.inforef = DGBTRF( LAPACK_COL_MAJOR,dgbtrf_obj.m,dgbtrf_obj.n,
                                 dgbtrf_obj.kl,dgbtrf_obj.ku,dgbtrf_obj.abref,
                                 dgbtrf_obj.ldab,dgbtrf_obj.ipivref);

    if( dgbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgbtrf \
                  is wrong\n",dgbtrf_obj.info );
    }
    if( dgbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgbtrf is wrong\n",
        dgbtrf_obj.inforef );
    }
    ipiv_diff = computeDiff_i( max(1,min(dgbtrf_obj.m,dgbtrf_obj.n)),
                                dgbtrf_obj.ipiv,dgbtrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: pivot computation in dgbtrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( (dgbtrf_obj.ldab * max(dgbtrf_obj.m,dgbtrf_obj.n) ),
                                             dgbtrf_obj.ab,dgbtrf_obj.abref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin gbtrf_float_parameters  class definition */
class gbtrf_float_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR 
      lapack_int m; // No of rows 
      lapack_int n; // No of Columns
      lapack_int kl; // No of subdiagonals within the band of A
      lapack_int ku; // No of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of the array ab

      /* Input/ Output parameters */
      float * ab,*abref; //The array ab contains the matrix A
      lapack_int * ipiv,*ipivref; // The pivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      gbtrf_float_parameters( int matrix_layout_i,lapack_int m_i,
                              lapack_int n_i,lapack_int kl_i,
                              lapack_int ku_i,lapack_int ldab_i);
      ~gbtrf_float_parameters (); 
    
};  /* end of gbtrf_float_parameters  class definition */

/* Constructor for gbtrf_float_parameters  */
gbtrf_float_parameters:: gbtrf_float_parameters ( int matrix_layout_i,
                                       lapack_int m_i,lapack_int n_i,
                                       lapack_int kl_i,lapack_int ku_i,
                                       lapack_int ldab_i) {
    int arr_size;
    m = m_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    ldab = ldab_i;
    matrix_layout = matrix_layout_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &ab, &abref, (ldab*max(m,n)));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,max(1,min(m,n)));

    if( (ab==NULL) || (abref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL)){
       gbtrf_free();
       printf(" gbtrf_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ab,abref,(ldab*max(m,n)));

   } /* end of Constructor  */

/* Destructor definition **/
gbtrf_float_parameters:: ~gbtrf_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbtrf_float_parameters object: destructor invoked. \n");
#endif
   gbtrf_free();
}

TEST(gbtrf,sgbtrf1) {

    /* LAPACKE sgbtrf prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbtrf) ( int matrix_layout ,lapack_int m ,
                                lapack_int n ,lapack_int kl ,lapack_int ku ,
                                float * ab ,lapack_int ldab ,
                                lapack_int * ipiv );

    Fptr_NL_LAPACKE_sgbtrf SGBTRF;
    float diff;
    int ipiv_diff;
    void *hModule,*dModule;

    gbtrf_float_parameters sgbtrf_obj(LAPACK_COL_MAJOR,456,300,120,130,420);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    SGBTRF = (Fptr_NL_LAPACKE_sgbtrf)dlsym(hModule,"LAPACKE_sgbtrf");
    if (NULL == SGBTRF)
    {
        printf("Could not get the symbol LAPACKE_sgbtrf. Exiting...\n");
        dlclose(hModule);
        dlclose(dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    sgbtrf_obj.info     = LAPACKE_sgbtrf( LAPACK_COL_MAJOR,sgbtrf_obj.m,
                           sgbtrf_obj.n,sgbtrf_obj.kl,sgbtrf_obj.ku,
                           sgbtrf_obj.ab,sgbtrf_obj.ldab,sgbtrf_obj.ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgbtrf_obj.inforef = SGBTRF( LAPACK_COL_MAJOR,sgbtrf_obj.m,sgbtrf_obj.n,
                                 sgbtrf_obj.kl,sgbtrf_obj.ku,sgbtrf_obj.abref,
                                 sgbtrf_obj.ldab,sgbtrf_obj.ipivref);

    if( sgbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgbtrf \
                  is wrong\n",sgbtrf_obj.info );
    }
    if( sgbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgbtrf is wrong\n",
        sgbtrf_obj.inforef );
    }
    ipiv_diff = computeDiff_i( max(1,min(sgbtrf_obj.m,sgbtrf_obj.n)),
                                sgbtrf_obj.ipiv,sgbtrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: pivot computation in sgbtrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( (sgbtrf_obj.ldab * max(sgbtrf_obj.m,sgbtrf_obj.n) ),
                                             sgbtrf_obj.ab,sgbtrf_obj.abref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin single precision complex_common_parameters  class definition */
class gbtrf_scomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //matrix storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR 
      lapack_int m; // No of rows 
      lapack_int n; // No of Columns
      lapack_int kl; // No of subdiagonals within the band of A
      lapack_int ku; // No of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of the array ab

      /* Input/ Output parameters */
      lapack_complex_float * ab,*abref; //The array ab contains the matrix A
      lapack_int *ipiv,*ipivref; // The pivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      gbtrf_scomplex_parameters ( int matrix_layout_i,lapack_int m_i,
                                  lapack_int n_i,lapack_int kl_i,
                                  lapack_int ku_i,lapack_int ldab_i);
      ~gbtrf_scomplex_parameters ();

};  /* end of gbtrf_scomplex_parameters class definition */

gbtrf_scomplex_parameters:: gbtrf_scomplex_parameters ( int matrix_layout_i,
                                       lapack_int m_i,lapack_int n_i,
                                       lapack_int kl_i,lapack_int ku_i,
                                       lapack_int ldab_i) {
    int arr_size;
    m = m_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    ldab = ldab_i;
    matrix_layout = matrix_layout_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ab, &abref, (ldab*max(m,n)));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,max(1,min(m,n)));

    if( (ab==NULL) || (abref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL)){
       gbtrf_free();
       printf(" gbtrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ab,abref,(ldab*max(m,n)));

   } /* end of Constructor  */

/* Destructor definition **/
gbtrf_scomplex_parameters:: ~gbtrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbtrf_scomplex_parameters object: destructor invoked. \n");
#endif
   gbtrf_free();
}

TEST(gbtrf,cgbtrf1) {

    /* LAPACKE cgbtrf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbtrf) ( int matrix_layout ,lapack_int m ,
                                lapack_int n ,lapack_int kl ,lapack_int ku ,
                                lapack_complex_float * ab ,lapack_int ldab ,
                                lapack_int * ipiv );

    Fptr_NL_LAPACKE_cgbtrf CGBTRF;
    float diff;
    int ipiv_diff;
    void *hModule,*dModule;

    gbtrf_scomplex_parameters cgbtrf_obj(LAPACK_ROW_MAJOR,456,300,120,130,420);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    CGBTRF = (Fptr_NL_LAPACKE_cgbtrf)dlsym(hModule,"LAPACKE_cgbtrf");
    if (NULL == CGBTRF)
    {
        printf("Could not get the symbol LAPACKE_cgbtrf. Exiting...\n");
        dlclose(hModule);
        dlclose(dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    cgbtrf_obj.info     = LAPACKE_cgbtrf( LAPACK_COL_MAJOR,cgbtrf_obj.m,
                           cgbtrf_obj.n,cgbtrf_obj.kl,cgbtrf_obj.ku,
                           cgbtrf_obj.ab,cgbtrf_obj.ldab,cgbtrf_obj.ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgbtrf_obj.inforef = CGBTRF( LAPACK_COL_MAJOR,cgbtrf_obj.m,cgbtrf_obj.n,
                                 cgbtrf_obj.kl,cgbtrf_obj.ku,cgbtrf_obj.abref,
                                 cgbtrf_obj.ldab,cgbtrf_obj.ipivref);

    if( cgbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgbtrf \
                  is wrong\n",cgbtrf_obj.info );
    }
    if( cgbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgbtrf is wrong\n",
        cgbtrf_obj.inforef );
    }
    ipiv_diff = computeDiff_i( max(1,min(cgbtrf_obj.m,cgbtrf_obj.n)),
                                cgbtrf_obj.ipiv,cgbtrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: pivot computation in cgbtrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( (cgbtrf_obj.ldab * max(cgbtrf_obj.m,cgbtrf_obj.n) ),
                                             cgbtrf_obj.ab,cgbtrf_obj.abref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin double precision complex_common_parameters  class definition */
class gbtrf_dcomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //matrix storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR 
      lapack_int m; // No of rows 
      lapack_int n; // No of Columns
      lapack_int kl; // No of subdiagonals within the band of A
      lapack_int ku; // No of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of the array ab

      /* Input/ Output parameters */
      lapack_complex_double * ab,*abref; //The array ab contains the matrix A
      lapack_int *ipiv,*ipivref; // The pivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      gbtrf_dcomplex_parameters( int matrix_layout_i,lapack_int m_i,
                                 lapack_int n_i,lapack_int kl_i,
                                 lapack_int ku_i,lapack_int ldab_i);
      ~gbtrf_dcomplex_parameters ();

};  /* end of gbtrf_dcomplex_parameters  class definition */

/* Destructor definition **/
gbtrf_dcomplex_parameters:: ~gbtrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbtrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   gbtrf_free();
}

/* Constructor definition **/
gbtrf_dcomplex_parameters:: gbtrf_dcomplex_parameters ( int matrix_layout_i,
                                       lapack_int m_i,lapack_int n_i,
                                       lapack_int kl_i,lapack_int ku_i,
                                       lapack_int ldab_i) {
    int arr_size;
    m = m_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    ldab = ldab_i;
    matrix_layout = matrix_layout_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ab, &abref, (ldab*max(m,n)));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (ab==NULL) || (abref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL)){
       gbtrf_free();
       printf(" gbtrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ab,abref,(ldab*max(m,n)));

   } /* end of Constructor  */

#if 0
TEST(gbtrf,dgbtrf1) {

    /* LAPACKE DGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbtrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_dgbtrf DGBTRF;  
    double diff;
    double_common_parameters dgbtrf_obj(512,1024);

    dgbtrf_obj.dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    dgbtrf_obj.hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == dgbtrf_obj.hModule) || (NULL == dgbtrf_obj.dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DGBTRF = (Fptr_NL_LAPACKE_dgbtrf)dlsym(dgbtrf_obj.hModule,"LAPACKE_dgbtrf");
    if (NULL == DGBTRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(dgbtrf_obj.hModule);
      dlclose(dgbtrf_obj.dModule);
      exit (-1);
    }
        /* Compute libflame's Lapacke o/p  */
    dgbtrf_obj.info     = LAPACKE_dgbtrf( LAPACK_COL_MAJOR,dgbtrf_obj.m,dgbtrf_obj.n,dgbtrf_obj.A,
                                                 dgbtrf_obj.lda,dgbtrf_obj.ipiv);
    
        /* Compute the reference o/p by invoking Netlib-Lapack's API */
        dgbtrf_obj.inforef = DGBTRF( LAPACK_COL_MAJOR,dgbtrf_obj.m,dgbtrf_obj.n,dgbtrf_obj.Aref,
                                              dgbtrf_obj.lda,dgbtrf_obj.ipivref);

    if( dgbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgbtrf is wrong\n",
                    dgbtrf_obj.info );
    }
    if( dgbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgbtrf is wrong\n",
        dgbtrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dgbtrf_obj.m*dgbtrf_obj.n,dgbtrf_obj.A,dgbtrf_obj.Aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}


TEST(gbtrf,sgbtrf1) {

    /* LAPACKE sgbtrf prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbtrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    float* a,lapack_int lda,lapack_int* ipiv );
                    
    Fptr_NL_LAPACKE_sgbtrf sgbtrf;  
    float diff;

    float_common_parameters sgbtrf_obj(900,900);

    sgbtrf_obj.dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    sgbtrf_obj.hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == sgbtrf_obj.hModule) || (NULL == sgbtrf_obj.dModule) )
    {
            printf("Load Library failed. Exiting ....\n");
            exit( -1 );
    }

    sgbtrf = (Fptr_NL_LAPACKE_sgbtrf)dlsym(sgbtrf_obj.hModule,"LAPACKE_sgbtrf");
    if (NULL == sgbtrf)
    {
          printf("Could not get the symbol. Exiting...\n");
          dlclose(sgbtrf_obj.hModule);
          dlclose(sgbtrf_obj.dModule);
          exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    sgbtrf_obj.info     = LAPACKE_sgbtrf( LAPACK_COL_MAJOR,sgbtrf_obj.m,sgbtrf_obj.n,sgbtrf_obj.A,
                                          sgbtrf_obj.lda,sgbtrf_obj.ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    sgbtrf_obj.inforef = sgbtrf( LAPACK_COL_MAJOR,sgbtrf_obj.m,sgbtrf_obj.n,sgbtrf_obj.Aref,
                                        sgbtrf_obj.lda,sgbtrf_obj.ipivref);

    if( sgbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgbtrf is wrong\n",
                    sgbtrf_obj.info );
    }
    if( sgbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgbtrf is wrong\n",
        sgbtrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( sgbtrf_obj.m*sgbtrf_obj.n,sgbtrf_obj.A,sgbtrf_obj.Aref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST(gbtrf,cgbtrf1) {

    /* LAPACKE sgbtrf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbtrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_cgbtrf cgbtrf;  
    float diff;

    scomplex_common_parameters cgbtrf_obj(900,900);

    cgbtrf_obj.dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    cgbtrf_obj.hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == cgbtrf_obj.hModule) || (NULL == cgbtrf_obj.dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    cgbtrf = (Fptr_NL_LAPACKE_cgbtrf)dlsym(cgbtrf_obj.hModule,"LAPACKE_cgbtrf");
    if (NULL == cgbtrf)
    {
        printf("Could not get the symbol. Exiting...\n");
        dlclose(cgbtrf_obj.hModule);
        dlclose(cgbtrf_obj.dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    cgbtrf_obj.info     = LAPACKE_cgbtrf( LAPACK_COL_MAJOR,cgbtrf_obj.m,cgbtrf_obj.n,cgbtrf_obj.A,
                                          cgbtrf_obj.lda,cgbtrf_obj.ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgbtrf_obj.inforef = cgbtrf( LAPACK_COL_MAJOR,cgbtrf_obj.m,cgbtrf_obj.n,cgbtrf_obj.Aref,
                                        cgbtrf_obj.lda,cgbtrf_obj.ipivref);

    if( cgbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgbtrf is wrong\n",
                    cgbtrf_obj.info );
    }
    if( cgbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgbtrf is wrong\n",
        cgbtrf_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( cgbtrf_obj.m*cgbtrf_obj.n,cgbtrf_obj.A,cgbtrf_obj.Aref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}


#endif

TEST(gbtrf,zgbtrf1) {

    /* LAPACKE zgbtrf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgbtrf) ( int matrix_layout ,lapack_int m ,
                                lapack_int n ,lapack_int kl ,lapack_int ku ,
                                lapack_complex_double * ab ,lapack_int ldab ,
                                lapack_int * ipiv );

    Fptr_NL_LAPACKE_zgbtrf ZGBTRF;
    double diff;
    int ipiv_diff;
    void *hModule,*dModule;

    gbtrf_dcomplex_parameters zgbtrf_obj(LAPACK_COL_MAJOR,256,300,120,130,400);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    ZGBTRF = (Fptr_NL_LAPACKE_zgbtrf)dlsym(hModule,"LAPACKE_zgbtrf");
    if (NULL == ZGBTRF)
    {
        printf("Could not get the symbol LAPACKE_zgbtrf. Exiting...\n");
        dlclose(hModule);
        dlclose(dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    zgbtrf_obj.info     = LAPACKE_zgbtrf( LAPACK_COL_MAJOR,zgbtrf_obj.m,
                              zgbtrf_obj.n,zgbtrf_obj.kl,zgbtrf_obj.ku,
                           zgbtrf_obj.ab,zgbtrf_obj.ldab,zgbtrf_obj.ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgbtrf_obj.inforef = ZGBTRF( LAPACK_COL_MAJOR,zgbtrf_obj.m,zgbtrf_obj.n,
                                 zgbtrf_obj.kl,zgbtrf_obj.ku,zgbtrf_obj.abref,
                                 zgbtrf_obj.ldab,zgbtrf_obj.ipivref);

    if( zgbtrf_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgbtrf is wrong\n",
                    zgbtrf_obj.info );
    }
    if( zgbtrf_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgbtrf is wrong\n",
        zgbtrf_obj.inforef );
    }
    ipiv_diff = computeDiff_i( max(1,min(zgbtrf_obj.m,zgbtrf_obj.n)),
                                zgbtrf_obj.ipiv,zgbtrf_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: pivot computation in zgbtrf1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( (zgbtrf_obj.ldab * max(zgbtrf_obj.m,zgbtrf_obj.n) ),
                                             zgbtrf_obj.ab,zgbtrf_obj.abref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}