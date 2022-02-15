#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define getri_free() \
       free (a   ); \
       free (aref); \
       free (ipiv); \
       free (ipivref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin getri_float_parameters  class definition */
class getri_float_parameters{

   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Input/ Output parameters */
      float *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      getri_float_parameters ( int matrix_layout_i,
                    lapack_int n_i,lapack_int lda_i);
              
      ~getri_float_parameters (); 
};  /* end of getri_float_parameters  class definition */


/* Constructor getri_float_parameters definition */
getri_float_parameters:: getri_float_parameters ( int matrix_layout_i,
                                       lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n getri :  n: %d,  lda: %d \n", n, lda);
#endif

    /* Memory allocation of the buffers */
    //  getri API applicable for input square matrix 'A'
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       printf(" getri_float_parameters object: malloc error. Exiting...\n");
       getri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,n*n);

   } /* end of Constructor  */

getri_float_parameters:: ~getri_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" getri_float_parameters object: destructor invoked. \n");
#endif
   getri_free();
}

//  Test fixture class definition
class sgetri_test  : public  ::testing::Test {
public:
   getri_float_parameters  *sgetri_obj;
   void SetUp();  
   void TearDown () { delete sgetri_obj; }
};


//TEST(getri,sgetri1) {
void sgetri_test::SetUp(){

    /* LAPACKE SGETRI prototype */
    typedef int (*Fptr_NL_LAPACKE_sgetri) ( int matrix_layout, lapack_int n,
                          float *a, lapack_int lda, const lapack_int *ipiv);
    Fptr_NL_LAPACKE_sgetri SGETRI;

     /* LAPACKE SGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetrf) ( int matrix_layout,lapack_int m,
				  lapack_int n, float* a,lapack_int lda,lapack_int* ipiv );
    Fptr_NL_LAPACKE_dgetrf SGETRF;

    sgetri_obj = new  getri_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda );
    idx = Circular_Increment_Index(idx);

    sgetri_obj->dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    sgetri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    
    ASSERT_TRUE(sgetri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgetri_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGETRI = (Fptr_NL_LAPACKE_sgetri)dlsym(sgetri_obj->hModule,"LAPACKE_sgetri");
    ASSERT_TRUE(SGETRI != NULL) << "failed to get the Netlib LAPACKE_sgetri symbol";

    SGETRF = (Fptr_NL_LAPACKE_dgetrf)dlsym(sgetri_obj->hModule,"LAPACKE_sgetrf");
    ASSERT_TRUE(SGETRF != NULL) << "failed to get the Netlib LAPACKE_dgetrf symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    sgetri_obj->inforef = SGETRF( sgetri_obj->matrix_layout, sgetri_obj->n,
                                           sgetri_obj->n, sgetri_obj->aref,
                                     sgetri_obj->lda, sgetri_obj->ipivref);

    sgetri_obj->inforef = SGETRI( sgetri_obj->matrix_layout, sgetri_obj->n,
                                          sgetri_obj->aref,sgetri_obj->lda,
                                    (const lapack_int *)sgetri_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    sgetri_obj->inforef = LAPACKE_sgetrf( sgetri_obj->matrix_layout, 
						sgetri_obj->n, sgetri_obj->n, sgetri_obj->a,
                                 sgetri_obj->lda, sgetri_obj->ipiv);
                                    
    sgetri_obj->info = LAPACKE_sgetri( sgetri_obj->matrix_layout,
					sgetri_obj->n, sgetri_obj->a, sgetri_obj->lda,
                           (const lapack_int *) sgetri_obj->ipiv);

    if( sgetri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgetri is wrong\n",
                    sgetri_obj->info );
    }
    if( sgetri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgetri is wrong\n",
        sgetri_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgetri_obj->diff =  computeDiff_s( sgetri_obj->n*sgetri_obj->n,
								    sgetri_obj->a,sgetri_obj->aref );
}

TEST_F(sgetri_test, sgetri1) {
    EXPECT_NEAR(0.0, sgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgetri_test, sgetri2) {
    EXPECT_NEAR(0.0, sgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgetri_test, sgetri3) {
    EXPECT_NEAR(0.0, sgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgetri_test, sgetri4) {
    EXPECT_NEAR(0.0, sgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin getri_double_parameters  class definition */
class getri_double_parameters{

   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Input/ Output parameters */
      double *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      getri_double_parameters ( int matrix_layout_i,
                    lapack_int n_i,lapack_int lda_i);
              
      ~getri_double_parameters (); 
};  /* end of getri_double_parameters  class definition */


/* Constructor getri_double_parameters definition */
getri_double_parameters:: getri_double_parameters ( int matrix_layout_i,
                                       lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n getri :  n: %d,  lda: %d \n", n, lda);
#endif

    /* Memory allocation of the buffers */
    //  getri API applicable for input square matrix 'A'
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       printf(" getri_double_parameters object: malloc error. Exiting...\n");
       getri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,n*n);

   } /* end of Constructor  */

getri_double_parameters:: ~getri_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" getri_double_parameters object: destructor invoked. \n");
#endif
   getri_free();
}

//  Test fixture class definition
class dgetri_test  : public  ::testing::Test {
public:
   getri_double_parameters  *dgetri_obj;
   void SetUp();  
   void TearDown () { delete dgetri_obj; }
};


//TEST(getri,dgetri1) {
void dgetri_test::SetUp(){

    /* LAPACKE DGETRI prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetri) ( int matrix_layout,
											lapack_int n, 
											double *a,
											lapack_int lda, 
											const lapack_int *ipiv);
    Fptr_NL_LAPACKE_dgetri DGETRI;

     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetrf) ( int matrix_layout,
											lapack_int m,
											lapack_int n,
											double* a,
											lapack_int lda,
											lapack_int* ipiv );
    Fptr_NL_LAPACKE_dgetrf DGETRF;

    dgetri_obj = new  getri_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda );
    idx = Circular_Increment_Index(idx);

    dgetri_obj->dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    dgetri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    
    ASSERT_TRUE(dgetri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgetri_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGETRI = (Fptr_NL_LAPACKE_dgetri)dlsym(dgetri_obj->hModule,"LAPACKE_dgetri");
    ASSERT_TRUE(DGETRI != NULL) << "failed to get the Netlib LAPACKE_dgetri symbol";

    DGETRF = (Fptr_NL_LAPACKE_dgetrf)dlsym(dgetri_obj->hModule,"LAPACKE_dgetrf");
    ASSERT_TRUE(DGETRF != NULL) << "failed to get the Netlib LAPACKE_dgetrf symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    dgetri_obj->inforef = DGETRF( dgetri_obj->matrix_layout,
											dgetri_obj->n,
                                            dgetri_obj->n, 
											dgetri_obj->aref,
											dgetri_obj->lda, 
											dgetri_obj->ipivref);

    dgetri_obj->inforef = DGETRI( dgetri_obj->matrix_layout,
											dgetri_obj->n,
											dgetri_obj->aref,
											dgetri_obj->lda,
                         (const lapack_int *)dgetri_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    dgetri_obj->inforef = LAPACKE_dgetrf( dgetri_obj->matrix_layout,
											dgetri_obj->n,
                                            dgetri_obj->n, 
											dgetri_obj->a,
											dgetri_obj->lda, 
											dgetri_obj->ipiv);
                                    
    dgetri_obj->info = LAPACKE_dgetri( dgetri_obj->matrix_layout,
											dgetri_obj->n,
											dgetri_obj->a,
											dgetri_obj->lda,
                                  (const lapack_int *) dgetri_obj->ipiv);

    if( dgetri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgetri is wrong\n",
                    dgetri_obj->info );
    }
    if( dgetri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgetri is wrong\n",
        dgetri_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgetri_obj->diff =  computeDiff_d( dgetri_obj->n*dgetri_obj->n,dgetri_obj->a,dgetri_obj->aref );
}

TEST_F(dgetri_test, dgetri1) {
    EXPECT_NEAR(0.0, dgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetri_test, dgetri2) {
    EXPECT_NEAR(0.0, dgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetri_test, dgetri3) {
    EXPECT_NEAR(0.0, dgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetri_test, dgetri4) {
    EXPECT_NEAR(0.0, dgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin getri_lapack_complex_float_parameters  class definition */
class getri_lapack_complex_float_parameters{

   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Input/ Output parameters */
      lapack_complex_float *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      getri_lapack_complex_float_parameters ( int matrix_layout_i,
                    lapack_int n_i,lapack_int lda_i);
              
      ~getri_lapack_complex_float_parameters (); 
};  /* end of getri_lapack_complex_float_parameters  class definition */


/* Constructor getri_lapack_complex_float_parameters definition */
getri_lapack_complex_float_parameters:: getri_lapack_complex_float_parameters ( int matrix_layout_i,
                                       lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n getri :  n: %d,  lda: %d \n", n, lda);
#endif

    /* Memory allocation of the buffers */
    //  getri API applicable for input square matrix 'A'
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       printf(" getri_lapack_complex_float_parameters object: malloc error. Exiting...\n");
       getri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,n*n);

   } /* end of Constructor  */

getri_lapack_complex_float_parameters:: ~getri_lapack_complex_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" getri_lapack_complex_float_parameters object: destructor invoked. \n");
#endif
   getri_free();
}

//  Test fixture class definition
class cgetri_test  : public  ::testing::Test {
public:
   getri_lapack_complex_float_parameters  *cgetri_obj;
   void SetUp();  
   void TearDown () { delete cgetri_obj; }
};


//TEST(getri,cgetri1) {
void cgetri_test::SetUp(){

    /* LAPACKE CGETRI prototype */
    typedef int (*Fptr_NL_LAPACKE_cgetri) ( int matrix_layout,
											lapack_int n, 
											lapack_complex_float *a,
											lapack_int lda, 
											const lapack_int *ipiv);
    Fptr_NL_LAPACKE_cgetri CGETRI;

     /* LAPACKE CGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgetrf) ( int matrix_layout,
											lapack_int m,
											lapack_int n,
											lapack_complex_float* a,
											lapack_int lda,
											lapack_int* ipiv );
    Fptr_NL_LAPACKE_cgetrf CGETRF;

    cgetri_obj = new  getri_lapack_complex_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda );
    idx = Circular_Increment_Index(idx);

    cgetri_obj->dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    cgetri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    
    ASSERT_TRUE(cgetri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgetri_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGETRI = (Fptr_NL_LAPACKE_cgetri)dlsym(cgetri_obj->hModule,"LAPACKE_cgetri");
    ASSERT_TRUE(CGETRI != NULL) << "failed to get the Netlib LAPACKE_cgetri symbol";

    CGETRF = (Fptr_NL_LAPACKE_cgetrf)dlsym(cgetri_obj->hModule,"LAPACKE_cgetrf");
    ASSERT_TRUE(CGETRF != NULL) << "failed to get the Netlib LAPACKE_cgetrf symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    cgetri_obj->inforef = CGETRF( cgetri_obj->matrix_layout,
											cgetri_obj->n,
                                            cgetri_obj->n, 
											cgetri_obj->aref,
											cgetri_obj->lda, 
											cgetri_obj->ipivref);

    cgetri_obj->inforef = CGETRI( cgetri_obj->matrix_layout,
												cgetri_obj->n,
												cgetri_obj->aref,
												cgetri_obj->lda,
                                    (const lapack_int *)cgetri_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    cgetri_obj->inforef = LAPACKE_cgetrf( 	cgetri_obj->matrix_layout,
											cgetri_obj->n,
                                            cgetri_obj->n, 
											cgetri_obj->a,
											cgetri_obj->lda, 
											cgetri_obj->ipiv);
                                    
    cgetri_obj->info = LAPACKE_cgetri( 	cgetri_obj->matrix_layout, 
										cgetri_obj->n,
										cgetri_obj->a, 
										cgetri_obj->lda,
										(const lapack_int *) cgetri_obj->ipiv);

    if( cgetri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgetri is wrong\n",
                    cgetri_obj->info );
    }
    if( cgetri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgetri is wrong\n",
        cgetri_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgetri_obj->diff =  computeDiff_c( cgetri_obj->n*cgetri_obj->n,cgetri_obj->a,cgetri_obj->aref );
}

TEST_F(cgetri_test, cgetri1) {
    EXPECT_NEAR(0.0, cgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgetri_test, cgetri2) {
    EXPECT_NEAR(0.0, cgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgetri_test, cgetri3) {
    EXPECT_NEAR(0.0, cgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgetri_test, cgetri4) {
    EXPECT_NEAR(0.0, cgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin getri_lapack_complex_double_parameters  class definition */
class getri_lapack_complex_double_parameters{

   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Input/ Output parameters */
      lapack_complex_double *a,*aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info,inforef;

   public: 
      getri_lapack_complex_double_parameters ( int matrix_layout_i,
                    lapack_int n_i,lapack_int lda_i);
              
      ~getri_lapack_complex_double_parameters (); 
};  /* end of getri_lapack_complex_double_parameters  class definition */


/* Constructor getri_lapack_complex_double_parameters definition */
getri_lapack_complex_double_parameters:: getri_lapack_complex_double_parameters ( int matrix_layout_i,
                                       lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n getri :  n: %d,  lda: %d \n", n, lda);
#endif

    /* Memory allocation of the buffers */
    //  getri API applicable for input square matrix 'A'
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       printf(" getri_lapack_complex_double_parameters object: malloc error. Exiting...\n");
       getri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,n*n);

   } /* end of Constructor  */

getri_lapack_complex_double_parameters:: ~getri_lapack_complex_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" getri_lapack_complex_double_parameters object: destructor invoked. \n");
#endif
   getri_free();
}

//  Test fixture class definition
class zgetri_test  : public  ::testing::Test {
public:
   getri_lapack_complex_double_parameters  *zgetri_obj;
   void SetUp();  
   void TearDown () { delete zgetri_obj; }
};


//TEST(getri,zgetri1) {
void zgetri_test::SetUp(){

    /* LAPACKE ZGETRI prototype */
    typedef int (*Fptr_NL_LAPACKE_zgetri) ( int matrix_layout,
											lapack_int n, 
											lapack_complex_double *a,
											lapack_int lda, 
											const lapack_int *ipiv);
    Fptr_NL_LAPACKE_zgetri ZGETRI;

     /* LAPACKE ZGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zgetrf) ( int matrix_layout,
											lapack_int m,
											lapack_int n,
											lapack_complex_double* a,
											lapack_int lda,
											lapack_int* ipiv );
    Fptr_NL_LAPACKE_zgetrf ZGETRF;

    zgetri_obj = new  getri_lapack_complex_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda );
    idx = Circular_Increment_Index(idx);

    zgetri_obj->dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    zgetri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    
    ASSERT_TRUE(zgetri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgetri_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGETRI = (Fptr_NL_LAPACKE_zgetri)dlsym(zgetri_obj->hModule,"LAPACKE_zgetri");
    ASSERT_TRUE(ZGETRI != NULL) << "failed to get the Netlib LAPACKE_zgetri symbol";

    ZGETRF = (Fptr_NL_LAPACKE_zgetrf)dlsym(zgetri_obj->hModule,"LAPACKE_zgetrf");
    ASSERT_TRUE(ZGETRF != NULL) << "failed to get the Netlib LAPACKE_zgetrf symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    zgetri_obj->inforef = ZGETRF( 	zgetri_obj->matrix_layout,
									zgetri_obj->n,
									zgetri_obj->n,
									zgetri_obj->aref,
                                    zgetri_obj->lda,
									zgetri_obj->ipivref);

    zgetri_obj->inforef = ZGETRI( 	zgetri_obj->matrix_layout, 
									zgetri_obj->n,
									zgetri_obj->aref,
									zgetri_obj->lda,
                                    (const lapack_int *)zgetri_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    zgetri_obj->inforef = LAPACKE_zgetrf( 	zgetri_obj->matrix_layout,
											zgetri_obj->n,
                                            zgetri_obj->n, 
											zgetri_obj->a,
											zgetri_obj->lda,
											zgetri_obj->ipiv);
                                    
    zgetri_obj->info = LAPACKE_zgetri( 	zgetri_obj->matrix_layout,
										zgetri_obj->n,
										zgetri_obj->a,
										zgetri_obj->lda,
										(const lapack_int *) zgetri_obj->ipiv);

    if( zgetri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgetri is wrong\n",
                    zgetri_obj->info );
    }
    if( zgetri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgetri is wrong\n",
        zgetri_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgetri_obj->diff =  computeDiff_z( zgetri_obj->n*zgetri_obj->n,zgetri_obj->a,zgetri_obj->aref );
}

TEST_F(zgetri_test, zgetri1) {
    EXPECT_NEAR(0.0, zgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgetri_test, zgetri2) {
    EXPECT_NEAR(0.0, zgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgetri_test, zgetri3) {
    EXPECT_NEAR(0.0, zgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgetri_test, zgetri4) {
    EXPECT_NEAR(0.0, zgetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}