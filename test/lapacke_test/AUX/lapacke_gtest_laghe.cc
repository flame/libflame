#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define laghe_free() \
if (a!=NULL)    free(a); \
if (aref!=NULL) free(aref);\
if (d!=NULL)    free(d); \
if (dref!=NULL)    free(dref); \
if (iseed!=NULL)  free(iseed);\
if (iseedref!=NULL) free(iseedref); \
if( hModule != NULL) dlclose(hModule); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin scomplex_common_parameters  class definition */
class laghe_scomplex_parameters{

   public:
	int bufsize;
	int bufsize_iseed;
	int matrix_layout;
	void *hModule, *dModule;
	float diff;
   /*input parameters */	
	lapack_int n;
	lapack_int k, lda;
	lapack_complex_float * a;
	float *d;
	/*Output Parameter*/
	lapack_int* iseed, *iseedref;
	lapack_complex_float *aref;
	float *dref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      laghe_scomplex_parameters (lapack_int n, lapack_int k, lapack_int lda);
      ~laghe_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
laghe_scomplex_parameters:: laghe_scomplex_parameters ( lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	n = n_i;
	lda  = lda_i;
	k = k_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n laghe scomplex: n: %d , k: %d, lda:%d \n", n, k, lda);
	#endif
	
	bufsize = lda*n;
	bufsize_iseed = 4;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&a, &aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&iseed, &iseedref, bufsize_iseed);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, n);
	if ((a==NULL) || (aref==NULL) ||
		(iseed==NULL) || (iseedref==NULL) ||
		(d==NULL) || (dref==NULL)){
		EXPECT_FALSE( true) << "laghe_scomplex_parameters object: malloc error.";
		laghe_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_int_buffer_pair_rand(iseed, iseedref, bufsize_iseed);
	lapacke_gtest_init_float_buffer_pair_rand(d, dref, n);

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
laghe_scomplex_parameters :: ~laghe_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" laghe_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   laghe_free();

}
/*  Test fixture class definition */
class claghe_test  : public  ::testing::Test {
public:
   laghe_scomplex_parameters  *claghe_obj;
   void SetUp();
   void TearDown () { delete claghe_obj; }
};

void claghe_test::SetUp(){

    /* LAPACKE claghe prototype */
    typedef int (*Fptr_NL_LAPACKE_claghe) (int matrix_layout , lapack_int n , lapack_int k , const float * d , lapack_complex_float * a , lapack_int lda , lapack_int * iseed);

    Fptr_NL_LAPACKE_claghe claghe;

    claghe_obj = new laghe_scomplex_parameters ( eig_paramslist[idx].n,
                           eig_paramslist[idx].k,
						   eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    claghe_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    claghe_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(claghe_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(claghe_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    claghe = (Fptr_NL_LAPACKE_claghe)dlsym(claghe_obj->hModule, "LAPACKE_claghe");
    ASSERT_TRUE(claghe != NULL) << "failed to get the Netlib LAPACKE_claghe symbol";
    

    claghe_obj->inforef = claghe(claghe_obj->matrix_layout, claghe_obj->n,
								claghe_obj->k, (const float*)claghe_obj->dref, claghe_obj->aref,
								claghe_obj->lda, claghe_obj->iseedref);

    /* Compute libflame's Lapacke o/p  */
    claghe_obj->info = LAPACKE_claghe(claghe_obj->matrix_layout, claghe_obj->n,	claghe_obj->k, (const float*)claghe_obj->d, claghe_obj->a,claghe_obj->lda, claghe_obj->iseed);

    if( claghe_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_claghe is wrong\n", claghe_obj->info );
    }
    if( claghe_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_claghe is wrong\n", 
        claghe_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    claghe_obj->diff =  computeDiff_c( claghe_obj->bufsize, 
                claghe_obj->a, claghe_obj->aref );

}

TEST_F(claghe_test, claghe1) {
    EXPECT_NEAR(0.0, claghe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(claghe_test, claghe2) {
    EXPECT_NEAR(0.0, claghe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(claghe_test, claghe3) {
    EXPECT_NEAR(0.0, claghe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(claghe_test, claghe4) {
    EXPECT_NEAR(0.0, claghe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin dcomplex_common_parameters  class definition */
class laghe_dcomplex_parameters{

   public:
	int bufsize;
	int bufsize_iseed;
	int matrix_layout;
	void *hModule, *dModule;
	double diff;
   /*input parameters */	
	lapack_int n;
	lapack_int k, lda;
	lapack_complex_double * a;
	double *d;
	/*Output Parameter*/
	lapack_int* iseed, *iseedref;
	lapack_complex_double *aref;
	double *dref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      laghe_dcomplex_parameters (lapack_int n, lapack_int k, lapack_int lda);
      ~laghe_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
laghe_dcomplex_parameters:: laghe_dcomplex_parameters ( lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	n = n_i;
	lda  = lda_i;
	k = k_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n laghe dcomplex: n: %d , k: %d, lda:%d \n", n, k, lda);
	#endif
	
	bufsize = lda*n;
	bufsize_iseed = 4;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&a, &aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&iseed, &iseedref, bufsize_iseed);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, n);
	if ((a==NULL) || (aref==NULL) ||
		(iseed==NULL) || (iseedref==NULL) ||
		(d==NULL) || (dref==NULL)){
		EXPECT_FALSE( true) << "laghe_dcomplex_parameters object: malloc error.";
		laghe_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_int_buffer_pair_rand(iseed, iseedref, bufsize_iseed);
	lapacke_gtest_init_double_buffer_pair_rand(d, dref, n);

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
laghe_dcomplex_parameters :: ~laghe_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" laghe_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   laghe_free();

}
/*  Test fixture class definition */
class zlaghe_test  : public  ::testing::Test {
public:
   laghe_dcomplex_parameters  *zlaghe_obj;
   void SetUp();
   void TearDown () { delete zlaghe_obj; }
};

void zlaghe_test::SetUp(){

    /* LAPACKE zlaghe prototype */
    typedef int (*Fptr_NL_LAPACKE_zlaghe) (int matrix_layout , lapack_int n , lapack_int k , const double * d , lapack_complex_double * a , lapack_int lda , lapack_int * iseed);

    Fptr_NL_LAPACKE_zlaghe zlaghe;

    zlaghe_obj = new laghe_dcomplex_parameters ( eig_paramslist[idx].n,
                           eig_paramslist[idx].k,
						   eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    zlaghe_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlaghe_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlaghe_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlaghe_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlaghe = (Fptr_NL_LAPACKE_zlaghe)dlsym(zlaghe_obj->hModule, "LAPACKE_zlaghe");
    ASSERT_TRUE(zlaghe != NULL) << "failed to get the Netlib LAPACKE_zlaghe symbol";
    

    zlaghe_obj->inforef = zlaghe(zlaghe_obj->matrix_layout, zlaghe_obj->n,
								zlaghe_obj->k, (const double*)zlaghe_obj->dref, zlaghe_obj->aref,
								zlaghe_obj->lda, zlaghe_obj->iseedref);

    /* Compute libflame's Lapacke o/p  */
    zlaghe_obj->info = LAPACKE_zlaghe(zlaghe_obj->matrix_layout, zlaghe_obj->n,	zlaghe_obj->k, (const double*)zlaghe_obj->d, zlaghe_obj->a,zlaghe_obj->lda, zlaghe_obj->iseed);

    if( zlaghe_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlaghe is wrong\n", zlaghe_obj->info );
    }
    if( zlaghe_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlaghe is wrong\n", 
        zlaghe_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlaghe_obj->diff =  computeDiff_z( zlaghe_obj->bufsize, 
                zlaghe_obj->a, zlaghe_obj->aref );

}

TEST_F(zlaghe_test, zlaghe1) {
    EXPECT_NEAR(0.0, zlaghe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlaghe_test, zlaghe2) {
    EXPECT_NEAR(0.0, zlaghe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlaghe_test, zlaghe3) {
    EXPECT_NEAR(0.0, zlaghe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlaghe_test, zlaghe4) {
    EXPECT_NEAR(0.0, zlaghe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
