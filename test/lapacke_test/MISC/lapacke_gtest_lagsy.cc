#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define lagsy_free() \
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
/* Begin float_common_parameters  class definition */
class lagsy_float_parameters{

   public:
	int bufsize;
	int bufsize_iseed;
	int matrix_layout;
	void *hModule, *dModule;
	float diff;
   /*input parameters */	
	lapack_int n;
	lapack_int k, lda;
	float * a, *d;
	/*Output Parameter*/
	lapack_int* iseed, *iseedref;
	float *aref, *dref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      lagsy_float_parameters (lapack_int n, lapack_int k, lapack_int lda);
      ~lagsy_float_parameters ();

};

/* Constructor definition  float_common_parameters */
lagsy_float_parameters:: lagsy_float_parameters ( lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	n = n_i;
	lda  = lda_i;
	k = k_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lagsy float: n: %d , k: %d, lda:%d \n", n, k, lda);
	#endif
	
	bufsize = lda*n;
	bufsize_iseed = 4;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&a, &aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&iseed, &iseedref, bufsize_iseed);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, n);
	if ((a==NULL) || (aref==NULL) ||
		(iseed==NULL) || (iseedref==NULL) ||
		(d==NULL) || (dref==NULL)){
		EXPECT_FALSE( true) << "lagsy_float_parameters object: malloc error.";
		lagsy_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_int_buffer_pair_rand(iseed, iseedref, bufsize_iseed);
	lapacke_gtest_init_float_buffer_pair_rand(d, dref, n);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lagsy_float_parameters :: ~lagsy_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lagsy_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lagsy_free();

}
/*  Test fixture class definition */
class slagsy_test  : public  ::testing::Test {
public:
   lagsy_float_parameters  *slagsy_obj;
   void SetUp();
   void TearDown () { delete slagsy_obj; }
};

void slagsy_test::SetUp(){

    /* LAPACKE slagsy prototype */
    typedef int (*Fptr_NL_LAPACKE_slagsy) (int matrix_layout , lapack_int n , lapack_int k , const float * d , float * a , lapack_int lda , lapack_int * iseed);

    Fptr_NL_LAPACKE_slagsy slagsy;

    slagsy_obj = new lagsy_float_parameters ( eig_paramslist[idx].n,
                           eig_paramslist[idx].k,
						   eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    slagsy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slagsy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slagsy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slagsy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slagsy = (Fptr_NL_LAPACKE_slagsy)dlsym(slagsy_obj->hModule, "LAPACKE_slagsy");
    ASSERT_TRUE(slagsy != NULL) << "failed to get the Netlib LAPACKE_slagsy symbol";
    

    slagsy_obj->inforef = slagsy(slagsy_obj->matrix_layout, slagsy_obj->n,
								slagsy_obj->k, (const float*)slagsy_obj->dref, slagsy_obj->aref,
								slagsy_obj->lda, slagsy_obj->iseedref);

    /* Compute libflame's Lapacke o/p  */
    slagsy_obj->info = LAPACKE_slagsy(slagsy_obj->matrix_layout, slagsy_obj->n,	slagsy_obj->k, (const float*)slagsy_obj->d, slagsy_obj->a,slagsy_obj->lda, slagsy_obj->iseed);

    if( slagsy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slagsy is wrong\n", slagsy_obj->info );
    }
    if( slagsy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slagsy is wrong\n", 
        slagsy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slagsy_obj->diff =  computeDiff_s( slagsy_obj->bufsize, 
                slagsy_obj->a, slagsy_obj->aref );

}

TEST_F(slagsy_test, slagsy1) {
    EXPECT_NEAR(0.0, slagsy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slagsy_test, slagsy2) {
    EXPECT_NEAR(0.0, slagsy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slagsy_test, slagsy3) {
    EXPECT_NEAR(0.0, slagsy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slagsy_test, slagsy4) {
    EXPECT_NEAR(0.0, slagsy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class lagsy_double_parameters{

   public:
	int bufsize;
	int bufsize_iseed;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int k, lda; 
	double * a, *d;
	/*Output Parameter*/
	lapack_int* iseed, *iseedref;
	double *aref, *dref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      lagsy_double_parameters (lapack_int n, lapack_int k, lapack_int lda);
      ~lagsy_double_parameters ();

};

/* Constructor definition  double_common_parameters */
lagsy_double_parameters:: lagsy_double_parameters ( lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	n = n_i;
	lda  = lda_i;
	k = k_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lagsy double: n: %d , k: %d, lda:%d \n", n, k, lda);
	#endif
	
	bufsize = lda*n;
	bufsize_iseed = 4;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&a, &aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&iseed, &iseedref, bufsize_iseed);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, n);
	if ((a==NULL) || (aref==NULL) ||\
		(iseed==NULL) || (iseedref==NULL) ||\
		(d==NULL) || (dref==NULL)){
		EXPECT_FALSE( true) << "lagsy_double_parameters object: malloc error.";
		lagsy_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_int_buffer_pair_rand(iseed, iseedref, bufsize_iseed);
	lapacke_gtest_init_double_buffer_pair_rand(d, dref, n);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lagsy_double_parameters :: ~lagsy_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lagsy_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lagsy_free();

}
/*  Test fixture class definition */
class dlagsy_test  : public  ::testing::Test {
public:
   lagsy_double_parameters  *dlagsy_obj;
   void SetUp();
   void TearDown () { delete dlagsy_obj; }
};

void dlagsy_test::SetUp(){

    /* LAPACKE dlagsy prototype */
    typedef int (*Fptr_NL_LAPACKE_dlagsy) (int matrix_layout , lapack_int n , lapack_int k , const double * d , double * a , lapack_int lda , lapack_int * iseed);

    Fptr_NL_LAPACKE_dlagsy dlagsy;

    dlagsy_obj = new lagsy_double_parameters ( eig_paramslist[idx].n,
                           eig_paramslist[idx].k,
						   eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    dlagsy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlagsy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlagsy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlagsy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlagsy = (Fptr_NL_LAPACKE_dlagsy)dlsym(dlagsy_obj->hModule, "LAPACKE_dlagsy");
    ASSERT_TRUE(dlagsy != NULL) << "failed to get the Netlib LAPACKE_dlagsy symbol";
    

    dlagsy_obj->inforef = dlagsy(dlagsy_obj->matrix_layout, dlagsy_obj->n,
								dlagsy_obj->k, (const double*)dlagsy_obj->dref, dlagsy_obj->aref,
								dlagsy_obj->lda, dlagsy_obj->iseedref);

    /* Compute libflame's Lapacke o/p  */
    dlagsy_obj->info = LAPACKE_dlagsy(dlagsy_obj->matrix_layout, dlagsy_obj->n,
	dlagsy_obj->k, (const double*)dlagsy_obj->d, dlagsy_obj->a,dlagsy_obj->lda, dlagsy_obj->iseed);

    if( dlagsy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlagsy is wrong\n", dlagsy_obj->info );
    }
    if( dlagsy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlagsy is wrong\n", 
        dlagsy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlagsy_obj->diff =  computeDiff_d( dlagsy_obj->bufsize, 
                dlagsy_obj->a, dlagsy_obj->aref );

}

TEST_F(dlagsy_test, dlagsy1) {
    EXPECT_NEAR(0.0, dlagsy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlagsy_test, dlagsy2) {
    EXPECT_NEAR(0.0, dlagsy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlagsy_test, dlagsy3) {
    EXPECT_NEAR(0.0, dlagsy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlagsy_test, dlagsy4) {
    EXPECT_NEAR(0.0, dlagsy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
