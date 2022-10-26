#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define hegst_free() \
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if (bp!=NULL)  free(bp);\
if (bpref!=NULL) free(bpref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if(lModule != NULL) dlclose(lModule);
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hegst_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_complex_float* ap;
	lapack_int itype;
	char uplo;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_float* bp;	
	lapack_complex_float *apref, *bpref;	
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_potrf, inforef_potrf;

   public:
      hegst_scomplex_parameters (int matrix_layout, lapack_int itype, char uplo, lapack_int n);
      ~hegst_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hegst_scomplex_parameters:: hegst_scomplex_parameters (int matrix_layout_i, lapack_int itype_i, char uplo_i, lapack_int n_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	itype = itype_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hegst scomplex: matrix_layout = %d, uplo:%c, n: %d lda: %d \n", matrix_layout, uplo, n, lda);
	#endif
	
	lda = ldb = n;
	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&bp, &bpref, bufsize);
	
	if ((ap==NULL) || (apref==NULL) || \
		(bp==NULL) || (bpref==NULL)) {
		EXPECT_FALSE( true) << "hegst_float_parameters object: malloc error.";
		hegst_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(ap, apref, bufsize, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(bp, bpref, bufsize, n, uplo);	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hegst_scomplex_parameters :: ~hegst_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hegst_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hegst_free();

}
/*  Test fixture class definition */
class chegst_test  : public  ::testing::Test {
public:
   hegst_scomplex_parameters  *chegst_obj;
   void SetUp();
   void TearDown () { delete chegst_obj; }
};

void chegst_test::SetUp(){

    /* LAPACKE cpotrf prototype */
	typedef int (*Fptr_NL_LAPACKE_cpotrf) (int matrix_layout, char uplo, lapack_int n,
											lapack_complex_float * bp , lapack_int ldb );
	Fptr_NL_LAPACKE_cpotrf cpotrf;
	/* LAPACKE chegst prototype */
    typedef int (*Fptr_NL_LAPACKE_chegst) (int matrix_layout, lapack_int itype, char uplo, lapack_int n,
										lapack_complex_float* ap, lapack_int lda, const lapack_complex_float* bp, lapack_int ldb);

    Fptr_NL_LAPACKE_chegst chegst;

    chegst_obj = new hegst_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].itype,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);
	
	/*potrf function call*/
	chegst_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chegst_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chegst_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chegst_obj->lModule != NULL) << "Netlib lapacke handle NULL";
	
    cpotrf = (Fptr_NL_LAPACKE_cpotrf)dlsym(chegst_obj->lModule, "LAPACKE_cpotrf");
    ASSERT_TRUE(cpotrf != NULL) << "failed to get the Netlib LAPACKE_chegst symbol";
	
	chegst_obj->inforef_potrf = cpotrf( chegst_obj->matrix_layout, chegst_obj->uplo,
								chegst_obj->n, chegst_obj->bpref, chegst_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    chegst_obj->info_potrf = LAPACKE_cpotrf( chegst_obj->matrix_layout, chegst_obj->uplo,
										chegst_obj->n,chegst_obj->bp,chegst_obj->ldb);

    if( chegst_obj->info_potrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chegst is wrong\n", chegst_obj->info_potrf );
    }
    if( chegst_obj->inforef_potrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chegst is wrong\n", 
        chegst_obj->inforef_potrf );
    }
	
	/*LAPACKE_chegst */

    chegst_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chegst_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chegst_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chegst_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chegst = (Fptr_NL_LAPACKE_chegst)dlsym(chegst_obj->hModule, "LAPACKE_chegst");
    ASSERT_TRUE(chegst != NULL) << "failed to get the Netlib LAPACKE_chegst symbol";
    

    chegst_obj->inforef = chegst( chegst_obj->matrix_layout, chegst_obj->itype, chegst_obj->uplo, chegst_obj->n, 
									chegst_obj->apref, chegst_obj->lda, chegst_obj->bpref, chegst_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    chegst_obj->info = LAPACKE_chegst( chegst_obj->matrix_layout, chegst_obj->itype, chegst_obj->uplo, chegst_obj->n,
										chegst_obj->ap, chegst_obj->lda, chegst_obj->bp, chegst_obj->ldb);

    if( chegst_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chegst is wrong\n", chegst_obj->info );
    }
    if( chegst_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chegst is wrong\n", 
        chegst_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chegst_obj->diff =  computeDiff_c( chegst_obj->bufsize, 
                chegst_obj->ap, chegst_obj->apref );

}

TEST_F(chegst_test, chegst1) {
    EXPECT_NEAR(0.0, chegst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chegst_test, chegst2) {
    EXPECT_NEAR(0.0, chegst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chegst_test, chegst3) {
    EXPECT_NEAR(0.0, chegst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chegst_test, chegst4) {
    EXPECT_NEAR(0.0, chegst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class hegst_dcomplex_parameters{
	
	 public:
	int bufsize;
	void *hModule, *dModule;
	void *lModule, *bModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* ap;
	lapack_int itype;
	char uplo;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_double* bp;	
	lapack_complex_double *apref, *bpref;	
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_potrf, inforef_potrf;

   public:
      hegst_dcomplex_parameters (int matrix_layout, lapack_int itype, char uplo, lapack_int n);
      ~hegst_dcomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hegst_dcomplex_parameters:: hegst_dcomplex_parameters (int matrix_layout_i, lapack_int itype_i, char uplo_i, lapack_int n_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	itype = itype_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hegst scomplex:  n: %d lda: %d \n", n, lda);
	#endif
	
	lda = ldb = n;
	bufsize = lda*n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&bp, &bpref, bufsize);	
	
	if ((ap==NULL) || (apref==NULL) || \
		(bp==NULL) || (bpref==NULL)) {
		EXPECT_FALSE( true) << "hegst_float_parameters object: malloc error.";
		hegst_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( ap, apref, bufsize, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( bp, bpref, bufsize, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hegst_dcomplex_parameters :: ~hegst_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hegst_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hegst_free();

}
/*  Test fixture class definition */
class zhegst_test  : public  ::testing::Test {
public:
   hegst_dcomplex_parameters  *zhegst_obj;
   void SetUp();
   void TearDown () { delete zhegst_obj; }
};

void zhegst_test::SetUp(){

    /* LAPACKE cpotrf prototype */
	typedef int (*Fptr_NL_LAPACKE_zpotrf) (int matrix_layout, char uplo, lapack_int n, 
											lapack_complex_double * bp , lapack_int ldb );
	Fptr_NL_LAPACKE_zpotrf zpotrf;
	/* LAPACKE zhegst prototype */
    typedef int (*Fptr_NL_LAPACKE_zhegst) (int matrix_layout, lapack_int itype, char uplo, lapack_int n,
										lapack_complex_double* ap, lapack_int lda, const lapack_complex_double* bp, lapack_int ldb);

    Fptr_NL_LAPACKE_zhegst zhegst;

    zhegst_obj = new hegst_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].itype,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);
	
	/*potrf function call*/
	zhegst_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhegst_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhegst_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhegst_obj->lModule != NULL) << "Netlib lapacke handle NULL";
	
    zpotrf = (Fptr_NL_LAPACKE_zpotrf)dlsym(zhegst_obj->lModule, "LAPACKE_zpotrf");
    ASSERT_TRUE(zpotrf != NULL) << "failed to get the Netlib LAPACKE_zhegst symbol";
	
	zhegst_obj->inforef_potrf = zpotrf( zhegst_obj->matrix_layout, zhegst_obj->uplo,
								zhegst_obj->n, zhegst_obj->bpref, zhegst_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    zhegst_obj->info_potrf = LAPACKE_zpotrf( zhegst_obj->matrix_layout, zhegst_obj->uplo,
										zhegst_obj->n,zhegst_obj->bp,zhegst_obj->ldb);

    if( zhegst_obj->info_potrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhegst is wrong\n", zhegst_obj->info_potrf );
    }
    if( zhegst_obj->inforef_potrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhegst is wrong\n", 
        zhegst_obj->inforef_potrf );
    }
	
	/*LAPACKE_zhegst */

    zhegst_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhegst_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhegst_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhegst_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhegst = (Fptr_NL_LAPACKE_zhegst)dlsym(zhegst_obj->hModule, "LAPACKE_zhegst");
    ASSERT_TRUE(zhegst != NULL) << "failed to get the Netlib LAPACKE_zhegst symbol";
    

    zhegst_obj->inforef = zhegst( zhegst_obj->matrix_layout, zhegst_obj->itype, zhegst_obj->uplo, zhegst_obj->n, 
									zhegst_obj->apref, zhegst_obj->lda, zhegst_obj->bpref, zhegst_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    zhegst_obj->info = LAPACKE_zhegst( zhegst_obj->matrix_layout, zhegst_obj->itype, zhegst_obj->uplo, zhegst_obj->n,
										zhegst_obj->ap, zhegst_obj->lda, zhegst_obj->bp, zhegst_obj->ldb);

    if( zhegst_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhegst is wrong\n", zhegst_obj->info );
    }
    if( zhegst_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhegst is wrong\n", 
        zhegst_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhegst_obj->diff =  computeDiff_z( zhegst_obj->bufsize, 
                zhegst_obj->ap, zhegst_obj->apref );

}

TEST_F(zhegst_test, zhegst1) {
    EXPECT_NEAR(0.0, zhegst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhegst_test, zhegst2) {
    EXPECT_NEAR(0.0, zhegst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhegst_test, zhegst3) {
    EXPECT_NEAR(0.0, zhegst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhegst_test, zhegst4) {
    EXPECT_NEAR(0.0, zhegst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

	