#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hegv_2stage_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (b!=NULL)  free(b);\
if (bref!=NULL)  free(bref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hegv_2stage_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	void *hModule, *dModule;
	float diff_a;
	float diff_w;
	float diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A, *b;	
	char uplo;
	char jobz;
	lapack_int lda;
	lapack_int ldb;
	lapack_int itype;
	/*Output Parameter*/
	float* w;	
	lapack_complex_float *Aref, *bref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hegv_2stage_scomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n, lapack_int lda);
      ~hegv_2stage_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hegv_2stage_scomplex_parameters:: hegv_2stage_scomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	lda = lda_i;
	ldb = lda;
	uplo = uplo_i;
	jobz = jobz_i;
	itype = itype_i;
	
	//if (jobz == 'U')
		jobz = 'N'; //'V' is not supported
	if (itype == 4)
		itype = 1;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hegv_2stage scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d lda: %d \n", matrix_layout, jobz,uplo, n, lda);
	#endif

	bufsize_a = lda*n;
	
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(b == NULL) || (bref == NULL)){
		EXPECT_FALSE( true) << "hegv_2stage_float_parameters object: malloc error.";
		hegv_2stage_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand(b, bref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hegv_2stage_scomplex_parameters :: ~hegv_2stage_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hegv_2stage_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hegv_2stage_free();
  

}
/*  Test fixture class definition */
class chegv_2stage_test  : public  ::testing::Test {
public:
   hegv_2stage_scomplex_parameters  *chegv_2stage_obj;
   void SetUp();
   void TearDown () { delete chegv_2stage_obj; }
};

void chegv_2stage_test::SetUp(){

    /* LAPACKE chegv_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_chegv_2stage) (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,lapack_complex_float *A, 
										lapack_int lda, lapack_complex_float* b, lapack_int ldb, float* w);

    Fptr_NL_LAPACKE_chegv_2stage chegv_2stage;

    chegv_2stage_obj = new hegv_2stage_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].nrhs,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].n,
							eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chegv_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chegv_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chegv_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chegv_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chegv_2stage = (Fptr_NL_LAPACKE_chegv_2stage)dlsym(chegv_2stage_obj->hModule, "LAPACKE_chegv_2stage");
    ASSERT_TRUE(chegv_2stage != NULL) << "failed to get the Netlib LAPACKE_chegv_2stage symbol";
    

    chegv_2stage_obj->inforef = chegv_2stage( chegv_2stage_obj->matrix_layout, chegv_2stage_obj->itype, chegv_2stage_obj->jobz, chegv_2stage_obj->uplo, 
								chegv_2stage_obj->n, chegv_2stage_obj->Aref, chegv_2stage_obj->lda, chegv_2stage_obj->bref, chegv_2stage_obj->ldb, chegv_2stage_obj->w);

    /* Compute libflame's Lapacke o/p  */
    chegv_2stage_obj->info = LAPACKE_chegv_2stage( chegv_2stage_obj->matrix_layout, chegv_2stage_obj->itype, chegv_2stage_obj->jobz, chegv_2stage_obj->uplo,
										chegv_2stage_obj->n,chegv_2stage_obj->A, chegv_2stage_obj->lda, chegv_2stage_obj->b, chegv_2stage_obj->ldb, chegv_2stage_obj->w);

    if( chegv_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chegv_2stage is wrong\n", chegv_2stage_obj->info );
    }
    if( chegv_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chegv_2stage is wrong\n", 
        chegv_2stage_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chegv_2stage_obj->diff_a =  computeDiff_c( chegv_2stage_obj->bufsize_a, chegv_2stage_obj->A, chegv_2stage_obj->Aref );
	chegv_2stage_obj->diff_b =  computeDiff_c( chegv_2stage_obj->bufsize_a, chegv_2stage_obj->b, chegv_2stage_obj->bref );
	chegv_2stage_obj->diff_w =  computeDiff_s( chegv_2stage_obj->n, chegv_2stage_obj->w, chegv_2stage_obj->wref );
}

TEST_F(chegv_2stage_test, chegv_2stage1) {
    EXPECT_NEAR(0.0, chegv_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chegv_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, chegv_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chegv_2stage_test, chegv_2stage2) {
    EXPECT_NEAR(0.0, chegv_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chegv_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, chegv_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chegv_2stage_test, chegv_2stage3) {
    EXPECT_NEAR(0.0, chegv_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chegv_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, chegv_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chegv_2stage_test, chegv_2stage4) {
    EXPECT_NEAR(0.0, chegv_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chegv_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, chegv_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hegv_2stage_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	void *hModule, *dModule;
	 double diff_a;
	double diff_w;
	double diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* A, *b;
	lapack_int itype;
	char uplo;
	char jobz;
	lapack_int lda;
	lapack_int ldb;
	/*Output Parameter*/
	double* w;	
	lapack_complex_double *Aref, *bref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hegv_2stage_dcomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n, lapack_int lda);
      ~hegv_2stage_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hegv_2stage_dcomplex_parameters:: hegv_2stage_dcomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	ldb = lda;
	uplo = uplo_i;
	jobz = jobz_i;
	itype = itype_i;
	
	//if (jobz == 'U')
		jobz = 'N'; //'V' is not supported
	if (itype == 4)
		itype = 1;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hegv_2stage scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d lda: %d \n", matrix_layout, jobz,uplo, n, lda);
	#endif

	bufsize_a = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(b == NULL) || (bref == NULL)){
		EXPECT_FALSE( true) << "hegv_2stage_double_parameters object: malloc error.";
		hegv_2stage_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(b, bref, bufsize_a);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hegv_2stage_dcomplex_parameters :: ~hegv_2stage_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hegv_2stage_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hegv_2stage_free();

}
/*  Test fixture class definition */ 
class zhegv_2stage_test  : public  ::testing::Test {
public:
   hegv_2stage_dcomplex_parameters  *zhegv_2stage_obj;
   void SetUp();
   void TearDown () { delete zhegv_2stage_obj; }
};

void zhegv_2stage_test::SetUp(){

    /* LAPACKE zhegv_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zhegv_2stage) (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,lapack_complex_double *A, 
										lapack_int lda, lapack_complex_double* b, lapack_int ldb, double* w );

    Fptr_NL_LAPACKE_zhegv_2stage zhegv_2stage;

    zhegv_2stage_obj = new hegv_2stage_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhegv_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhegv_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhegv_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhegv_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhegv_2stage = (Fptr_NL_LAPACKE_zhegv_2stage)dlsym(zhegv_2stage_obj->hModule, "LAPACKE_zhegv_2stage");
    ASSERT_TRUE(zhegv_2stage != NULL) << "failed to get the Netlib LAPACKE_zhegv_2stage symbol";
    

    zhegv_2stage_obj->inforef = zhegv_2stage( zhegv_2stage_obj->matrix_layout, zhegv_2stage_obj->itype, zhegv_2stage_obj->jobz, zhegv_2stage_obj->uplo,
								zhegv_2stage_obj->n, zhegv_2stage_obj->Aref, zhegv_2stage_obj->lda, zhegv_2stage_obj->bref, zhegv_2stage_obj->ldb, zhegv_2stage_obj->wref);

    /* Compute libflame's Lapacke o/p  */
    zhegv_2stage_obj->info = LAPACKE_zhegv_2stage( zhegv_2stage_obj->matrix_layout, zhegv_2stage_obj->itype, zhegv_2stage_obj->jobz, zhegv_2stage_obj->uplo,
										zhegv_2stage_obj->n,zhegv_2stage_obj->A, zhegv_2stage_obj->lda, zhegv_2stage_obj->b, zhegv_2stage_obj->ldb, zhegv_2stage_obj->w);

    if( zhegv_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhegv_2stage is wrong\n", zhegv_2stage_obj->info );
    }
    if( zhegv_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhegv_2stage is wrong\n", 
        zhegv_2stage_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhegv_2stage_obj->diff_a =  computeDiff_z( zhegv_2stage_obj->bufsize_a, zhegv_2stage_obj->A, zhegv_2stage_obj->Aref );
	zhegv_2stage_obj->diff_b =  computeDiff_z( zhegv_2stage_obj->bufsize_a, zhegv_2stage_obj->b, zhegv_2stage_obj->bref );
	zhegv_2stage_obj->diff_w =  computeDiff_d( zhegv_2stage_obj->n, zhegv_2stage_obj->w, zhegv_2stage_obj->wref );
}

TEST_F(zhegv_2stage_test, zhegv_2stage1) {
    EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhegv_2stage_test, zhegv_2stage2) {
    EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhegv_2stage_test, zhegv_2stage3) {
    EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhegv_2stage_test, zhegv_2stage4) {
    EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, zhegv_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}
