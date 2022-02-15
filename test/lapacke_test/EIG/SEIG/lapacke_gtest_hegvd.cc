#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hegvd_free() \
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
class hegvd_scomplex_parameters{

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
      hegvd_scomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n, lapack_int lda);
      ~hegvd_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hegvd_scomplex_parameters:: hegvd_scomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	ldb = lda;
	uplo = uplo_i;
	jobz = jobz_i;
	itype = itype_i;
	
	if (jobz == 'U')
		jobz = 'N';
	if (itype == 4)
		itype = 1;
	
	//#if LAPACKE_TEST_VERBOSE
	printf(" \n hegvd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d lda: %d \n", matrix_layout, jobz,uplo, n, lda);
//	#endif

	bufsize_a = lda*n;
	
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(b == NULL) || (bref == NULL)){
		EXPECT_FALSE( true) << "hegvd_float_parameters object: malloc error.";
		hegvd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand(b, bref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hegvd_scomplex_parameters :: ~hegvd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hegvd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hegvd_free();
  

}
/*  Test fixture class definition */
class chegvd_test  : public  ::testing::Test {
public:
   hegvd_scomplex_parameters  *chegvd_obj;
   void SetUp();
   void TearDown () { delete chegvd_obj; }
};

void chegvd_test::SetUp(){

    /* LAPACKE chegvd prototype */
    typedef int (*Fptr_NL_LAPACKE_chegvd) (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,lapack_complex_float *A, 
										lapack_int lda, lapack_complex_float* b, lapack_int ldb, float* w);

    Fptr_NL_LAPACKE_chegvd chegvd;

    chegvd_obj = new hegvd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].nrhs,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].n,
							eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chegvd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chegvd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chegvd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chegvd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chegvd = (Fptr_NL_LAPACKE_chegvd)dlsym(chegvd_obj->hModule, "LAPACKE_chegvd");
    ASSERT_TRUE(chegvd != NULL) << "failed to get the Netlib LAPACKE_chegvd symbol";
    

    chegvd_obj->inforef = chegvd( chegvd_obj->matrix_layout, chegvd_obj->itype, chegvd_obj->jobz, chegvd_obj->uplo, 
								chegvd_obj->n, chegvd_obj->Aref, chegvd_obj->lda, chegvd_obj->bref, chegvd_obj->ldb, chegvd_obj->w);

    /* Compute libflame's Lapacke o/p  */
    chegvd_obj->info = LAPACKE_chegvd( chegvd_obj->matrix_layout, chegvd_obj->itype, chegvd_obj->jobz, chegvd_obj->uplo,
										chegvd_obj->n,chegvd_obj->A, chegvd_obj->lda, chegvd_obj->b, chegvd_obj->ldb, chegvd_obj->w);

    if( chegvd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chegvd is wrong\n", chegvd_obj->info );
    }
    if( chegvd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chegvd is wrong\n", 
        chegvd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chegvd_obj->diff_a =  computeDiff_c( chegvd_obj->bufsize_a, chegvd_obj->A, chegvd_obj->Aref );
	chegvd_obj->diff_b =  computeDiff_c( chegvd_obj->bufsize_a, chegvd_obj->b, chegvd_obj->bref );
	chegvd_obj->diff_w =  computeDiff_s( chegvd_obj->n, chegvd_obj->w, chegvd_obj->wref );
}

TEST_F(chegvd_test, chegvd1) {
    EXPECT_NEAR(0.0, chegvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chegvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, chegvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chegvd_test, chegvd2) {
    EXPECT_NEAR(0.0, chegvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chegvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, chegvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chegvd_test, chegvd3) {
    EXPECT_NEAR(0.0, chegvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chegvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, chegvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chegvd_test, chegvd4) {
    EXPECT_NEAR(0.0, chegvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chegvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, chegvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hegvd_dcomplex_parameters{

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
      hegvd_dcomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n, lapack_int lda);
      ~hegvd_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hegvd_dcomplex_parameters:: hegvd_dcomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	ldb = lda;
	uplo = uplo_i;
	jobz = jobz_i;
	itype = itype_i;
	
	if (jobz == 'U')
		jobz = 'N';
	if (itype == 4)
		itype = 1;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hegvd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d lda: %d \n", matrix_layout, jobz,uplo, n, lda);
	#endif

	bufsize_a = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(b == NULL) || (bref == NULL)){
		EXPECT_FALSE( true) << "hegvd_double_parameters object: malloc error.";
		hegvd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(b, bref, bufsize_a);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hegvd_dcomplex_parameters :: ~hegvd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hegvd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hegvd_free();

}
/*  Test fixture class definition */ 
class zhegvd_test  : public  ::testing::Test {
public:
   hegvd_dcomplex_parameters  *zhegvd_obj;
   void SetUp();
   void TearDown () { delete zhegvd_obj; }
};

void zhegvd_test::SetUp(){

    /* LAPACKE zhegvd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhegvd) (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,lapack_complex_double *A, 
										lapack_int lda, lapack_complex_double* b, lapack_int ldb, double* w );

    Fptr_NL_LAPACKE_zhegvd zhegvd;

    zhegvd_obj = new hegvd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhegvd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhegvd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhegvd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhegvd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhegvd = (Fptr_NL_LAPACKE_zhegvd)dlsym(zhegvd_obj->hModule, "LAPACKE_zhegvd");
    ASSERT_TRUE(zhegvd != NULL) << "failed to get the Netlib LAPACKE_zhegvd symbol";
    

    zhegvd_obj->inforef = zhegvd( zhegvd_obj->matrix_layout, zhegvd_obj->itype, zhegvd_obj->jobz, zhegvd_obj->uplo,
								zhegvd_obj->n, zhegvd_obj->Aref, zhegvd_obj->lda, zhegvd_obj->bref, zhegvd_obj->ldb, zhegvd_obj->wref);

    /* Compute libflame's Lapacke o/p  */
    zhegvd_obj->info = LAPACKE_zhegvd( zhegvd_obj->matrix_layout, zhegvd_obj->itype, zhegvd_obj->jobz, zhegvd_obj->uplo,
										zhegvd_obj->n,zhegvd_obj->A, zhegvd_obj->lda, zhegvd_obj->b, zhegvd_obj->ldb, zhegvd_obj->w);

    if( zhegvd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhegvd is wrong\n", zhegvd_obj->info );
    }
    if( zhegvd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhegvd is wrong\n", 
        zhegvd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhegvd_obj->diff_a =  computeDiff_z( zhegvd_obj->bufsize_a, zhegvd_obj->A, zhegvd_obj->Aref );
	zhegvd_obj->diff_b =  computeDiff_z( zhegvd_obj->bufsize_a, zhegvd_obj->b, zhegvd_obj->bref );
	zhegvd_obj->diff_w =  computeDiff_d( zhegvd_obj->n, zhegvd_obj->w, zhegvd_obj->wref );
}

TEST_F(zhegvd_test, zhegvd1) {
    EXPECT_NEAR(0.0, zhegvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, zhegvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhegvd_test, zhegvd2) {
    EXPECT_NEAR(0.0, zhegvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, zhegvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhegvd_test, zhegvd3) {
    EXPECT_NEAR(0.0, zhegvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, zhegvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhegvd_test, zhegvd4) {
    EXPECT_NEAR(0.0, zhegvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, zhegvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}
