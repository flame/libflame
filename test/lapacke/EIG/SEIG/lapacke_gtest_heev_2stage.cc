#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define heev_2stage_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class heev_2stage_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff_a;
	float diff_w;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A;	
	char uplo;
	char jobz;
	lapack_int lda;
	/*Output Parameter*/
	float* w;	
	lapack_complex_float *Aref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      heev_2stage_scomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int lda);
      ~heev_2stage_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
heev_2stage_scomplex_parameters:: heev_2stage_scomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	uplo = uplo_i;
	jobz = jobz_i;
	
//	if (jobz == 'U')
		jobz = 'N'; // 'V' is not supported
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n heev_2stage scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d lda: %d \n", matrix_layout, jobz,uplo, n, lda);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL)){
		EXPECT_FALSE( true) << "heev_2stage_float_parameters object: malloc error.";
		heev_2stage_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	/*initialize output matrix by 0 */
	for(i=0;i<n;i++) {
		w[i] = 0;
		wref[i] = w[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
heev_2stage_scomplex_parameters :: ~heev_2stage_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heev_2stage_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heev_2stage_free();

}
/*  Test fixture class definition */
class cheev_2stage_test  : public  ::testing::Test {
public:
   heev_2stage_scomplex_parameters  *cheev_2stage_obj;
   void SetUp();
   void TearDown () { delete cheev_2stage_obj; }
};

void cheev_2stage_test::SetUp(){

    /* LAPACKE cheev_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_cheev_2stage) (int matrix_layout, char jobz, char uplo, lapack_int n,lapack_complex_float *A, 
											lapack_int lda, float* w);

    Fptr_NL_LAPACKE_cheev_2stage cheev_2stage;

    cheev_2stage_obj = new heev_2stage_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    cheev_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cheev_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cheev_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cheev_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cheev_2stage = (Fptr_NL_LAPACKE_cheev_2stage)dlsym(cheev_2stage_obj->hModule, "LAPACKE_cheev_2stage");
    ASSERT_TRUE(cheev_2stage != NULL) << "failed to get the Netlib LAPACKE_cheev_2stage symbol";
    

    cheev_2stage_obj->inforef = cheev_2stage( cheev_2stage_obj->matrix_layout, cheev_2stage_obj->jobz, cheev_2stage_obj->uplo,
								cheev_2stage_obj->n, cheev_2stage_obj->Aref, cheev_2stage_obj->lda, cheev_2stage_obj->wref);

    /* Compute libflame's Lapacke o/p  */
    cheev_2stage_obj->info = LAPACKE_cheev_2stage( cheev_2stage_obj->matrix_layout, cheev_2stage_obj->jobz, cheev_2stage_obj->uplo,
										cheev_2stage_obj->n,cheev_2stage_obj->A, cheev_2stage_obj->lda, cheev_2stage_obj->w);

    if( cheev_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cheev_2stage is wrong\n", cheev_2stage_obj->info );
    }
    if( cheev_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cheev_2stage is wrong\n", 
        cheev_2stage_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cheev_2stage_obj->diff_a =  computeDiff_c( cheev_2stage_obj->bufsize, cheev_2stage_obj->A, cheev_2stage_obj->Aref );
	cheev_2stage_obj->diff_w =  computeDiff_s( cheev_2stage_obj->n, cheev_2stage_obj->w, cheev_2stage_obj->wref );
}

TEST_F(cheev_2stage_test, cheev_2stage1) {
    EXPECT_NEAR(0.0, cheev_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheev_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheev_2stage_test, cheev_2stage2) {
    EXPECT_NEAR(0.0, cheev_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheev_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheev_2stage_test, cheev_2stage3) {
    EXPECT_NEAR(0.0, cheev_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheev_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheev_2stage_test, cheev_2stage4) {
    EXPECT_NEAR(0.0, cheev_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheev_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class heev_2stage_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	 double diff_a;
	double diff_w;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* A;	
	char uplo;
	char jobz;
	lapack_int lda;
	/*Output Parameter*/
	double* w;	
	lapack_complex_double *Aref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      heev_2stage_dcomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int lda);
      ~heev_2stage_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
heev_2stage_dcomplex_parameters:: heev_2stage_dcomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	uplo = uplo_i;
	jobz = jobz_i;
	
	//if (jobz == 'U')
		jobz = 'N'; // 'V' is not supported
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n heev_2stage scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d lda: %d \n", matrix_layout, jobz,uplo, n, lda);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL)){
		EXPECT_FALSE( true) << "heev_2stage_double_parameters object: malloc error.";
		heev_2stage_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	/*initialize output matrix by 0 */
	for(i=0;i<n;i++) {
		w[i] = 0;
		wref[i] = w[i];
	}

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
heev_2stage_dcomplex_parameters :: ~heev_2stage_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heev_2stage_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heev_2stage_free();

}
/*  Test fixture class definition */
class zheev_2stage_test  : public  ::testing::Test {
public:
   heev_2stage_dcomplex_parameters  *zheev_2stage_obj;
   void SetUp();
   void TearDown () { delete zheev_2stage_obj; }
};

void zheev_2stage_test::SetUp(){

    /* LAPACKE zheev_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zheev_2stage) (int matrix_layout, char jobz, char uplo, lapack_int n,lapack_complex_double *A, 
											lapack_int lda, double* w);

    Fptr_NL_LAPACKE_zheev_2stage zheev_2stage;

    zheev_2stage_obj = new heev_2stage_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zheev_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zheev_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zheev_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zheev_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zheev_2stage = (Fptr_NL_LAPACKE_zheev_2stage)dlsym(zheev_2stage_obj->hModule, "LAPACKE_zheev_2stage");
    ASSERT_TRUE(zheev_2stage != NULL) << "failed to get the Netlib LAPACKE_zheev_2stage symbol";
    

    zheev_2stage_obj->inforef = zheev_2stage( zheev_2stage_obj->matrix_layout, zheev_2stage_obj->jobz, zheev_2stage_obj->uplo,
								zheev_2stage_obj->n, zheev_2stage_obj->Aref, zheev_2stage_obj->lda, zheev_2stage_obj->wref);

    /* Compute libflame's Lapacke o/p  */
    zheev_2stage_obj->info = LAPACKE_zheev_2stage( zheev_2stage_obj->matrix_layout, zheev_2stage_obj->jobz, zheev_2stage_obj->uplo,
										zheev_2stage_obj->n,zheev_2stage_obj->A, zheev_2stage_obj->lda, zheev_2stage_obj->w);

    if( zheev_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zheev_2stage is wrong\n", zheev_2stage_obj->info );
    }
    if( zheev_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zheev_2stage is wrong\n", 
        zheev_2stage_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zheev_2stage_obj->diff_a =  computeDiff_z( zheev_2stage_obj->bufsize, zheev_2stage_obj->A, zheev_2stage_obj->Aref );
	zheev_2stage_obj->diff_w =  computeDiff_d( zheev_2stage_obj->n, zheev_2stage_obj->w, zheev_2stage_obj->wref );
}

TEST_F(zheev_2stage_test, zheev_2stage1) {
    EXPECT_NEAR(0.0, zheev_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheev_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheev_2stage_test, zheev_2stage2) {
    EXPECT_NEAR(0.0, zheev_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheev_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheev_2stage_test, zheev_2stage3) {
    EXPECT_NEAR(0.0, zheev_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheev_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheev_2stage_test, zheev_2stage4) {
    EXPECT_NEAR(0.0, zheev_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheev_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}
