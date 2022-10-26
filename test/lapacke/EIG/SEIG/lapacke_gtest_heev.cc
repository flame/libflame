#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define heev_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class heev_scomplex_parameters{

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
      heev_scomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int lda);
      ~heev_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
heev_scomplex_parameters:: heev_scomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	uplo = uplo_i;
	jobz = jobz_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n heev scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d lda: %d \n", matrix_layout, jobz,uplo, n, lda);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL)){
		EXPECT_FALSE( true) << "heev_float_parameters object: malloc error.";
		heev_free();
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
heev_scomplex_parameters :: ~heev_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heev_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heev_free();

}
/*  Test fixture class definition */
class cheev_test  : public  ::testing::Test {
public:
   heev_scomplex_parameters  *cheev_obj;
   void SetUp();
   void TearDown () { delete cheev_obj; }
};

void cheev_test::SetUp(){

    /* LAPACKE cheev prototype */
    typedef int (*Fptr_NL_LAPACKE_cheev) (int matrix_layout, char jobz, char uplo, lapack_int n,lapack_complex_float *A, 
											lapack_int lda, float* w);

    Fptr_NL_LAPACKE_cheev cheev;

    cheev_obj = new heev_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    cheev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cheev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cheev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cheev_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cheev = (Fptr_NL_LAPACKE_cheev)dlsym(cheev_obj->hModule, "LAPACKE_cheev");
    ASSERT_TRUE(cheev != NULL) << "failed to get the Netlib LAPACKE_cheev symbol";
    

    cheev_obj->inforef = cheev( cheev_obj->matrix_layout, cheev_obj->jobz, cheev_obj->uplo,
								cheev_obj->n, cheev_obj->Aref, cheev_obj->lda, cheev_obj->wref);

    /* Compute libflame's Lapacke o/p  */
    cheev_obj->info = LAPACKE_cheev( cheev_obj->matrix_layout, cheev_obj->jobz, cheev_obj->uplo,
										cheev_obj->n,cheev_obj->A, cheev_obj->lda, cheev_obj->w);

    if( cheev_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cheev is wrong\n", cheev_obj->info );
    }
    if( cheev_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cheev is wrong\n", 
        cheev_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cheev_obj->diff_a =  computeDiff_c( cheev_obj->bufsize, cheev_obj->A, cheev_obj->Aref );
	cheev_obj->diff_w =  computeDiff_s( cheev_obj->n, cheev_obj->w, cheev_obj->wref );
}

TEST_F(cheev_test, cheev1) {
    EXPECT_NEAR(0.0, cheev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheev_test, cheev2) {
    EXPECT_NEAR(0.0, cheev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheev_test, cheev3) {
    EXPECT_NEAR(0.0, cheev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheev_test, cheev4) {
    EXPECT_NEAR(0.0, cheev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class heev_dcomplex_parameters{

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
      heev_dcomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int lda);
      ~heev_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
heev_dcomplex_parameters:: heev_dcomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	uplo = uplo_i;
	jobz = jobz_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n heev scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d lda: %d \n", matrix_layout, jobz,uplo, n, lda);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL)){
		EXPECT_FALSE( true) << "heev_double_parameters object: malloc error.";
		heev_free();
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
heev_dcomplex_parameters :: ~heev_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heev_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heev_free();

}
/*  Test fixture class definition */
class zheev_test  : public  ::testing::Test {
public:
   heev_dcomplex_parameters  *zheev_obj;
   void SetUp();
   void TearDown () { delete zheev_obj; }
};

void zheev_test::SetUp(){

    /* LAPACKE zheev prototype */
    typedef int (*Fptr_NL_LAPACKE_zheev) (int matrix_layout, char jobz, char uplo, lapack_int n,lapack_complex_double *A, 
											lapack_int lda, double* w);

    Fptr_NL_LAPACKE_zheev zheev;

    zheev_obj = new heev_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zheev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zheev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zheev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zheev_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zheev = (Fptr_NL_LAPACKE_zheev)dlsym(zheev_obj->hModule, "LAPACKE_zheev");
    ASSERT_TRUE(zheev != NULL) << "failed to get the Netlib LAPACKE_zheev symbol";
    

    zheev_obj->inforef = zheev( zheev_obj->matrix_layout, zheev_obj->jobz, zheev_obj->uplo,
								zheev_obj->n, zheev_obj->Aref, zheev_obj->lda, zheev_obj->wref);

    /* Compute libflame's Lapacke o/p  */
    zheev_obj->info = LAPACKE_zheev( zheev_obj->matrix_layout, zheev_obj->jobz, zheev_obj->uplo,
										zheev_obj->n,zheev_obj->A, zheev_obj->lda, zheev_obj->w);

    if( zheev_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zheev is wrong\n", zheev_obj->info );
    }
    if( zheev_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zheev is wrong\n", 
        zheev_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zheev_obj->diff_a =  computeDiff_z( zheev_obj->bufsize, zheev_obj->A, zheev_obj->Aref );
	zheev_obj->diff_w =  computeDiff_d( zheev_obj->n, zheev_obj->w, zheev_obj->wref );
}

TEST_F(zheev_test, zheev1) {
    EXPECT_NEAR(0.0, zheev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheev_test, zheev2) {
    EXPECT_NEAR(0.0, zheev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheev_test, zheev3) {
    EXPECT_NEAR(0.0, zheev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheev_test, zheev4) {
    EXPECT_NEAR(0.0, zheev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}
